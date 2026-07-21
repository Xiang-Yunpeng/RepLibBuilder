import os
import sys
import argparse
import subprocess
import shutil
from collections import defaultdict
from multiprocessing import Pool

# Dynamically import TEDictionary from step 1
BASE_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
if BASE_DIR not in sys.path:
    sys.path.append(BASE_DIR)

from modules.rating_module.normalizer import TEDictionary

# ---------------------------------------------------------------------------
# Classification policy for clustering
# ---------------------------------------------------------------------------
# The 8 classes the dictionary recognizes that PARTICIPATE in clustering.
# Anything the dictionary normalizes to a class OUTSIDE this set is treated as
# non-TE / structural (Satellite, Simple_repeat, Low_complexity, snRNA, tRNA,
# rRNA, scRNA, ncRNA, RNA, Segmental, Other, ARTEFACT, ...) and is HELD OUT
# from clustering entirely: such sequences are passed through to the final
# library unchanged. This avoids chimeric clusters, preserves distinct
# functional families (e.g. tRNA isoacceptors that can exceed 80% identity),
# and sidesteps CD-HIT-EST's alignment-coverage model, which is ill-suited to
# tandem repeats. The whitelist (rather than a blacklist) is robust to the
# dictionary's open-ended Global pass-through rules.
TE_CLUSTER_CLASSES = frozenset({
    "DNA", "LTR", "SINE", "LINE", "RC", "PLE", "RETROPOSON", "UNKNOWN"
})

# Single bin name used when stratification is disabled.
_COMBINED_BIN = "ALL"


class CDHITWrapper:
    def __init__(self, input_fasta, output_dir, te_dict,
                 c=0.8, aL=0.8, aS=0.8, threads=4, stratify=True,
                 dedup_unknown=False, dedup_id=0.8, dedup_cov=0.8, dedup_aL=0.0):
        self.input_fasta = os.path.abspath(input_fasta)
        self.output_dir = os.path.abspath(output_dir)
        self.te_dict = te_dict
        self.stratify = bool(stratify)

        # 80/80/80 Rule parameters (always user-settable, in either mode)
        self.c = float(c)
        self.aL = float(aL)
        self.aS = float(aS)
        self.threads = int(threads)
        self.n = self._calculate_word_length(self.c)

        self.output_fasta = os.path.join(self.output_dir, "step2_clustered_final.fa")
        self.tmp_dir = os.path.join(self.output_dir, "tmp_clustering")
        self.holdout_fasta = os.path.join(self.tmp_dir, "holdout_nonTE.fa")

        # Dictionaries for tracking statistics
        self.stats_before = {'DB': defaultdict(int), 'Class': defaultdict(int)}
        self.stats_after  = {'DB': defaultdict(int), 'Class': defaultdict(int)}

        # Per-bin sequence counts (filled during partitioning). Used to skip
        # CD-HIT for bins with <2 sequences (clustering 1 sequence is a no-op).
        self.bin_counts = defaultdict(int)

        # Optional post-clustering Unknown de-duplication (cd-hit-est-2d):
        # drop fully-Unknown consensus that a classified consensus already
        # covers. OFF by default; conservative 80/80 thresholds.
        self.dedup_unknown = bool(dedup_unknown)
        self.dedup_id = float(dedup_id)
        self.dedup_cov = float(dedup_cov)
        self.dedup_aL = float(dedup_aL)
        self.n_dropped_unknown = 0
        if self.dedup_unknown:
            if not (0.75 <= self.dedup_id <= 1.0):
                raise ValueError(
                    "--dedup_unknown_id must be within [0.75, 1.0] for cd-hit-est-2d "
                    f"(word-length floor); got {self.dedup_id}."
                )
            if not (0.0 <= self.dedup_cov <= 1.0):
                raise ValueError(
                    f"--dedup_unknown_cov must be within [0.0, 1.0]; got {self.dedup_cov}."
                )
            if not (0.0 <= self.dedup_aL <= 1.0):
                raise ValueError(
                    f"--dedup_unknown_aL must be within [0.0, 1.0]; got {self.dedup_aL}."
                )

    def _calculate_word_length(self, c_threshold):
        """
        Dynamically calculates the optimal word length (-n) based on
        the identity threshold (-c) according to CD-HIT official guidelines.
        """
        if c_threshold >= 0.95: return 10
        elif c_threshold >= 0.90: return 8
        elif c_threshold >= 0.88: return 7
        elif c_threshold >= 0.85: return 6
        elif c_threshold >= 0.80: return 5
        else: return 4

    def _parse_header(self, header):
        """
        Extracts (software, normalized_class_UPPER) from the header.

        This returns the TRUE normalized class instead of force-collapsing
        everything into {LTR, LINE, SINE, DNA, UNKNOWN}.
        The caller (_bin_for) decides whether the class is one of the 8
        clusterable classes or a held-out non-TE class.
        """
        clean_header = header[1:].split()[0]

        if '#' in clean_header:
            id_part, raw_class = clean_header.split('#', 1)
        else:
            id_part, raw_class = clean_header, "Unknown"

        # Isolate Software (DB) from the standardized "Software_Tag_ID" form
        software = id_part.split('_', 1)[0]

        # Get standardized class from TEDictionary
        norm_class, _ = self.te_dict.normalize(software, raw_class)
        return software, norm_class.upper()

    def _bin_for(self, norm_class_upper):
        """
        Map a normalized class to its clustering bin, or None if the sequence
        is non-TE and should be held out.
          - non-TE class          -> None (held out)
          - stratified clustering  -> the class itself (one of the 8 bins)
          - non-stratified         -> a single combined bin
        """
        if norm_class_upper not in TE_CLUSTER_CLASSES:
            return None
        return norm_class_upper if self.stratify else _COMBINED_BIN

    @staticmethod
    def _write_wrapped(handle, header, seq, width=60):
        """Write one FASTA record with the sequence wrapped at `width` bp."""
        handle.write(f"{header}\n")
        for i in range(0, len(seq), width):
            handle.write(seq[i:i + width] + "\n")

    def _route_record(self, header, seq, file_handles, holdout_handle):
        """
        Record before-stats and write a single FASTA record to its destination:
        a per-bin temp file (clusterable) or the held-out file (non-TE).
        Returns 1 if the record was held out, else 0.
        """
        software, norm_class = self._parse_header(header)
        self.stats_before['DB'][software] += 1
        self.stats_before['Class'][norm_class] += 1

        bin_name = self._bin_for(norm_class)
        if bin_name is None:
            # Non-TE: held out (passed through unchanged), wrapped at 60 bp.
            self._write_wrapped(holdout_handle, header, seq)
            return 1

        if bin_name not in file_handles:
            tmp_path = os.path.join(self.tmp_dir, f"split_{bin_name}.fa")
            file_handles[bin_name] = open(tmp_path, 'w')
        self._write_wrapped(file_handles[bin_name], header, seq)
        self.bin_counts[bin_name] += 1
        return 0

    def _partition(self):
        """
        Streams the input FASTA, routing each record to a clusterable per-bin
        file (TE classes + Unknown) or the held-out file (non-TE). Records the
        'before' statistics. Returns (list_of_bin_names, has_holdout).
        """
        mode = "stratified per-class" if self.stratify else "single combined pool"
        print(f"  --> Partitioning sequences ({mode}); non-TE sequences are held out from clustering...")
        os.makedirs(self.tmp_dir, exist_ok=True)

        self.bin_counts = defaultdict(int)
        file_handles = {}
        holdout_handle = open(self.holdout_fasta, 'w')
        holdout_count = 0

        header = ""
        seq_buffer = []

        with open(self.input_fasta, 'r') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    if header:
                        holdout_count += self._route_record(
                            header, "".join(seq_buffer), file_handles, holdout_handle
                        )
                    header = line
                    seq_buffer = []
                else:
                    seq_buffer.append(line)
            if header:  # Last record
                holdout_count += self._route_record(
                    header, "".join(seq_buffer), file_handles, holdout_handle
                )

        for fh in file_handles.values():
            fh.close()
        holdout_handle.close()

        has_holdout = holdout_count > 0
        if not has_holdout and os.path.exists(self.holdout_fasta):
            os.remove(self.holdout_fasta)

        return list(file_handles.keys()), has_holdout

    def _run_single_cdhit(self, args):
        """Worker function for CD-HIT execution (one bin)."""
        bin_name, in_file, out_file, threads = args

        if not os.path.exists(in_file):
            return True, bin_name

        cmd = [
            "cd-hit-est", "-i", in_file, "-o", out_file,
            "-c", str(self.c), "-n", str(self.n),
            "-aL", str(self.aL), "-aS", str(self.aS),
            "-T", str(threads), "-M", "0"
        ]

        try:
            subprocess.check_call(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            return True, bin_name
        except subprocess.CalledProcessError:
            return False, bin_name

    def _run_cdhit(self, cluster_bins):
        """
        Executes cd-hit-est for each clusterable bin. Runs bins in parallel via
        a multiprocessing Pool when stratified (multiple bins); a single
        in-process job (all threads) when there is only one bin.

        Bins with fewer than 2 sequences are NOT sent to cd-hit-est (clustering
        a single sequence is a no-op): their split file is copied straight to
        the clustered file so the merge step picks it up uniformly. This makes
        0- and 1-sequence classes safe regardless of cd-hit's behaviour on tiny
        inputs, and avoids spawning a pointless subprocess.
        """
        if not cluster_bins:
            print("  --> No clusterable (TE) sequences found; nothing to cluster.")
            return

        to_cluster = [b for b in cluster_bins if self.bin_counts.get(b, 0) >= 2]
        passthrough = [b for b in cluster_bins if self.bin_counts.get(b, 0) < 2]

        # Single-sequence bins: copy split -> clustered (no clustering needed)
        for bin_name in passthrough:
            src = os.path.join(self.tmp_dir, f"split_{bin_name}.fa")
            dst = os.path.join(self.tmp_dir, f"clustered_{bin_name}.fa")
            if os.path.exists(src):
                shutil.copyfile(src, dst)
        if passthrough:
            print(f"  --> {len(passthrough)} bin(s) had a single sequence; "
                  f"passed through without clustering ({', '.join(sorted(passthrough))}).")

        num_jobs = len(to_cluster)
        if num_jobs == 0:
            print("  --> No bin has >= 2 sequences; nothing to cluster with CD-HIT.")
            return

        print(f"  --> Executing CD-HIT-EST on {num_jobs} bin(s) "
              f"(Parameters: -c {self.c} -n {self.n} -aL {self.aL} -aS {self.aS})")

        # Smart thread allocation: avoid CD-HIT internal threads clashing with
        # the Python worker processes.
        threads_per_job = max(1, self.threads // num_jobs) if self.threads > num_jobs else 1
        pool_size = max(1, min(num_jobs, self.threads))

        tasks = []
        for bin_name in to_cluster:
            in_file = os.path.join(self.tmp_dir, f"split_{bin_name}.fa")
            out_file = os.path.join(self.tmp_dir, f"clustered_{bin_name}.fa")
            tasks.append((bin_name, in_file, out_file, threads_per_job))

        if pool_size == 1:
            # Single bin (e.g. non-stratified mode): run in-process with all threads.
            results = [self._run_single_cdhit(t) for t in tasks]
        else:
            with Pool(processes=pool_size) as pool:
                results = pool.map(self._run_single_cdhit, tasks)

        for success, bin_name in results:
            if not success:
                print(f"[Error] CD-HIT failed for bin: {bin_name}")
                sys.exit(1)

    @staticmethod
    def _iter_fasta_records(path):
        """Yield (header_line, sequence_string) for each record in a FASTA file."""
        header = None
        seq_buf = []
        with open(path, 'r') as f:
            for line in f:
                line = line.rstrip("\n")
                if not line:
                    continue
                if line.startswith(">"):
                    if header is not None:
                        yield header, "".join(seq_buf)
                    header = line
                    seq_buf = []
                else:
                    seq_buf.append(line)
            if header is not None:
                yield header, "".join(seq_buf)

    def _run_cdhit_est_2d(self, db1, db2, out_file):
        """Single-pairwise cd-hit-est-2d: keep db2 (Unknown) records NOT covered
        by any db1 (classified) record.

        -c / -aS act on the SHORTER sequence (typically the Unknown query), so
        this drops Unknown consensus that are >= dedup_cov covered by and
        >= dedup_id identical to a classified consensus. -aL (coverage of the
        LONGER sequence) defaults to 0.0 = off, so by default only the Unknown's
        own coverage is constrained; raise dedup_aL for a stricter, more
        symmetric containment test. `out_file` = retained Unknown;
        `out_file`.clstr records the covering rep.
        """
        cmd = [
            "cd-hit-est-2d", "-i", db1, "-i2", db2, "-o", out_file,
            "-c", str(self.dedup_id),
            "-aS", str(self.dedup_cov), "-aL", str(self.dedup_aL),
            "-n", str(self._calculate_word_length(self.dedup_id)),
            "-T", str(self.threads), "-M", "0", "-d", "0"
        ]
        subprocess.check_call(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    @staticmethod
    def _parse_clstr_coverage(clstr_path, unknown_ids):
        """Best-effort parse of a cd-hit-est-2d .clstr into
        {dropped_unknown_id: (covering_rep_id, identity)}. Non-critical: the
        drop-set itself is derived from the kept output, not from here.
        """
        cover = {}
        if not os.path.exists(clstr_path):
            return cover
        rep = None
        members = []

        def flush():
            for mid, ident in members:
                if mid in unknown_ids:
                    cover[mid] = (rep if rep is not None else "NA", ident)

        with open(clstr_path, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith(">Cluster"):
                    flush()
                    rep, members = None, []
                    continue
                if not line or ">" not in line:
                    continue
                name = line.split(">", 1)[1].split("...", 1)[0].strip()
                if line.endswith("*"):
                    rep = name
                else:
                    ident = "NA"
                    if "%" in line:
                        ident = line.rsplit("/", 1)[-1].rstrip("%").strip()
                    members.append((name, ident))
            flush()
        return cover

    def _dedup_unknown_global(self, cluster_bins):
        """Optional post-clustering step (cd-hit-est-2d): drop fully-Unknown
        consensus that a classified consensus already covers.

        Re-splits the clustered bins by normalized class (same _parse_header
        used for binning) into a classified reference and an Unknown query set,
        runs cd-hit-est-2d, and returns the set of dropped Unknown record IDs.
        Products/logs go to <output_dir>/dedup/ (NOT tmp_dir, which run()
        deletes).
        """
        dedup_dir = os.path.join(self.output_dir, "dedup")
        os.makedirs(dedup_dir, exist_ok=True)
        classified_fa = os.path.join(dedup_dir, "classified.fa")
        unknown_fa = os.path.join(dedup_dir, "unknown.fa")
        kept_fa = os.path.join(dedup_dir, "unknown_kept.fa")

        # 1. Re-split clustered bins by class.
        seen_unknown = set()
        dup_ids = False
        n_classified = 0
        with open(classified_fa, 'w') as fc, open(unknown_fa, 'w') as fu:
            for bin_name in cluster_bins:
                clustered_file = os.path.join(self.tmp_dir, f"clustered_{bin_name}.fa")
                if not os.path.exists(clustered_file):
                    continue
                for header, seq in self._iter_fasta_records(clustered_file):
                    _, norm_class = self._parse_header(header)
                    rec_id = header[1:].split()[0]
                    if norm_class == "UNKNOWN":
                        if rec_id in seen_unknown:
                            dup_ids = True
                        seen_unknown.add(rec_id)
                        self._write_wrapped(fu, header, seq)
                    else:
                        self._write_wrapped(fc, header, seq)
                        n_classified += 1

        # 2. Guard: skip if either side is empty (prevents cd-hit-est-2d crash).
        if not seen_unknown or n_classified == 0:
            print(f"  --> Unknown-dedup: nothing to dedup "
                  f"({len(seen_unknown)} Unknown, {n_classified} classified consensus).")
            return frozenset()
        if dup_ids:
            print("  --> [Warning] Duplicate Unknown record IDs detected; "
                  "the ID-keyed drop-set may be imprecise.")

        # 3. cd-hit-est-2d: classified = reference (db1), Unknown = query (db2).
        try:
            self._run_cdhit_est_2d(classified_fa, unknown_fa, kept_fa)
        except Exception as e:
            print(f"  --> [Warning] Unknown-dedup step failed ({e}); "
                  "proceeding WITHOUT dedup (full clustered library retained).")
            return frozenset()

        # 4. Dropped = all Unknown IDs - retained Unknown IDs.
        kept_ids = {h[1:].split()[0] for h, _ in self._iter_fasta_records(kept_fa)}
        dropped_ids = seen_unknown - kept_ids

        # 5. Traceability log (best-effort; never fatal).
        try:
            cover = self._parse_clstr_coverage(kept_fa + ".clstr", seen_unknown)
        except Exception:
            cover = {}
        # Defensive: cd-hit-est-2d (v4.8.1) seeds clusters only from db1
        # (classified), so every dropped Unknown is covered by a classified rep.
        # Warn if a future version clusters db2 (Unknown) internally.
        n_unknown_rep = sum(1 for rid in dropped_ids
                            if str(cover.get(rid, ("", ""))[0]).endswith("#Unknown"))
        if n_unknown_rep:
            print(f"  --> [Warning] {n_unknown_rep} dropped Unknown(s) appear covered by "
                  "another Unknown, not a classified consensus; check dropped_unknown.tsv.")
        tsv = os.path.join(dedup_dir, "dropped_unknown.tsv")
        with open(tsv, 'w') as ft:
            ft.write("dropped_unknown_id\tcovering_classified_rep\tidentity\n")
            for rid in sorted(dropped_ids):
                rep, ident = cover.get(rid, ("NA", "NA"))
                ft.write(f"{rid}\t{rep}\t{ident}\n")

        print(f"  --> Unknown-dedup: dropped {len(dropped_ids):,} of {len(seen_unknown):,} "
              f"Unknown consensus (>= {self.dedup_id} id / {self.dedup_cov} cov). Log: {tsv}")
        return frozenset(dropped_ids)

    def _merge_and_after_stats(self, cluster_bins, has_holdout, dropped_ids=frozenset()):
        """Merges clustered bins + held-out sequences and calculates final stats.

        Records whose sequence ID is in ``dropped_ids`` (redundant Unknown
        consensus removed by the optional cd-hit-est-2d dedup step) are skipped
        entirely, along with their sequence lines, and excluded from stats.
        """
        print("  --> Merging clustered bins and held-out sequences; compiling statistics...")

        with open(self.output_fasta, 'w') as fout:
            # 1. Clustered bins
            for bin_name in cluster_bins:
                clustered_file = os.path.join(self.tmp_dir, f"clustered_{bin_name}.fa")
                if not os.path.exists(clustered_file):
                    continue
                skip = False
                with open(clustered_file, 'r') as fin:
                    for line in fin:
                        if line.startswith(">"):
                            rec_id = line.strip()[1:].split()[0]
                            skip = rec_id in dropped_ids
                            if skip:
                                continue
                            software, out_class = self._parse_header(line.strip())
                            self.stats_after['DB'][software] += 1
                            self.stats_after['Class'][out_class] += 1
                            fout.write(line)
                        elif not skip:
                            fout.write(line)

            # 2. Held-out non-TE sequences (passed through unchanged)
            if has_holdout and os.path.exists(self.holdout_fasta):
                with open(self.holdout_fasta, 'r') as fin:
                    for line in fin:
                        fout.write(line)
                        if line.startswith(">"):
                            software, out_class = self._parse_header(line.strip())
                            self.stats_after['DB'][software] += 1
                            self.stats_after['Class'][out_class] += 1

    def _print_report(self):
        """Prints a formatted comparison of Before and After counts."""
        width = 64
        mode = "Stratified (per-class)" if self.stratify else "Single combined pool"

        print("\n" + "=" * width)
        print(f"{'CD-HIT Clustering Summary':^{width}}")
        print("=" * width)
        print(f"Mode  : {mode}")
        print(f"Params: -c {self.c}  -n {self.n}  -aL {self.aL}  -aS {self.aS}  -T {self.threads}")
        if self.dedup_unknown:
            extra = f" / {self.dedup_aL} covL" if self.dedup_aL > 0 else ""
            print(f"Unknown-dedup: dropped {self.n_dropped_unknown:,} redundant Unknown "
                  f"consensus (>= {self.dedup_id} id / {self.dedup_cov} cov{extra} vs a classified consensus)")

        print("\n--- By Database Source ---")
        print(f"{'Source':<22} | {'Before':>10} | {'After':>10}")
        print("-" * width)
        all_dbs = sorted(set(self.stats_before['DB'].keys()) | set(self.stats_after['DB'].keys()))
        for db in all_dbs:
            before = self.stats_before['DB'].get(db, 0)
            after = self.stats_after['DB'].get(db, 0)
            print(f"{db:<22} | {before:>10,} | {after:>10,}")

        # Clustered TE classes (the 8)
        clustered = [(c, n) for c, n in self.stats_before['Class'].items()
                     if c in TE_CLUSTER_CLASSES]
        print("\n--- By TE Class (clustered) ---")
        print(f"{'Class':<22} | {'Before':>10} | {'After':>10}")
        print("-" * width)
        for cls, before in sorted(clustered, key=lambda x: x[1], reverse=True):
            after = self.stats_after['Class'].get(cls, 0)
            print(f"{cls:<22} | {before:>10,} | {after:>10,}")

        # Held-out non-TE / structural classes (before == after by definition)
        held = [(c, n) for c, n in self.stats_before['Class'].items()
                if c not in TE_CLUSTER_CLASSES]
        if held:
            print("\n--- Non-TE / structural (held out, passed through unchanged) ---")
            print(f"{'Class':<22} | {'Count':>10}")
            print("-" * width)
            for cls, before in sorted(held, key=lambda x: x[1], reverse=True):
                print(f"{cls:<22} | {before:>10,}")

        total_before = sum(self.stats_before['DB'].values())
        total_after = sum(self.stats_after['DB'].values())
        print("-" * width)
        print(f"{'TOTAL SEQUENCES':<22} | {total_before:>10,} | {total_after:>10,}")
        print("=" * width)

    def run(self):
        mode = "Stratified" if self.stratify else "Non-stratified"
        print(f"\n[Module: CDHITWrapper] Initializing Clustering ({mode})...")

        cluster_bins, has_holdout = self._partition()
        self._run_cdhit(cluster_bins)
        dropped_unknown_ids = self._dedup_unknown_global(cluster_bins) if self.dedup_unknown else frozenset()
        self.n_dropped_unknown = len(dropped_unknown_ids)
        self._merge_and_after_stats(cluster_bins, has_holdout, dropped_unknown_ids)
        self._print_report()

        if os.path.exists(self.tmp_dir):
            shutil.rmtree(self.tmp_dir)

        print(f"\n[Success] Final clustered library saved to: {self.output_fasta}\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Stratified CD-HIT clustering for TE libraries.")
    parser.add_argument("-i", "--input", required=True, help="Input filtered FASTA")
    parser.add_argument("-o", "--out_dir", required=True, help="Output directory for clustered results")
    parser.add_argument("-m", "--mapping", required=True, help="Path to te_mapping_rules.tsv")
    parser.add_argument("-c", type=float, default=0.8, help="Sequence identity threshold (default: 0.8)")
    parser.add_argument("-aL", type=float, default=0.8, help="Alignment coverage for the longer sequence (default: 0.8)")
    parser.add_argument("-aS", type=float, default=0.8, help="Alignment coverage for the shorter sequence (default: 0.8)")
    parser.add_argument("-t", "--threads", type=int, default=8, help="Number of threads (default: 8)")
    parser.add_argument("--no_stratify", action="store_true",
                        help="Disable stratified (per-class) clustering; cluster all TE sequences "
                             "in a single CD-HIT run. Non-TE sequences are held out in either mode.")
    parser.add_argument("--dedup_unknown", action="store_true",
                        help="Opt-in: after clustering, drop fully-Unknown consensus that a classified "
                             "consensus covers (cd-hit-est-2d, single pairwise). Off by default.")
    parser.add_argument("--dedup_unknown_id", type=float, default=0.8,
                        help="Identity threshold (-c) for Unknown-dedup, on the shorter/Unknown seq (default: 0.8).")
    parser.add_argument("--dedup_unknown_cov", type=float, default=0.8,
                        help="Coverage threshold (-aS) on the Unknown seq for Unknown-dedup (default: 0.8).")
    parser.add_argument("--dedup_unknown_aL", type=float, default=0.0,
                        help="Coverage threshold (-aL) on the LONGER seq for Unknown-dedup "
                             "(default: 0.0 = off; raise for stricter containment).")

    args = parser.parse_args()

    te_dict = TEDictionary(args.mapping)
    cdhit_module = CDHITWrapper(
        input_fasta=args.input,
        output_dir=args.out_dir,
        te_dict=te_dict,
        c=args.c, aL=args.aL, aS=args.aS, threads=args.threads,
        stratify=not args.no_stratify,
        dedup_unknown=args.dedup_unknown,
        dedup_id=args.dedup_unknown_id,
        dedup_cov=args.dedup_unknown_cov,
        dedup_aL=args.dedup_unknown_aL
    )
    cdhit_module.run()
