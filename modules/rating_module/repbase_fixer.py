import os
import re
import sys
from collections import Counter
from utils.common import check_file, ensure_dir
from .normalizer import format_te_label


class RepbaseFixer:
    def __init__(self, input_fasta, output_dir, te_dict):
        """
        Args:
            input_fasta (str): Path to raw Repbase .ref/.fa file.
            output_dir (str): Output directory (e.g. 01.public_db/Repbase/).
            te_dict (TEDictionary): Pre-loaded TE dictionary instance for label
                                    normalization. Must contain Repbase rules.
        """
        self.input_fasta = check_file(input_fasta, "Repbase raw input")
        self.output_dir = ensure_dir(output_dir)
        self.te_dict = te_dict

        self.output_fasta = os.path.join(self.output_dir, "Repbase_custom.fixed.unique.fa")
        self.dropped_log = os.path.join(self.output_dir, "dropped_sequences.log")
        self.unmatched_log = os.path.join(self.output_dir, "unmatched_labels.log")

        # Pre-compiled patterns
        self.re_simple_repeat = re.compile(r'^>(.+)@([0-9]+)$')

    def _parse_header(self, header_line):
        """
        Parse a raw Repbase header line into (core_id, raw_label).

        Repbase native format is TAB-separated:
            >SeqID<TAB>Classification<TAB>Species [<TAB>...]

        Special cases:
          1. Simple repeat shorthand:  >NAME@123   (no TAB, label = Simple_repeat)
          2. Already-fixed input:      >Repbase_xxx#Class/Superfamily
          3. Whitespace-only fallback: >NAME Classification Species
                                       (older Repbase dumps using spaces)
          4. Bare ID, no fields:       >NAME
        """
        # Strip leading '>'
        body = header_line[1:].strip()

        # Defensive: drop pre-existing Repbase_ prefix
        if body.startswith("Repbase_"):
            body = body[len("Repbase_"):]

        # Case 1: Simple repeat shorthand (NAME@count)
        m_sr = self.re_simple_repeat.match(">" + body)
        if m_sr:
            core_id = f"{m_sr.group(1)}_{m_sr.group(2)}"
            return core_id, "Simple_repeat"

        # Case 2: Already-fixed header carrying #Class/Subclass
        if "#" in body and "\t" not in body:
            core_id, raw_label = body.split("#", 1)
            # Strip any trailing whitespace-padded species annotation
            raw_label = raw_label.split()[0] if raw_label else ""
            return core_id, raw_label

        # Case 3 (primary path): TAB-separated Repbase native format
        if "\t" in body:
            fields = body.split("\t")
            core_id = fields[0].strip()
            raw_label = fields[1].strip() if len(fields) > 1 else ""
            return core_id, raw_label

        # Case 4: Whitespace-separated legacy dump
        # Heuristic: split on runs of whitespace, take first token as ID,
        # second token (if any) as label. Only use when no TAB is present.
        if re.search(r'\s', body):
            parts = re.split(r'\s+', body, maxsplit=2)
            core_id = parts[0]
            raw_label = parts[1] if len(parts) > 1 else ""
            return core_id, raw_label

        # Case 5: Bare ID, no classification at all
        return body, ""

    def _build_new_header(self, core_id, norm_class, norm_superfam):
        """Construct the standardized RepeatMasker-compatible header."""
        clean_id = f"Repbase_{core_id}"
        # Shared single-source-of-truth formatter: produces "Class/Superfamily",
        # "Class/Unknown", or "Unknown" exactly as the previous inline logic did.
        label = format_te_label(norm_class, norm_superfam)
        return f">{clean_id}#{label}"

    def _apply_dedup(self, header, id_counts):
        """Append _dupN suffix to the ID portion if collisions occur."""
        # header format: >Repbase_xxx#Class/Subclass  OR  >Repbase_xxx
        if "#" in header:
            head, tail = header.split("#", 1)
            tail = "#" + tail
        else:
            head, tail = header, ""

        count = id_counts[head] + 1
        id_counts[head] = count

        if count > 1:
            head = f"{head}_dup{count}"
        return head + tail

    def run(self):
        print(f"\n[Module: RepbaseFixer] Processing")
        print(f"Input:  {self.input_fasta}")
        print(f"Output: {self.output_fasta}")

        id_counts = Counter()
        unmatched_labels = Counter()
        dropped_records = []  # list of (core_id, raw_label, reason)

        # Stats
        total_records = 0
        kept_records = 0
        dropped_count = 0
        unmatched_count = 0  # records that fell through to Repbase_FINAL_catchall

        # Buffered record-level processing: we must decide to keep or drop a
        # record BEFORE writing any of its lines, otherwise sequence bodies
        # could leak from a dropped header into the previous kept record.
        try:
            with open(self.input_fasta, 'r') as f_in, open(self.output_fasta, 'w') as f_out:

                current_header = None      # the new (rewritten) header line
                current_seq_lines = []     # buffered sequence body
                current_keep = False       # whether to flush this record

                def flush():
                    nonlocal kept_records
                    if current_header is not None and current_keep:
                        f_out.write(current_header + "\n")
                        for s in current_seq_lines:
                            f_out.write(s + "\n")
                        kept_records += 1

                for raw_line in f_in:
                    line = raw_line.rstrip("\r\n")
                    if not line:
                        continue

                    if line.startswith(">"):
                        # Flush the previous record (if any) before starting a new one
                        flush()
                        current_header = None
                        current_seq_lines = []
                        current_keep = False

                        total_records += 1

                        # Parse and normalize
                        core_id, raw_label = self._parse_header(line)
                        norm_class, norm_superfam = self.te_dict.normalize(
                            "Repbase", raw_label if raw_label else "Unknown"
                        )

                        # DROP filter: gene fragments, host genes, non-TE artifacts
                        if norm_class == "DROP":
                            dropped_records.append((core_id, raw_label, "gene_fragment_or_artifact"))
                            dropped_count += 1
                            current_keep = False
                            continue

                        # Track entries that hit the Repbase_FINAL_catchall
                        # (i.e. no specific rule matched; both fields are Unknown
                        # AND the original label was non-empty and non-Unknown).
                        if (norm_class.lower() == "unknown"
                                and norm_superfam.lower() == "unknown"
                                and raw_label
                                and raw_label.lower() != "unknown"):
                            unmatched_labels[raw_label] += 1
                            unmatched_count += 1

                        # Build the rewritten header and apply deduplication
                        new_header = self._build_new_header(core_id, norm_class, norm_superfam)
                        new_header = self._apply_dedup(new_header, id_counts)

                        current_header = new_header
                        current_keep = True
                    else:
                        # Sequence line: only retain if we plan to keep this record
                        if current_keep:
                            current_seq_lines.append(line)

                # Flush the last record
                flush()

        except Exception as e:
            print(f"[Error] Internal Python error while processing file: {e}")
            sys.exit(1)

        # Write the dropped-sequences log
        with open(self.dropped_log, 'w') as f:
            f.write("# Sequences dropped by RepbaseFixer (gene fragments / non-TE artifacts)\n")
            f.write("# Triggered by DROP/DROP rules in te_mapping_rules.tsv\n")
            f.write("# Columns: core_id\traw_label\treason\n")
            for core_id, raw_label, reason in dropped_records:
                f.write(f"{core_id}\t{raw_label}\t{reason}\n")

        # Write the unmatched-labels log
        with open(self.unmatched_log, 'w') as f:
            f.write("# Repbase labels that fell through to the FINAL catch-all rule\n")
            f.write("# These were normalized to Unknown/Unknown.\n")
            f.write("# Add explicit rules to te_mapping_rules.tsv if any of these matter.\n")
            f.write("# Columns: count\traw_label\n")
            for label, n in unmatched_labels.most_common():
                f.write(f"{n}\t{label}\n")

        # Summary report
        print("\n[RepbaseFixer Summary]")
        print(f"  Total records read:        {total_records}")
        print(f"  Kept and written:          {kept_records}")
        print(f"  Dropped (DROP/DROP rules): {dropped_count}")
        print(f"  Unmatched (->Unknown):     {unmatched_count} records, "
              f"{len(unmatched_labels)} distinct labels")
        print(f"  Dropped log:    {self.dropped_log}")
        print(f"  Unmatched log:  {self.unmatched_log}")

        # ID uniqueness verification (Python-native, no shell)
        print("Verifying ID uniqueness...")
        seen = Counter()
        with open(self.output_fasta, 'r') as f:
            for line in f:
                if line.startswith(">"):
                    sid = line[1:].split("#", 1)[0].strip()
                    seen[sid] += 1
        dups = [(sid, n) for sid, n in seen.items() if n > 1]
        if dups:
            print(f"[Warning] {len(dups)} duplicate IDs remain after dedup:")
            for sid, n in dups[:10]:
                print(f"  {sid}\t{n}")
        else:
            print("Verification passed: All sequence IDs are unique. [Check]")

        return self.output_fasta
