import os
import csv
import subprocess
import sys

from .normalizer import format_te_label


# ---------------------------------------------------------------------------
# Module-level helpers and constants
# ---------------------------------------------------------------------------

def _sf_compatible(sf_a, sf_b):
    """Are two confident, non-Unknown superfamily strings biologically compatible?

    Returns True when:
      - The strings are identical (EXACT match), OR
      - They are in a Dfam-canonical parent-child relationship: under the
        Dfam naming convention 'X' is the parent of 'X-Y' (e.g. 'CMC' is
        the parent of 'CMC-EnSpm', 'MULE' is the parent of 'MULE-MuDR',
        'CMC-Chapaev' is the parent of 'CMC-Chapaev-3'). Order does not
        matter: parent-vs-child and child-vs-parent both return True.

    Returns False for sibling-level disagreement ('CMC-EnSpm' vs
    'CMC-Transib') and completely unrelated superfamilies ('TcMar-Tc1'
    vs 'PiggyBac').

    Both inputs are expected to be already-normalized Dfam canonical
    superfamily labels (case-sensitive). 'Unknown' should be filtered
    BEFORE calling this function; the caller is responsible for treating
    Unknown as 'no information' rather than a competing label.
    """
    if sf_a == sf_b:
        return True
    if sf_a.startswith(sf_b + "-") or sf_b.startswith(sf_a + "-"):
        return True
    return False


# Autonomous TE classes that SHOULD carry detectable catalytic / structural
# protein domains (transposase / reverse transcriptase / integrase / RepHel).
# Used by Track B to exempt SHORT (<1000 bp) truncated copies of these
# classes from penalty: such fragments legitimately lack a TEsorter signal,
# so the absence of a domain is not held against them.
#
# Retroposon is INTENTIONALLY EXCLUDED. Retroposons (SVA, L1-derived,
# sno-derived, RTE-derived, etc.) are mobilized in trans by their parent
# autonomous element and are domain-less by definition, so the absence of a
# TEsorter domain is never meaningful evidence for or against them.
_DOMAIN_BEARING_CLASSES = frozenset({"LTR", "LINE", "DNA", "RC", "PLE"})


# ---------------------------------------------------------------------------
# TEsorter subprocess and parsing
# ---------------------------------------------------------------------------

def run_tesorter(fasta_file, outdir, threads, db, prefix):
    """
    Executes TEsorter via subprocess for domain-based diagnosis.
    Captures stdout and stderr to a log file for easier debugging if it fails.
    """
    tmp_dir = os.path.join(outdir, f"tmp_{prefix}_tesorter")
    os.makedirs(tmp_dir, exist_ok=True)
    abs_fasta = os.path.abspath(fasta_file)
    
    cmd = [
        "TEsorter", abs_fasta,
        "-db", db,
        "-p", str(threads),
        "-pre", prefix,
        "-tmp", tmp_dir
    ]
    
    log_file = os.path.join(outdir, f"{prefix}_tesorter.log")
    print(f"[{prefix}] Running TEsorter diagnosis...")
    
    with open(log_file, "w") as log_f:
        try:
            # Route both stdout and stderr to the log file instead of DEVNULL
            subprocess.check_call(cmd, cwd=outdir, stdout=log_f, stderr=subprocess.STDOUT)
            print(f"[{prefix}] TEsorter diagnosis completed.")
        except subprocess.CalledProcessError:
            print(f"\n[Error] TEsorter failed to execute for {prefix}.", file=sys.stderr)
            print(f"Detailed error log saved to: {log_file}", file=sys.stderr)
            
            # Attempt to read and print the last 10 lines of the log for quick diagnosis
            try:
                with open(log_file, "r") as read_log:
                    lines = read_log.readlines()
                    if lines:
                        print("-" * 50, file=sys.stderr)
                        print("--- Last 10 lines of TEsorter log ---", file=sys.stderr)
                        for line in lines[-10:]:
                            print(line.strip(), file=sys.stderr)
                        print("-" * 50, file=sys.stderr)
            except Exception:
                pass # Fail silently if the log file cannot be read
                
            sys.exit(1)


def parse_tesorter_cls(cls_path):
    """
    Parses TEsorter's .cls.tsv output to extract class and superfamily.
    Only Order (col 2) and Superfamily (col 3) are retained; the Clade
    column is intentionally dropped because RLB's dictionary patterns
    operate on the 2-segment 'Order/Superfamily' form.
    """
    tesorter_info = {}
    if not os.path.exists(cls_path):
        return tesorter_info
        
    with open(cls_path, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        for row in reader:
            if not row or row[0].startswith('#'):
                continue
            
            full_id = row[0]
            # Strip suffix if present to match the core sequence ID
            pure_id = full_id.split('#')[0] if '#' in full_id else full_id
            order = row[1]
            superfamily = row[2]
            tesorter_info[pure_id] = f"{order}/{superfamily}"
            
    return tesorter_info


# ---------------------------------------------------------------------------
# Main diagnosis + scoring + repair
# ---------------------------------------------------------------------------

def diagnose_score_and_repair(software_name, seq_info, tesorter_cls_path, te_dict):
    """
    Applies Track A/B logic for the Confidence Index and strict rules for Active Repair.
    Ensures high-resolution preservation for unmapped but valid TE labels.

    Track A scoring (parent-child aware):
        Once Class is confirmed equal between Denovo and TEsorter, the
        superfamily outcomes are scored as follows:
          (1) Both sides confident AND compatible at sf          -> +3 Perfect
                Compatibility means either:
                  (a) Exact string match.
                  (b) Hierarchical parent-child relationship under
                      the Dfam naming convention ('X' parent of
                      'X-Y'). Common when EDTA's DT* parent rollup
                      meets TEsorter's confident child clade, or
                      when a source library's Repbase-derived child
                      meets TEsorter's parent-only resolution.
                Both routes represent positive mutual evidence with
                no contradiction, so they share the strongest tier.
          (2) At least one side sf=Unknown (no refutation)       -> +2 Asymmetric
                Two sub-cases share this tier:
                  (a) Denovo confident sf, TEsorter sf=Unknown.
                      TEsorter cannot detect a domain (e.g.
                      truncated copy), so it cannot refute the
                      Denovo claim. The repair stage already
                      keeps the Denovo superfamily here via the
                      Protection Shield.
                  (b) Both sides sf=Unknown.
                      Class-level consensus only; the sf layer
                      is untestable on either side. Treated as
                      'neither side can refute the other'.
                Condition simplifies to: ts_sf_unknown is True.
          (3) Denovo sf=Unknown, TEsorter sf confident           -> +1 Fuzzy
                TEsorter unilaterally rescues; weaker consensus.
          (4) Both confident but INCOMPATIBLE at sf              -> +1 Fuzzy
                Sibling-level disagreement (e.g. 'CMC-EnSpm' vs
                'CMC-Transib') or completely unrelated
                superfamilies (e.g. 'TcMar-Tc1' vs 'PiggyBac').
                Both sides confident but the claims contradict.

    Track B scoring:
        For sequences WITHOUT a TEsorter hit, scoring is structural only:
        short SINE/MITEs get a small positive exemption, short truncated
        copies of domain-bearing classes (LTR, LINE, DNA, RC, PLE) are
        exempted at 0, and everything else is retained as Unverified at 0.
        No sequence is penalised for lacking a domain, and none is dropped.

    Parent-child guard in repair branch:
        The Confirmed-repair branch would otherwise overwrite the source
        superfamily with TEsorter's whenever they differed. That
        silently DOWNGRADES resolution when TEsorter only resolved
        to the parent (e.g. orig='CMC-EnSpm' overwritten by ts='CMC').
        A parent check routes such downgrade cases into
        the Protection Shield instead of overwriting.

        Relatedly, formatted_label uses the normalized
        Class (orig_norm_class) instead of the pre-normalization
        raw_class prefix. Otherwise, any source where the dictionary
        REMAPPED the Class layer (e.g. EDTA 'MITE/DTM' -> DNA, HiTE
        'DNA/Helitron' -> RC, 'LINE/Penelope' -> PLE, 'DNAnona/*' ->
        RC) would produce FASTA headers whose Class disagreed with the
        internal new_norm_class field.

    Parent-child aware Track A:
        Track A's 'both sides confident' branch uses
        _sf_compatible() which treats EXACT and PARENT_CHILD as the
        same Perfect (+3) tier, since both are positive mutual
        evidence with no contradiction. This eliminates the
        systematic false-negative where EDTA's DT* parent rollup vs
        TEsorter's confident child (e.g. EDTA 'DNA/DTM' vs TEsorter
        'TIR/MuDR_Mutator', both normalizing into the same Dfam
        lineage) would otherwise be scored as Fuzzy +1 despite being
        biologically consistent.

        Track B's truncation-exemption class set is centralized into
        the module-level _DOMAIN_BEARING_CLASSES constant and covers
        RC (Helitron) and PLE (Penelope-like) alongside LTR/LINE/DNA.
        Retroposon is intentionally excluded (domain-less by design).
    """
    tesorter_info = parse_tesorter_cls(tesorter_cls_path)
    
    repair_plan = {}
    repair_stats = {
        'Confirmed': 0, 'Corrected': 0, 'Recovered': 0, 
        'Exempted': 0, 'Unverified': 0
    }
    
    total_score = 0
    scored_count = 0 
    total_sequences = len(seq_info)
    
    score_stats = {
        'TrackA_Perfect_Plus3': 0,
        'TrackA_Asymmetric_Plus2': 0,  # Covers (a) Denovo confident sf + TEsorter sf=Unknown, (b) both sf=Unknown
        'TrackA_Fuzzy_Plus1': 0,
        'TrackA_Miss_Minus1': 0,
        'TrackA_Conflict_Minus2': 0,
        'TrackB_SINE_Exempt_Plus0_5': 0,
        'TrackB_Trunc_Exempt_0': 0,
        'TrackB_Unknown_0': 0,
    }

    for seq_id, info in seq_info.items():
        # These are now safely preserved by the updated normalizer.py fallback mechanism
        orig_norm_class = info['norm_class']
        orig_norm_superfam = info['norm_superfam']
        raw_class = info['raw_class']
        seq_len = info['length']
        
        raw_class_upper = raw_class.upper()
        is_mite = "MITE" in raw_class_upper
        
        new_norm_class = orig_norm_class
        new_norm_superfam = orig_norm_superfam
        
        # Default initialization: the STANDARDIZED form of the source's own
        # label. Branches that do not override this (Unverified, Exempted,
        # and Confirmed-via-Protection-Shield) therefore emit a
        # canonical label instead of leaking the raw, pre-normalization string
        # into the repaired FASTA and the final build library. The normalized
        # default is applied to every branch here, not only the
        # Confirmed-overwrite branch.
        formatted_label = format_te_label(orig_norm_class, orig_norm_superfam)
        tag = "Unverified"
        score_changed = False
        
        has_tesorter = seq_id in tesorter_info
        is_orig_unknown = orig_norm_class.lower() == "unknown"
        
        # --- Dual Engine: Track Scoring and Active Diagnosis ---
        if has_tesorter:
            tesorter_raw = tesorter_info[seq_id]
            ts_class, ts_superfam = te_dict.normalize('TEsorter', tesorter_raw)
            is_ts_valid = ts_class.lower() != "unknown"
            
            # 1. Track A Scoring (Confidence Index)
            if not is_orig_unknown and orig_norm_class == ts_class:
                # Class layer agrees. Now arbitrate at the superfamily layer,
                # treating "Unknown" as missing information rather than a
                # competing label.
                orig_sf_unknown = orig_norm_superfam.lower() == "unknown"
                ts_sf_unknown   = ts_superfam.lower() == "unknown"

                if (not orig_sf_unknown) and (not ts_sf_unknown):
                    # Both sides have confident sf-level evidence.
                    # Compatibility includes both EXACT match and
                    # Dfam-canonical parent-child relationship.
                    if _sf_compatible(orig_norm_superfam, ts_superfam):
                        # Perfect (+3): exact string match OR hierarchical
                        # parent-child agreement. Common cases routed here
                        # (previously misclassified as Fuzzy +1):
                        #   - EDTA 'DNA/DTM' (->DNA/MULE) vs TEsorter
                        #     'TIR/MuDR_Mutator' (->DNA/MULE-MuDR)
                        #   - HiTE 'DNA/CMC-EnSpm' vs TEsorter 'TIR/CMC'
                        #   - Repbase 'MuDR-N1_AT' (->DNA/MULE-MuDR) vs
                        #     TEsorter 'TIR/MULE' (->DNA/MULE)
                        total_score += 3
                        score_stats['TrackA_Perfect_Plus3'] += 1
                    else:
                        # Fuzzy (+1): sibling-level disagreement (e.g.
                        # 'CMC-EnSpm' vs 'CMC-Transib') or unrelated
                        # superfamilies (e.g. 'TcMar-Tc1' vs 'PiggyBac').
                        # Both sides confident but the claims contradict.
                        total_score += 1
                        score_stats['TrackA_Fuzzy_Plus1'] += 1
                elif ts_sf_unknown:
                    # Asymmetric (+2): TEsorter has no sf-level evidence,
                    # so it cannot refute whatever the other side claims
                    # (or fails to claim). Covers both:
                    #   (a) Denovo confident sf + TEsorter sf=Unknown
                    #   (b) Both sides sf=Unknown
                    # The Protection Shield in the repair branch keeps the
                    # Denovo sf in case (a) and leaves the label untouched
                    # in case (b), so scoring stays consistent with repair.
                    total_score += 2
                    score_stats['TrackA_Asymmetric_Plus2'] += 1
                else:
                    # Fuzzy (+1): Denovo sf=Unknown + TEsorter confident sf.
                    # TEsorter unilaterally rescues (weaker consensus,
                    # since Denovo gave up).
                    total_score += 1
                    score_stats['TrackA_Fuzzy_Plus1'] += 1
                score_changed = True
            elif is_orig_unknown and is_ts_valid:
                total_score -= 1
                score_stats['TrackA_Miss_Minus1'] += 1
                score_changed = True
            elif is_orig_unknown and not is_ts_valid:
                # Edge Case Protection: Both sides are Unknown. No points added/subtracted.
                pass
            elif orig_norm_class != ts_class:
                total_score -= 2
                score_stats['TrackA_Conflict_Minus2'] += 1
                score_changed = True
                
            # 2. Strict Active Repair Logic (The 4 Conditions)
            if is_orig_unknown and is_ts_valid:
                # Cond 4: Recovered (Original Unknown -> Domain Identified)
                tag = "Recovered"
                new_norm_class, new_norm_superfam = ts_class, ts_superfam
                formatted_label = format_te_label(ts_class, ts_superfam)
                
            elif not is_orig_unknown and is_ts_valid and orig_norm_class != ts_class:
                # Cond 1: Corrected (Major class conflict -> Forced Override)
                tag = "Corrected"
                new_norm_class, new_norm_superfam = ts_class, ts_superfam
                formatted_label = format_te_label(ts_class, ts_superfam)
                
            elif not is_orig_unknown and orig_norm_class == ts_class:
                # Cond 2 & 3: Confirmed (Perfect class match, evaluate subclass)
                tag = "Confirmed"
                is_ts_sf_present = ts_superfam.lower() not in ["unknown", ""]

                if is_ts_sf_present and orig_norm_superfam != ts_superfam:
                    # Parent-vs-child guard:
                    # Dfam canonical naming convention encodes hierarchy as
                    # 'PARENT' vs 'PARENT-CHILD' (e.g. CMC vs CMC-EnSpm,
                    # MULE vs MULE-MuDR, hAT vs hAT-Charlie). If TEsorter
                    # only resolved to the parent while the source library
                    # already pinned down a child, the unconditional
                    # overwrite of the original code was DOWNGRADING
                    # resolution. The check below routes that case into
                    # the Protection Shield instead.
                    # - Sibling-level disagreement (e.g. CMC-EnSpm vs
                    #   CMC-Transib) and unrelated-superfamily disagreement
                    #   still fall through to overwrite, preserving the
                    #   original "TEsorter wins on disagreement" behaviour.
                    is_ts_parent_of_orig = orig_norm_superfam.startswith(ts_superfam + "-")
                    if not is_ts_parent_of_orig:
                        new_norm_superfam = ts_superfam
                        # Use normalized Class, not raw_class prefix.
                        # Otherwise sources whose dictionary remaps the
                        # Class layer (MITE->DNA, DNA/Helitron->RC,
                        # LINE/Penelope->PLE, DNAnona/*->RC, snoRNA->
                        # Retroposon, Evirus/*->LTR, ...) end up writing
                        # a Class in the FASTA header that disagrees with
                        # new_norm_class and is non-canonical to
                        # RepeatMasker.
                        formatted_label = format_te_label(orig_norm_class, ts_superfam)
                    # else: Protection Shield (TEsorter is parent, keep child)
                else:
                    # Protection Shield: Retain original high-resolution string if TEsorter lacks subclass
                    pass
            else:
                # Catch-all for "Unknown == Unknown" collisions or other unverified states
                tag = "Unverified"
                
        else:
            # 1. Track B Scoring & Matrix Diagnosis (No Domains Detected)
            is_sine = orig_norm_class.upper() == "SINE"
            
            if (is_sine or is_mite) and seq_len < 1000:
                # Matrix D: Exempted (Short structural TEs)
                total_score += 0.5  
                score_stats['TrackB_SINE_Exempt_Plus0_5'] += 1
                score_changed = True
                tag = "Exempted"
            elif orig_norm_class in _DOMAIN_BEARING_CLASSES and seq_len < 1000:
                # Truncated copies of domain-bearing autonomous classes
                # legitimately lack a TEsorter signal, so they are exempted
                # from penalty and simply retained as Unverified.
                score_stats['TrackB_Trunc_Exempt_0'] += 1
                tag = "Unverified"
            else:
                score_stats['TrackB_Unknown_0'] += 1
                tag = "Unverified"

        if score_changed:
            scored_count += 1
            
        repair_stats[tag] += 1
        
        # Final cleanup for edge case outputs
        if formatted_label.lower() in ["unknown/unknown", "unknown"]:
            formatted_label = "Unknown"
            
        # Construct the finalized diagnostic tagged header
        new_header = f"{software_name}_{tag}_{seq_id}#{formatted_label}"
        
        repair_plan[seq_id] = {
            'new_header': new_header,
            'new_class': new_norm_class,
            'new_superfam': new_norm_superfam,
            'tag': tag
        }

    # Calculate final Confidence Metrics
    norm_score = (total_score / total_sequences) * 100 if total_sequences > 0 else 0
    mean_scored = (total_score / scored_count) if scored_count > 0 else 0
    
    confidence_metrics = {
        "Absolute_Score": round(total_score, 2),
        "Normalized_Score": round(norm_score, 2),
        "Mean_Scored": round(mean_scored, 2),
        "Scored_Count": scored_count,
        "Total_Count": total_sequences,
        "Detailed_Stats": score_stats
    }
        
    return repair_plan, repair_stats, confidence_metrics
