import os
import json
from .normalizer import TEDictionary, process_and_normalize_fasta
from .metrics import calculate_basic_metrics
from .tesorter_scorer import run_tesorter, diagnose_score_and_repair

def run_evaluate_pipeline(args, dirs):
    """
    Main entry point for evaluating de novo TE libraries.
    Distributes outputs to the proper subdirectories in 02.denovo_db.
    """
    te_dict = TEDictionary(args.mapping)
    
    inputs = {
        "HiTE": getattr(args, 'hite', None),
        "EDTA": getattr(args, 'edta', None),
        "RM2": getattr(args, 'rm2', None)
    }
    
    final_report = {}
    repaired_files = {}

    for software, fasta_path in inputs.items():
        if fasta_path and os.path.exists(fasta_path):
            print(f"\n---> Analyzing {software} Library...")
            
            # Route to correct output subdirectory
            out_subdir = dirs.get(f"denovo_{software.lower()}")
            clean_fasta_path = os.path.join(out_subdir, f"{software}_clean.fa")
            
            # Step 1: Normalization
            print(f"  [{software}] 1. Normalizing labels and cleaning FASTA headers...")
            seq_info = process_and_normalize_fasta(fasta_path, clean_fasta_path, software, te_dict)
            
            # Step 2: Pre-Repair Metrics
            print(f"  [{software}] 2. Calculating Pre-Repair Baseline Metrics...")
            pre_repair_metrics = calculate_basic_metrics(seq_info)
            
            # Step 3: Run TEsorter
            print(f"  [{software}] 3. Executing TEsorter for domain-based diagnosis...")
            run_tesorter(clean_fasta_path, out_subdir, args.threads, args.te_db, software)
            cls_path = os.path.join(out_subdir, f"{software}.cls.tsv")
            
            # Step 4: Active Repair & Scoring
            print(f"  [{software}] 4. Applying Active Repair Rules & Track Scoring...")
            repair_plan, repair_stats, confidence_metrics = diagnose_score_and_repair(
                software, seq_info, cls_path, te_dict
            )
            
            # Step 5: Write Repaired FASTA
            repaired_fasta_path = os.path.join(out_subdir, f"{software}_repaired.fa")
            print(f"  [{software}] 5. Generating finalized sequence library: {repaired_fasta_path}")
            
            repaired_seq_info = {}
            with open(clean_fasta_path, 'r') as fin, open(repaired_fasta_path, 'w') as fout:
                for line in fin:
                    if line.startswith(">"):
                        raw_id = line.strip()[1:]
                        pure_id = raw_id.split('#')[0] if '#' in raw_id else raw_id
                        
                        if pure_id in repair_plan:
                            plan = repair_plan[pure_id]
                            fout.write(f">{plan['new_header']}\n")
                            repaired_seq_info[plan['new_header']] = {
                                'length': seq_info[pure_id]['length'],
                                'norm_class': plan['new_class'],
                                'norm_superfam': plan['new_superfam'],
                                'tag': plan['tag'] 
                            }
                        else:
                            fout.write(line)
                    else:
                        fout.write(line)
                        
            # Register path for the merger module
            repaired_files[software.lower()] = repaired_fasta_path
            
            # Step 6: Post-Repair Metrics
            print(f"  [{software}] 6. Calculating Post-Repair Health Metrics...")
            post_repair_metrics = calculate_basic_metrics(repaired_seq_info)
            
            total_seqs = pre_repair_metrics["Yield"]["Total_Sequences"]
            corrected_count = repair_stats.get("Corrected", 0)
            repair_rate = round((corrected_count / total_seqs) * 100, 2) if total_seqs > 0 else 0.0

            final_report[software] = {
                "1_Confidence_Index": confidence_metrics,
                "2_Repair_Statistics": {
                    "Counts": repair_stats,
                    "Repair_Rate_Percent": repair_rate
                },
                "3_Pre_Repair_Yield": pre_repair_metrics["Yield"],
                "4_Post_Repair_Yield": post_repair_metrics["Yield"],
                "5_Post_Repair_Lengths": post_repair_metrics["Lengths"],
                "6_Payload_Score_Percent": post_repair_metrics["Payload_Score"]
            }

    # Step 7: Output Report
    if final_report:
        report_path = os.path.join(dirs['report'], "evaluation_report.json")
        with open(report_path, 'w') as f:
            json.dump(final_report, f, indent=4)
        print(f"\n[Info] De novo diagnostics saved to: {report_path}")

    return repaired_files
