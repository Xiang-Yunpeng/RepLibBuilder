import os
import json
from utils.common import ensure_dir

class RawLibraryMerger:
    def __init__(self, inputs, output_dir, report_dir):
        """
        Initializes the pre-clustering merger for Step 1.
        Args:
            inputs (dict): Dictionary containing paths to processed fasta files:
                           {'repbase': path, 'dfam': path, 'hite': path, 'edta': path, 'rm2': path}
            output_dir (str): Path to '03.merge_db/'
            report_dir (str): Path to '04.report/' for dumping merge statistics
        """
        self.inputs = inputs
        self.output_dir = ensure_dir(output_dir)
        self.report_dir = ensure_dir(report_dir)
        
        self.merged_raw = os.path.join(self.output_dir, "step1_combined_raw.fa")
        self.stats_json = os.path.join(self.report_dir, "merge_stats.json")

    def run(self):
        """Executes the raw merge and gathers basic tagging statistics."""
        print("\n" + "="*60)
        print("[Module: RawLibraryMerger] Assembling Combined Library")
        print("="*60)
        
        stats = {
            'Public_DB': {'Repbase': 0, 'Dfam': 0},
            'Denovo_DB': {
                'HiTE': {'Confirmed': 0, 'Corrected': 0, 'Recovered': 0, 'Exempted': 0, 'Unverified': 0},
                'EDTA': {'Confirmed': 0, 'Corrected': 0, 'Recovered': 0, 'Exempted': 0, 'Unverified': 0},
                'RM2':  {'Confirmed': 0, 'Corrected': 0, 'Recovered': 0, 'Exempted': 0, 'Unverified': 0}
            },
            'Total_Sequences': 0
        }

        with open(self.merged_raw, 'w') as fout:
            # 1. Process Public Databases
            for pub in ['repbase', 'dfam']:
                path = self.inputs.get(pub)
                if path and os.path.exists(path):
                    print(f"  --> Appending {pub.capitalize()}...")
                    count = 0
                    with open(path, 'r') as fin:
                        for line in fin:
                            fout.write(line)
                            if line.startswith(">"):
                                count += 1
                    stats['Public_DB'][pub.capitalize()] = count
                    stats['Total_Sequences'] += count

            # 2. Process Denovo Databases
            for denovo in ['hite', 'edta', 'rm2']:
                path = self.inputs.get(denovo)
                db_name = 'RM2' if denovo == 'rm2' else denovo.upper()
                if denovo == 'hite': db_name = 'HiTE'
                
                if path and os.path.exists(path):
                    print(f"  --> Appending {db_name} and parsing diagnostic tags...")
                    with open(path, 'r') as fin:
                        for line in fin:
                            fout.write(line)
                            if line.startswith(">"):
                                stats['Total_Sequences'] += 1
                                
                                # Expected header format: >Software_Tag_ID#Class/Superfamily
                                # Example: >RM2_Confirmed_rnd-5...
                                header = line[1:].strip()
                                parts = header.split('_', 2)
                                
                                if len(parts) >= 2:
                                    tag = parts[1]
                                    if tag in stats['Denovo_DB'][db_name]:
                                        stats['Denovo_DB'][db_name][tag] += 1

        # Dump merge statistics to JSON for the HTML reporter
        with open(self.stats_json, 'w') as f:
            json.dump(stats, f, indent=4)
            
        print(f"\n[Raw Merge Complete]")
        print(f"Total Sequences Pooled: {stats['Total_Sequences']}")
        print(f"Output saved to: {self.merged_raw}")
        print(f"Merge statistics saved to: {self.stats_json}")
        
        return self.merged_raw
