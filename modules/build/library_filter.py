import os
import sys
import argparse

# Dynamically import TEDictionary from step 1
BASE_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
if BASE_DIR not in sys.path:
    sys.path.append(BASE_DIR)

from modules.rating_module.normalizer import TEDictionary

class LibraryFilter:
    def __init__(self, input_fasta, output_fasta, te_dict, filter_db=None, filter_tag=None, filter_length=None):
        self.input_fasta = input_fasta
        self.output_fasta = output_fasta
        self.te_dict = te_dict
        
        # Initialize storage structures for filtering conditions
        self.drop_tags = set()
        self.drop_db_class = {}  # Format: { 'RM2': set(['unknown', 'sine']), 'EDTA': set(['unknown']) }
        self.length_limits = {}  # Format: { 'EDTA': { 'sine': {'min': 0, 'max': 1000} } }
        
        self.stats = {
            'Total': 0, 
            'Passed': 0, 
            'Filtered_by_Tag': 0, 
            'Filtered_by_DB': 0, 
            'Filtered_by_Length': 0
        }
        
        self._parse_arguments(filter_db, filter_tag, filter_length)

    def _parse_arguments(self, filter_db, filter_tag, filter_length):
        """Parses user input strings into highly efficient O(1) lookup dictionaries."""
        # 1. Parse Tag (e.g., "Unverified,Exempted")
        if filter_tag:
            self.drop_tags = set([t.strip() for t in filter_tag.split(',') if t.strip()])

        # 2. Parse DB:Class (e.g., "RM2:Unknown,EDTA:Unknown")
        if filter_db:
            for item in filter_db.split(','):
                parts = item.strip().split(':')
                if len(parts) == 2:
                    db, te_class = parts[0], parts[1]
                    if db not in self.drop_db_class:
                        self.drop_db_class[db] = set()
                    self.drop_db_class[db].add(te_class.lower())

        # 3. Parse Length Limits (e.g., "EDTA:SINE:max:1000,RM2:DNA:min:50")
        if filter_length:
            for item in filter_length.split(','):
                parts = item.strip().split(':')
                if len(parts) == 4:
                    db, te_class, limit_type, value = parts[0], parts[1], parts[2], int(parts[3])
                    if db not in self.length_limits:
                        self.length_limits[db] = {}
                    if te_class.lower() not in self.length_limits[db]:
                        self.length_limits[db][te_class.lower()] = {}
                    self.length_limits[db][te_class.lower()][limit_type.lower()] = value

    def _read_fasta_generator(self):
        """Yields FASTA records one by one to avoid memory exhaustion on huge libraries."""
        header = ""
        seq_buffer = []
        with open(self.input_fasta, 'r') as f:
            for line in f:
                line = line.strip()
                if not line: 
                    continue
                if line.startswith(">"):
                    if header:
                        yield header, "".join(seq_buffer)
                    header = line
                    seq_buffer = []
                else:
                    seq_buffer.append(line)
            if header:
                yield header, "".join(seq_buffer)

    def _parse_header(self, header):
        """
        Extracts key diagnostic info from the standardized diagnostic header.
        """
        # Example: >RM2_Confirmed_rnd-5#SINE/tRNA-V
        clean_header = header[1:].split()[0]  
        
        if '#' in clean_header:
            id_part, raw_class = clean_header.split('#', 1)
        else:
            id_part, raw_class = clean_header, "Unknown"

        # Split at the first underscore only to isolate Software from the rest
        parts = id_part.split('_', 1)  
        software = parts[0]
        tag = None
        
        # Handle software containing tags
        if software in ['RM2', 'EDTA', 'HiTE'] and len(parts) > 1:
            tag = parts[1].split('_', 1)[0]  # Extracts Confirmed, Unverified, etc.
                
        # Utilize TEDictionary to translate complex IDs (e.g., MITE/DTC -> DNA)
        norm_class, _ = self.te_dict.normalize(software, raw_class)
        
        return software, tag, norm_class.lower()

    def run(self):
        """Executes the filtering cascade."""
        print(f"\n[Module: LibraryFilter] Initializing...")
        
        with open(self.output_fasta, 'w') as fout:
            for header, seq in self._read_fasta_generator():
                self.stats['Total'] += 1
                seq_len = len(seq)
                software, tag, norm_class = self._parse_header(header)
                
                # --- Condition 1: Tag Filtering ---
                if tag and tag in self.drop_tags:
                    self.stats['Filtered_by_Tag'] += 1
                    continue
                    
                # --- Condition 2: DB + Class Filtering ---
                if software in self.drop_db_class and norm_class in self.drop_db_class[software]:
                    self.stats['Filtered_by_DB'] += 1
                    continue
                    
                # --- Condition 3: Length Filtering ---
                if software in self.length_limits and norm_class in self.length_limits[software]:
                    limits = self.length_limits[software][norm_class]
                    if 'max' in limits and seq_len > limits['max']:
                        self.stats['Filtered_by_Length'] += 1
                        continue
                    if 'min' in limits and seq_len < limits['min']:
                        self.stats['Filtered_by_Length'] += 1
                        continue

                # Passed all filters - Write output elegantly formatted
                fout.write(f"{header}\n")
                for i in range(0, seq_len, 50):
                    fout.write(seq[i:i+50] + "\n")
                self.stats['Passed'] += 1

        print("-" * 50)
        print("Filtering Summary:")
        print(f"  Total Input Sequences : {self.stats['Total']:,}")
        print(f"  Filtered by Tag       : {self.stats['Filtered_by_Tag']:,}")
        print(f"  Filtered by DB/Class  : {self.stats['Filtered_by_DB']:,}")
        print(f"  Filtered by Length    : {self.stats['Filtered_by_Length']:,}")
        print(f"  Passed & Saved        : {self.stats['Passed']:,}")
        print("-" * 50)
        print(f"Filtered library saved to: {self.output_fasta}\n")

if __name__ == "__main__":
    # For standalone testing if needed
    parser = argparse.ArgumentParser(description="Filter merged TE library based on Step 1 diagnostic tags and classes.")
    parser.add_argument("-i", "--input", required=True, help="Input merged raw FASTA")
    parser.add_argument("-o", "--output", required=True, help="Output filtered FASTA")
    parser.add_argument("-m", "--mapping", required=True, help="Path to te_mapping_rules.tsv")
    parser.add_argument("--filter_db", help="e.g., RM2:Unknown,EDTA:Unknown")
    parser.add_argument("--filter_tag", help="e.g., Unverified,Exempted")
    parser.add_argument("--filter_length", help="e.g., EDTA:SINE:max:1000,RM2:DNA:min:50")
    
    args = parser.parse_args()
    te_dict = TEDictionary(args.mapping)
    filter_module = LibraryFilter(args.input, args.output, te_dict, args.filter_db, args.filter_tag, args.filter_length)
    filter_module.run()
