import os
import sys
import re
from utils.common import run_command, check_file, ensure_dir
from .normalizer import format_te_label

class DfamFetcher:
    def __init__(self, dfam_db, tax_id, build_rm_path, output_dir, te_dict):
        """
        Args:
            dfam_db (str): Path to Dfam database file (h5).
            tax_id (str): NCBI Taxonomy ID.
            build_rm_path (str): [Obsolete] Maintained for argument compatibility only. Ignored.
            output_dir (str): Output directory.
            te_dict (TEDictionary): Pre-loaded TE dictionary used to normalize the
                                    Dfam Type/SubType labels through the same gate
                                    as every other source (identity mapping for
                                    canonical Dfam labels).
        """
        self.dfam_db = dfam_db
        self.tax_id = str(tax_id)
        self.te_dict = te_dict
        
        # Dynamically resolve the internal famdb.py path
        current_module_dir = os.path.dirname(os.path.abspath(__file__))
        parent_module_dir = os.path.dirname(current_module_dir)
        internal_famdb_path = os.path.join(parent_module_dir, 'famdb', 'famdb.py')
        self.famdb_path = check_file(internal_famdb_path, "Internal famdb.py script")
        
        self.output_dir = ensure_dir(output_dir)
        self.famdb_python = sys.executable

        self.output_fasta = os.path.join(self.output_dir, f"Dfam_taxID_{self.tax_id}_clean.fa")
        self.temp_embl = os.path.join(self.output_dir, "temp_Dfam.embl")
        self.temp_raw_fa = os.path.join(self.output_dir, "temp_Dfam_raw.fa")

    def _embl_to_fasta(self, embl_file, fasta_file):
        """
        100% Python port of RepeatMasker's buildRMLibFromEMBL.pl.
        Pixel-perfect output matching the original Perl script, with Dfam_ prefix added.
        """
        with open(embl_file, 'r') as f_in, open(fasta_file, 'w') as f_out:
            record_lines = []
            
            def process_record(lines):
                if not lines: return
                
                rm_id = ""
                rm_type = ""
                rm_subtype = ""
                rm_desc = ""
                rm_species_list = []
                rm_search_stages = []
                rm_buffer_stages = []
                seq_buffer = []
                in_sq = False
                
                for line in lines:
                    # Strip ONLY line breaks to preserve any structural leading/trailing spaces
                    line_str = line.strip('\r\n')
                    if not line_str: continue
                    
                    if in_sq:
                        # Only keep A-Z, a-z (removes spaces, block numbers, etc.)
                        clean_seq = re.sub(r'[^A-Za-z]', '', line_str)
                        seq_buffer.append(clean_seq)
                        continue
                        
                    if line_str.startswith('ID '):
                        parts = line_str.split()
                        if len(parts) > 1:
                            rm_id = parts[1].strip(';')
                    elif line_str.startswith('DE '):
                        # Append multiple lines of descriptions perfectly
                        new_desc = line_str[3:].strip()
                        rm_desc = f"{rm_desc} {new_desc}" if rm_desc else new_desc
                    else:
                        m_type = re.match(r'^CC\s+Type:\s*(.*)', line_str)
                        if m_type:
                            rm_type = m_type.group(1).strip()
                            continue
                            
                        m_subtype = re.match(r'^CC\s+SubType:\s*(.*)', line_str)
                        if m_subtype:
                            rm_subtype = m_subtype.group(1).strip()
                            continue
                            
                        m_species = re.match(r'^CC\s+Species:\s*(.*)', line_str)
                        if m_species:
                            sp_raw = m_species.group(1).strip()
                            for sp in sp_raw.replace(';', ' ').split(','):
                                sp_clean = sp.strip()
                                if sp_clean: rm_species_list.append(sp_clean)
                            continue
                            
                        m_stages = re.match(r'^CC\s+SearchStages:\s*(.*)', line_str)
                        if m_stages:
                            st_raw = m_stages.group(1).strip()
                            for st in st_raw.split(','):
                                st_clean = st.strip()
                                if st_clean: rm_search_stages.append(st_clean)
                            continue
                            
                        m_buffer = re.match(r'^CC\s+BufferStages:\s*(.*)', line_str)
                        if m_buffer:
                            bs_raw = m_buffer.group(1).strip()
                            for bs in bs_raw.split(','):
                                bs_clean = bs.strip()
                                if bs_clean: rm_buffer_stages.append(bs_clean)
                            continue
                            
                        if line_str.startswith('SQ '):
                            in_sq = True
                
                if not rm_id:
                    return
                
                # 1. Normalize the Type/SubType through the dictionary so that
                #    the Dfam classification label passes through the same gate
                #    as every other source. For canonical Dfam labels this is an
                #    identity mapping, but it guarantees a uniform
                #    "Class/Superfamily" format across all libraries. The
                #    structural buffer records below are NOT TE class labels and
                #    are intentionally left untouched.
                raw_label = f"{rm_type}/{rm_subtype}" if rm_subtype else rm_type
                norm_class, norm_superfam = self.te_dict.normalize(
                    "Dfam", raw_label if raw_label else "Unknown"
                )
                type_str = f"#{format_te_label(norm_class, norm_superfam)}"
                    
                # 2. Recreate Perl's Species logic (Includes trailing space per item)
                species_str = ""
                for sp in rm_species_list:
                    species_str += f"@{sp} "
                    
                # 3. Recreate Perl's Search Stages string handling
                stage_str = "[S:"
                for st in rm_search_stages:
                    stage_str += f"{st},"
                stage_str += "]"
                stage_str = stage_str.replace(",]", "]")
                
                # 4. Retain original case (Do not force .upper() as Perl doesn't either)
                full_seq = "".join(seq_buffer)
                
                def write_fasta(header_id, h_type, h_species, h_stages, h_desc, sequence):
                    # Ensure the ID starts with 'Dfam_'
                    if not header_id.startswith("Dfam_"):
                        header_id = f"Dfam_{header_id}"
                        
                    # Pixel-perfect Perl string interpolation recreation:
                    header_line = f">{header_id}{h_type} {h_species} {h_stages} {h_desc}"
                    f_out.write(header_line + "\n")
                    
                    # Wrap sequence at 50 base pairs
                    for i in range(0, len(sequence), 50):
                        f_out.write(sequence[i:i+50] + "\n")
                
                # Write standard sequence
                write_fasta(rm_id, type_str, species_str, stage_str, rm_desc, full_seq)
                
                # Process and write buffer sequences
                stage_hash = {}
                for bs in rm_buffer_stages:
                    m_range = re.match(r'(\d+)\[(\d+)-(\d+)\]', bs)
                    if m_range:
                        stage, start, end = m_range.groups()
                        key = f"{start}-{end}"
                        if key not in stage_hash: stage_hash[key] = []
                        stage_hash[key].append(stage)
                    else:
                        m_full = re.match(r'(\d+)', bs)
                        if m_full:
                            stage = m_full.group(1)
                            key = "full"
                            if key not in stage_hash: stage_hash[key] = []
                            stage_hash[key].append(stage)
                
                for key, stages in stage_hash.items():
                    b_st_joined = ",".join(stages)
                    b_stages_str = f"[S:{b_st_joined}]"
                    
                    if key == "full":
                        b_type = "#buffer"
                        b_seq = full_seq
                    else:
                        start, end = map(int, key.split('-'))
                        b_type = f"_{start}_{end}#buffer"
                        # Translate Perl's 1-based index to Python's 0-based
                        b_seq = full_seq[start-1 : end]
                        
                    write_fasta(rm_id, b_type, species_str, b_stages_str, rm_desc, b_seq)

            for line in f_in:
                if line.strip() == '//':
                    process_record(record_lines)
                    record_lines = []
                else:
                    record_lines.append(line)
                    
            if record_lines:
                process_record(record_lines)

    def run(self):
        """Execute extraction, conversion, and cleaning process"""
        print(f"\n[Module: DfamFetcher] Processing TaxID: {self.tax_id}")
        print(f"Using Python: {self.famdb_python}")
        print(f"Using Internal FamDB: {self.famdb_path}")

        cmd_famdb = f"{self.famdb_python} {self.famdb_path} -i {self.dfam_db} families -f embl --curated -ad {self.tax_id} > {self.temp_embl}"
        
        run_command(cmd_famdb, shell=True, description=f"Extracting EMBL data from Dfam")

        if not os.path.exists(self.temp_embl) or os.path.getsize(self.temp_embl) == 0:
            raise FileNotFoundError(f"famdb failed to generate valid data. Please check if TaxID {self.tax_id} exists in your Dfam database.")

        print(f"==> Converting EMBL to FASTA (Internal Python Parser)...")
        self._embl_to_fasta(self.temp_embl, self.temp_raw_fa)

        awk_filter = r"""awk 'BEGIN{RS=">";ORS=""} NR>1 { if ($1 !~ /#Unknown/) { print ">"$0 } }'"""
        
        cmd_clean = f"{awk_filter} {self.temp_raw_fa} > {self.output_fasta}"
        run_command(cmd_clean, shell=True, description="Removing Unknown sequences from Dfam")

        if os.path.exists(self.temp_embl): os.remove(self.temp_embl)
        if os.path.exists(self.temp_raw_fa): os.remove(self.temp_raw_fa)

        print(f"Dfam data extraction completed: {self.output_fasta}")
        return self.output_fasta
