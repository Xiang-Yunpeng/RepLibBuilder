import os
import re


def format_te_label(norm_class, norm_superfam):
    """
    Build the canonical RepeatMasker-style classification label from a
    normalized (class, superfamily) pair.

    This is the SINGLE SOURCE OF TRUTH for label formatting, shared by every
    source module (Repbase / Dfam / HiTE / EDTA / RM2) so that all libraries
    emit an identical label format:

        - superfamily known  ->  "Class/Superfamily"
        - only class known   ->  "Class/Unknown"
        - neither known      ->  "Unknown"
    """
    if norm_superfam and norm_superfam.lower() != "unknown":
        return f"{norm_class}/{norm_superfam}"
    if norm_class and norm_class.lower() != "unknown":
        return f"{norm_class}/Unknown"
    return "Unknown"


class TEDictionary:
    def __init__(self, dict_path):
        """
        Initializes the TE Dictionary with rules specific to different software.
        """
        self.rules = {'EDTA': [], 'RM2': [], 'HiTE': [], 'TEsorter': [], 'Repbase': [], 'Dfam': [], 'Global': []}
        self.load_dictionary(dict_path)

    def load_dictionary(self, dict_path):
        """
        Parses the mapping dictionary file (e.g., te_mapping_rules.tsv) and compiles regex patterns.
        """
        if not os.path.exists(dict_path):
            raise FileNotFoundError(f"Mapping dictionary not found: {dict_path}")
            
        with open(dict_path, 'r') as f:
            for line in f:
                if line.strip() == "" or line.startswith('#'):
                    continue
                parts = line.strip().split('\t')
                if len(parts) >= 4:
                    software, pattern, target_class, target_superfamily = parts[0], parts[1], parts[2], parts[3]
                    if software in self.rules:
                        self.rules[software].append({
                            'pattern': re.compile(pattern),
                            'class': target_class,
                            'superfamily': target_superfamily
                        })

    def normalize(self, software, raw_id):
        """
        Normalizes the raw TE classification into standard class and superfamily based on rules.
        Incorporates a dynamic fallback mechanism to preserve unmapped valid labels.
        """
        # 1. Check software-specific rules
        for rule in self.rules.get(software, []):
            match = rule['pattern'].match(raw_id)
            if match:
                t_class = match.expand(rule['class']) if '\\' in rule['class'] else rule['class']
                t_superfamily = match.expand(rule['superfamily']) if '\\' in rule['superfamily'] else rule['superfamily']
                return t_class, t_superfamily
        
        # 2. Check global rules
        for rule in self.rules['Global']:
            match = rule['pattern'].match(raw_id)
            if match:
                t_class = match.expand(rule['class']) if '\\' in rule['class'] else rule['class']
                t_superfamily = match.expand(rule['superfamily']) if '\\' in rule['superfamily'] else rule['superfamily']
                return t_class, t_superfamily
        
        # 3. Dynamic Fallback Mechanism (Protection Shield)
        # If not found in dictionary, attempt to parse and preserve the raw format 
        # instead of blindly converting to "Unknown".
        raw_id_clean = raw_id.strip()
        if not raw_id_clean or raw_id_clean.lower() == "unknown":
            return "Unknown", "Unknown"
            
        if '/' in raw_id_clean:
            # e.g., "LTR/Gypsy" -> class="LTR", superfamily="Gypsy".
            # Coerce empty segments to "Unknown" so a malformed label with a
            # trailing or leading slash (e.g. "LINE/" produced from an empty
            # TEsorter Superfamily column) never yields an empty class or
            # superfamily. This guarantees the invariant that normalize() never
            # returns an empty string, which keeps downstream label formatting
            # and the Track A "no-refutation" scoring (Unknown -> +2 Asymmetric)
            # correct without any special-casing of "".
            parts = raw_id_clean.split('/', 1)
            t_class = parts[0].strip() or "Unknown"
            t_superfamily = parts[1].strip() or "Unknown"
            return t_class, t_superfamily
        else:
            # e.g., "LINE" -> class="LINE", superfamily="Unknown"
            return raw_id_clean, "Unknown"


def process_and_normalize_fasta(input_fasta, output_fasta, software_name, te_dict):
    """
    Cleans FASTA headers, calculates sequence lengths, and normalizes TE
    classification.

    The header written to `output_fasta` is rewritten into the STANDARDIZED
    ">ID#Class/Superfamily" form *before* the file is handed to TEsorter, so
    that the de novo libraries enter evaluation already normalized and the
    clean FASTA never carries a raw, non-standard label.

    The original (pre-normalization) label is still retained in seq_info as
    'raw_class', because the downstream Track B heuristics depend on it
    (e.g. detecting "MITE" in the source's own label).

    Returns an in-memory dictionary containing all required sequence info.
    """
    seq_info = {}
    current_id = ""
    current_length = 0
    current_raw_class = "Unknown"
    current_norm_class = "Unknown"
    current_norm_superfam = "Unknown"

    with open(input_fasta, 'r') as fin, open(output_fasta, 'w') as fout:
        for line in fin:
            if line.startswith(">"):
                # Persist the previous record before parsing the new header
                if current_id:
                    seq_info[current_id] = {
                        'length': current_length,
                        'raw_class': current_raw_class,
                        'norm_class': current_norm_class,
                        'norm_superfam': current_norm_superfam
                    }

                # Take the first whitespace-delimited token for every source:
                # this is the canonical "ID#Class/Superfamily" part and uniformly
                # drops trailing descriptions (e.g. RM2's RepeatScout annotation).
                clean_header = line.strip().split()[0]

                # Extract ID and raw class from the cleaned header
                header_no_arrow = clean_header[1:]
                if '#' in header_no_arrow:
                    current_id, current_raw_class = header_no_arrow.split('#', 1)
                else:
                    current_id = header_no_arrow
                    current_raw_class = "Unknown"

                # Normalize through the dictionary and write a STANDARDIZED header
                current_norm_class, current_norm_superfam = te_dict.normalize(
                    software_name, current_raw_class
                )
                std_label = format_te_label(current_norm_class, current_norm_superfam)
                fout.write(f">{current_id}#{std_label}\n")

                current_length = 0
            else:
                fout.write(line)
                current_length += len(line.strip())

        # Handle the very last sequence in the file
        if current_id:
            seq_info[current_id] = {
                'length': current_length,
                'raw_class': current_raw_class,
                'norm_class': current_norm_class,
                'norm_superfam': current_norm_superfam
            }

    return seq_info
