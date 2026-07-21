import subprocess
import sys
import os
import shutil

def check_dependency(tool_name):
    """
    Check if a tool exists in the system PATH.
    If not found, print error and exit.
    """
    if shutil.which(tool_name) is None:
        print(f"[Error] Dependency not found: {tool_name}")
        print(f"Please ensure it is installed and added to your PATH.")
        sys.exit(1)

def check_file(filepath, description):
    """
    Check if a file exists.
    If not found, print error and exit.
    Returns absolute path.
    """
    if not filepath:
        print(f"[Error] Path not provided for {description}.")
        sys.exit(1)
        
    if not os.path.isfile(filepath):
        print(f"[Error] File not found for {description}: {filepath}")
        sys.exit(1)
    return os.path.abspath(filepath)

def ensure_dir(directory):
    """
    Ensure directory exists, create if missing.
    Returns absolute path.
    """
    if not os.path.exists(directory):
        try:
            os.makedirs(directory)
        except OSError as e:
            print(f"[Error] Cannot create directory {directory}: {e}")
            sys.exit(1)
    return os.path.abspath(directory)

def run_command(cmd, shell=False, description=None):
    """
    Generic wrapper for subprocess.run.
    
    Args:
        cmd (list or str): Command line list (preferred) or string.
        shell (bool): Run via shell (required for pipes/wildcards).
        description (str): Description to print.
    """
    if description:
        print(f"==> {description}...")
    
    try:
        # Use executable='/bin/bash' for shell=True to support process substitution if needed
        subprocess.run(
            cmd, 
            shell=shell, 
            check=True, 
            executable='/bin/bash' if shell else None
        )
    except subprocess.CalledProcessError as e:
        print(f"\n[Fatal Error] Command failed (Exit Code: {e.returncode})")
        print(f"Command: {cmd if isinstance(cmd, str) else ' '.join(cmd)}")
        sys.exit(1)
    except Exception as e:
        print(f"\n[Fatal Error] Unexpected error: {e}")
        sys.exit(1)

def append_normalized_fasta(input_path, output_handle, id_prefix=None):
    """
    Read input_path, normalize sequence (Upper case + 60bp wrapping),
    and append to output_handle.
    
    Args:
        input_path (str): Path to input FASTA.
        output_handle (file object): Opened file object for writing.
        id_prefix (str): Optional prefix to add to sequence IDs (e.g., 'repbase').
                         Result: >prefix_originalID
    """
    if not input_path or not os.path.exists(input_path):
        print(f"[Warning] File missing or empty path, skipping: {input_path}")
        return

    try:
        with open(input_path, 'r') as f_in:
            header = None
            seq_buffer = []
            
            for line in f_in:
                line = line.strip()
                if not line: continue
                
                if line.startswith(">"):
                    # Write previous sequence
                    if header:
                        full_seq = "".join(seq_buffer).upper()
                        output_handle.write(f"{header}\n")
                        for i in range(0, len(full_seq), 60):
                            output_handle.write(full_seq[i:i+60] + "\n")
                    
                    # Process new header
                    # Remove '>' and any trailing whitespace
                    raw_id_line = line[1:].strip()
                    
                    if id_prefix:
                        # Add prefix: >prefix_ID
                        # Note: We keep the rest of the description line if present
                        header = f">{id_prefix}_{raw_id_line}"
                    else:
                        header = line
                    
                    seq_buffer = []
                else:
                    seq_buffer.append(line)
            
            # Write final sequence
            if header and seq_buffer:
                full_seq = "".join(seq_buffer).upper()
                output_handle.write(f"{header}\n")
                for i in range(0, len(full_seq), 60):
                    output_handle.write(full_seq[i:i+60] + "\n")
                    
    except Exception as e:
        print(f"[Error] Failed processing file {input_path}: {e}")
        sys.exit(1)
