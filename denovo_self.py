import os
import sys
import subprocess
import gzip
import shutil  # 修正：匯入 shutil 模組

def print_help():
    """
    Print help message for script usage.
    """
    help_message = """
    Usage: python denovo_self.py --indir <input_directory> --outdir <output_directory> [--threads <number_of_threads>]

    Options:
    --indir         Path to the directory containing input files (e.g., trimmed FASTQ files).
    --outdir        Path to the output directory where results will be saved.
    --threads       Number of threads/CPUs to use (default: 1).
    -h, --help      Display this help message.

    Description:
    This script processes FASTQ files for de novo assembly by:
    1. Converting FASTQ files to FASTA format.
    2. Renaming FASTA files and organizing them.
    3. Generating a population map file (popmap.tsv).
    4. Running the `denovo_map.pl` pipeline from Stacks.

    Additional functionality:
    - If output files already exist, the script will prompt the user to overwrite or skip processing.
    """
    print(help_message)

def check_file_and_ask(file_path, description):
    """
    Check if a file exists and ask the user whether to overwrite or skip.

    Args:
        file_path (str): Path to the file.
        description (str): Description of the file for user reference.

    Returns:
        bool: True if the step should proceed (overwrite), False to skip.
    """
    if os.path.exists(file_path):
        while True:
            user_input = input(f"{description} already exists at {file_path}. Do you want to overwrite it? (y/n): ").strip().lower()
            if user_input == 'y':
                print(f"Overwriting {description}...")
                return True
            elif user_input == 'n':
                print(f"Skipping {description}...")
                return False
            else:
                print("Invalid input. Please enter 'y' to overwrite or 'n' to skip.")
    return True  # Proceed if the file does not exist

def convert_and_rename(input_dir, output_dir):
    """
    Convert input files to FASTA format and rename them.

    Args:
        input_dir (str): Directory containing input files.
        output_dir (str): Directory to save the converted files.
    """
    output_fasta_dir = os.path.join(output_dir, "fasta")

    # Check if the entire directory exists
    if os.path.exists(output_fasta_dir):
        while True:
            user_input = input(f"The directory {output_fasta_dir} already exists. Do you want to overwrite it? (y/n): ").strip().lower()
            if user_input == 'y':
                print(f"Overwriting the directory {output_fasta_dir}...")
                shutil.rmtree(output_fasta_dir)  # Remove the existing directory
                os.makedirs(output_fasta_dir, exist_ok=True)
                break
            elif user_input == 'n':
                print(f"Skipping FASTA conversion because the directory {output_fasta_dir} already exists.")
                return  # Skip conversion
            else:
                print("Invalid input. Please enter 'y' to overwrite or 'n' to skip.")
    else:
        os.makedirs(output_fasta_dir, exist_ok=True)

    # Process R1 and R2 files
    paired_files = {}
    for file_name in os.listdir(input_dir):
        if file_name.endswith(".fastq.gz"):  # Check for gzipped FASTQ files
            if "_R1_trimmed" in file_name:
                sample_name = handle_sample_naming(file_name, "_R1_trimmed")
                paired_files.setdefault(sample_name, [None, None])[0] = file_name
            elif "_R2_trimmed" in file_name:
                sample_name = handle_sample_naming(file_name, "_R2_trimmed")
                paired_files.setdefault(sample_name, [None, None])[1] = file_name

    # Convert and rename
    for sample_name, (r1_file, r2_file) in paired_files.items():
        if r1_file:
            convert_to_fasta(input_dir, r1_file, output_fasta_dir, sample_name, "1")
        if r2_file:
            convert_to_fasta(input_dir, r2_file, output_fasta_dir, sample_name, "2")

def handle_sample_naming(file_name, suffix):
    """
    Handle naming correction for any prefix, ensuring single-digit IDs are zero-padded.

    Args:
        file_name (str): Original file name.
        suffix (str): Suffix to remove for naming.

    Returns:
        tuple: Original sample name and corrected sample name.
    """
    # 定義所有可能的樣點名稱前綴
    prefixes = ["WHX", "BT", "MK", "TI", "DF", "KT"]
    # 移除檔案名中的後綴，取得樣本名稱
    sample_name = file_name.split(suffix)[0]

    for prefix in prefixes:
        if sample_name.startswith(prefix):
            # 判斷是否需要補零
            numeric_part = sample_name[len(prefix):]
            if numeric_part.isdigit() and len(numeric_part) == 1:  # 單位數需補零
                corrected_name = f"{prefix}0{numeric_part}"
                print(f"Original: {sample_name}, Corrected: {corrected_name}")
                return corrected_name
            print(f"Original: {sample_name}, Corrected: {sample_name}")
            return sample_name  # 無需補零的直接返回
    raise ValueError(f"Unexpected sample prefix in file name: {file_name}")


def convert_to_fasta(input_dir, file_name, output_fasta_dir, sample_name, read_type):
    """
    Convert a single FASTQ file to FASTA format.

    Args:
        input_dir (str): Directory containing input files.
        file_name (str): Name of the file to convert.
        output_fasta_dir (str): Directory to save the converted file.
        sample_name (str): Sample name for renaming.
        read_type (str): "1" for R1, "2" for R2.
    """
    input_file = os.path.join(input_dir, file_name)
    fasta_file = os.path.join(output_fasta_dir, f"{sample_name}.{read_type}.fasta")
    print(f"Converting {file_name} to {fasta_file}...")
    try:
        with gzip.open(input_file, 'rt') as f_in, open(fasta_file, 'w') as f_out:
            for i, line in enumerate(f_in):
                if i % 4 == 0:  # FASTQ header line
                    f_out.write(f">{line[1:]}")  # Convert @ to >
                elif i % 4 == 1:  # Sequence line
                    f_out.write(line)
        print(f"Converted {file_name} to {fasta_file}")
    except Exception as e:
        print(f"Error converting {file_name}: {e}")

def generate_popmap(fasta_dir, output_path):
    """
    Generate a population map file.

    Args:
        fasta_dir (str): Directory containing FASTA files.
        output_path (str): Path to save the population map.
    """
    prefixes = ["WHX", "BT", "MK", "TI", "DF", "KT"]
    if check_file_and_ask(output_path, "Population map"):
        print("Generating population map...")
        samples = []
        for file_name in os.listdir(fasta_dir):
            if file_name.endswith(".fasta"):
                sample_name = file_name.split('.')[0]
                samples.append(sample_name)
        # 按順序排列樣本
        sorted_samples = sorted(samples, key=lambda x: prefixes.index(x[:3]))
        with open(output_path, 'w') as popmap:
            for sample in sorted_samples:
                popmap.write(f"{sample}\tpop1\n")
        print(f"Population map generated at {output_path}")

def run_denovo_map(output_dir, threads, exe_path, bin_path):
    """
    Run the `denovo_map.pl` pipeline.

    Args:
        output_dir (str): Directory containing the pipeline output.
        threads (int): Number of threads/CPUs to use.
        exe_path (str): Path to the `denovo_map.pl` executable.
        bin_path (str): Path to the directory containing Stacks executables.
    """
    log_file = os.path.join(output_dir, "denovo_map.log")
    if check_file_and_ask(log_file, "denovo_map.log"):
        print("Running denovo_map.pl...")
        try:
            subprocess.run([
                exe_path,
                "-M", "3", "-n", "2", "-T", str(threads),
                "--samples", os.path.join(output_dir, "fasta"),
                "--popmap", os.path.join(output_dir, "popmap.tsv"),
                "--out-path", os.path.join(output_dir, "stacks_output"),
                "--paired",
                "-e", bin_path
            ], check=True)
            print("denovo_map.pl completed successfully.")
        except subprocess.CalledProcessError as e:
            print(f"Error running denovo_map.pl: {e}")
            sys.exit(1)

def main():
    """
    Main function to parse arguments and execute the pipeline.
    """
    import argparse
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("--indir", required=True, help="Path to the input directory.")
    parser.add_argument("--outdir", required=True, help="Path to the output directory.")
    parser.add_argument("--threads", type=int, default=1, help="Number of threads to use (default: 1).")
    parser.add_argument("--help", "-h", action="store_true", help="Display help message.")
    args = parser.parse_args()

    if args.help:
        print_help()
        sys.exit(0)

    input_dir = args.indir
    output_dir = args.outdir
    threads = args.threads
    exe_path = "/Users/weiwei19510gmail.com/Documents/bioinformation/stacks-2.68/scripts/denovo_map.pl"
    bin_path = "/Users/weiwei19510gmail.com/Documents/bioinformation/stacks-2.68/"

    if not os.path.exists(input_dir):
        print(f"Error: Input directory '{input_dir}' does not exist.")
        sys.exit(1)
    os.makedirs(output_dir, exist_ok=True)

    # Step 1: Convert and rename files
    convert_and_rename(input_dir, output_dir)

    # Step 2: Generate popmap.tsv
    popmap_path = os.path.join(output_dir, "popmap.tsv")
    generate_popmap(os.path.join(output_dir, "fasta"), popmap_path)

    # Ask if user wants to proceed with denovo_map.pl
    user_input = input("Do you want to proceed with running denovo_map.pl? (y/n): ").strip().lower()
    if user_input == 'y':
        run_denovo_map(output_dir, threads, exe_path, bin_path)
    else:
        print("Skipping denovo_map.pl execution.")

if __name__ == "__main__":
    main()
