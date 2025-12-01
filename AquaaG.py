import subprocess
import pandas as pd
import re
import os
import shutil
import sys
import time
import yaml
import argparse
from Bio import Entrez
from termcolor import colored
from tqdm import tqdm
import logging
import glob

# --- LOGGING SETUP ---
class StepFilter(logging.Filter):
    def filter(self, record):
        return record.msg.startswith(("[Step", "Processing species", "Found species"))

def setup_logging(output_dir):
    logger = logging.getLogger(__name__)
    if logger.hasHandlers():
        logger.handlers.clear()
    logger.setLevel(logging.INFO)
    log_file = os.path.join(output_dir, "pipeline_main.log")
    file_handler = logging.FileHandler(log_file, mode='w')
    file_handler.setLevel(logging.INFO)
    file_formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    file_handler.setFormatter(file_formatter)
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(logging.INFO)
    console_handler.addFilter(StepFilter())
    logger.addHandler(file_handler)
    logger.addHandler(console_handler)
    return logger

# --- CONFIGURATION & PARAMETER HANDLING ---
def load_config(config_file):
    try:
        with open(config_file, 'r') as file:
            config = yaml.safe_load(file)
        required_fields = ['email', 'organism_type', 'quast_threads', 'busco_lineage']
        missing = [f for f in required_fields if f not in config]
        if missing:
            raise ValueError(f"Missing required fields in config: {', '.join(missing)}")
        org_type = config['organism_type']
        if org_type not in ['PK', 'EK']:
            raise ValueError("organism_type must be 'PK' or 'EK'")
        if org_type == 'PK' and 'prokka_kingdom' not in config:
            raise ValueError("Missing 'prokka_kingdom' for PK pipeline.")
        if org_type == 'EK':
            if 'braker_params' not in config:
                raise ValueError("Missing 'braker_params' section for EK pipeline.")
            if 'cpus' not in config.get('braker_params', {}):
                 raise ValueError("Missing required 'cpus' field in braker_params")
        config.setdefault('quast_params', {})
        config.setdefault('prokka_params', {})
        config.setdefault('busco_params', {})
        return config
    except Exception as e:
        print(f"FATAL ERROR loading configuration file '{config_file}': {e}", file=sys.stderr)
        raise

# --- CORE UTILITY FUNCTIONS ---
def format_params(params_dict):
    params = []
    for key, value in params_dict.items():
        param = f"--{key.replace('_', '-')}"
        if isinstance(value, bool) and value:
            params.append(param)
        elif not isinstance(value, bool):
            params.append(f"{param} {value}")
    return " ".join(params)

def run_command_logged(command, logger, cwd=None, error_message="Command failed"):
    logger.info(f"Executing command in '{cwd or os.getcwd()}'\n$ {command}")
    try:
        result = subprocess.run(command, shell=True, check=True, text=True, capture_output=True, cwd=cwd, executable="/bin/bash")
        if result.stdout: logger.info(f"STDOUT:\n{result.stdout.strip()}")
        if result.stderr: logger.warning(f"STDERR:\n{result.stderr.strip()}")
    except subprocess.CalledProcessError as e:
        logger.error(f"FATAL ERROR: {error_message}\n--- FAILED COMMAND ---\n{command}\n--- STDOUT ---\n{e.stdout}\n--- STDERR ---\n{e.stderr}")
        sys.exit(1)

def fetch_reference_urls(species_name, logger):
    logger.info(f"Searching NCBI for reference files for '{species_name}'...")
    search_queries = [ f'"{species_name}"[Organism] AND "reference genome"[Filter] AND "latest refseq"[Filter]', f'"{species_name}"[Organism] AND "representative genome"[Filter] AND "latest refseq"[Filter]', f'"{species_name}"[Organism] AND latest[Filter]' ]
    for i, term in enumerate(search_queries):
        try:
            handle = Entrez.esearch(db="assembly", term=term, retmax=5, sort="relevance")
            record = Entrez.read(handle); handle.close()
            if not record["IdList"]: continue
            summary_handle = Entrez.esummary(db="assembly", id=",".join(record["IdList"]), report="full")
            summaries = Entrez.read(summary_handle); summary_handle.close()
            for summary in summaries["DocumentSummarySet"]["DocumentSummary"]:
                ftp_path = summary.get("FtpPath_RefSeq") or summary.get("FtpPath_GenBank")
                if not ftp_path: continue
                ftp_path = ftp_path.replace("ftp://", "https://")
                base_name = os.path.basename(ftp_path)
                fasta_url = f"{ftp_path}/{base_name}_genomic.fna.gz"
                gff_url = f"{ftp_path}/{base_name}_genomic.gff.gz"
                protein_url = f"{ftp_path}/{base_name}_protein.faa.gz"
                logger.info(f"Success! Found reference: {summary['AssemblyAccession']}")
                return fasta_url, gff_url, protein_url
        except Exception as e:
            logger.error(f"Entrez search failed: {e}. Waiting before retry...")
            time.sleep(5)
            continue
    logger.error(f"Could not find a suitable reference genome for '{species_name}'.")
    return None, None, None

def download_file_and_unzip(url, output_dir, logger, max_retries=3):
    if not url: return None
    filename = os.path.join(output_dir, os.path.basename(url))
    final_filename = filename.removesuffix('.gz')
    if os.path.exists(final_filename):
        logger.info(f"File already exists: {final_filename}")
        return final_filename
    for attempt in range(max_retries):
        try:
            run_command_logged(f"wget --show-progress -O {filename} {url}", logger, error_message=f"Failed to download {url}")
            run_command_logged(f"gunzip -f {filename}", logger, error_message=f"Failed to decompress {filename}")
            return final_filename
        except SystemExit:
            logger.warning(f"Download attempt {attempt + 1} failed.")
            if attempt + 1 < max_retries: time.sleep(10)
            else: logger.error(f"All download attempts failed for {url}.")
    return None

def fetch_assembly_data(organism_query, logger):
    logger.info(f"Fetching assembly metadata for: '{organism_query}'")
    cleaned_query = re.sub(r'\s*\(.*\)\s*', '', organism_query).strip()
    logger.info(f"Using cleaned query for shell command: '{cleaned_query}'")

    cols = "AssemblyAccession,SpeciesName,AssemblyName,AssemblyStatus,SubmitterOrganization,TotalLength,ContigN50"
    filtered_query = f'{cleaned_query} AND "latest"[Filter]'
    command = f'esearch -db assembly -query "{filtered_query}" | esummary | xtract -pattern DocumentSummary -element {cols}'
    try:
        result = subprocess.run(command, shell=True, text=True, capture_output=True, check=True, timeout=300)
        if not result.stdout.strip():
            logger.warning("No data returned from NCBI for the query.")
            return pd.DataFrame()
        cols_list = cols.split(',')
        num_cols = len(cols_list)
        raw_data = [line.split("\t") for line in result.stdout.strip().split("\n")]
        padded_data = []
        for i, row in enumerate(raw_data):
            if len(row) != num_cols:
                logger.warning(f"Record {i+1} has missing data (expected {num_cols}, got {len(row)}). Padding with 'N/A'.")
                while len(row) < num_cols:
                    row.append("N/A")
            padded_data.append(row)
        df = pd.DataFrame(padded_data, columns=cols_list)
        logger.info(f"Fetched metadata for {len(df)} assembly records.")
        return df
    except subprocess.TimeoutExpired:
        logger.error(f"NCBI query for '{organism_query}' timed out after 5 minutes and was cancelled.")
        return pd.DataFrame()
    except subprocess.CalledProcessError as e:
        logger.error(f"NCBI query failed: {e.stderr}")
        return pd.DataFrame()

def filter_for_india(df, logger):
    indian_cities = [ "IND", "Indian", "india", "India", "Agartala", "Ahmedabad", "Aizawl", "Ajmer", "Allahabad", "Amritsar", "Anand", "Avikanagar", "Aurangabad", "Amravati", "Bangalore", "Bareilly", "Batinda", "Belgavi", "Banglore", "Bengaluru", "Bhopal", "Berhampur", "Bhagalpur", "Bhilai", "Bhubaneswar", "Bhubaneshwar", "Bilaspur", "Bibinagar", "Calicut", "Chandigarh", "Chennai", "Cochin", "Dehradun", "Deoghar", "Delhi", "Dindigul", "Dhanbad", "Durgapur", "Faridabad", "Gangtok", "Gandhinagar", "Goa", "Greater Noida", "Gorakhpur", "Gurugram", "Guwahati", "Hyderabad", "Telangana", "Imphal", "Indore", "Itanagar", "Jaipur", "Jabalpur", "Jammu", "Jamshedpur", "Jalandhar", "Jodhpur", "Jorhat", "Kalaburagi", "Kalpakkam", "Kerala", "Kalyani", "Kanpur", "Karnal", "Kangra", "Kannur", "Kashipur", "Kochi", "Kolkata", "Kolkata, Hindupur", "Leh", "Lucknow", "Malda", "Manesar", "Manipur", "Mangalagiri", "Mandi", "Mumbai", "Mysuru", "Nagpur", "Nadiad", "Nellore", "Navi Mumbai", "New Delhi", "Noida", "Patiala", "Pilani", "Puducherry", "Pune", "Prayagraj", "Bareli", "Raichur", "Raipur", "Rajendranagar", "Rajkot", "Rishikesh", "Rohtak", "Roorkee", "Ropar", "Rourkela", "Sambalpur", "Secunderabad", "Sikkim", "Shibpur", "Shillong", "Shimla", "Silchar", "Sonepat", "Sri City", "Srinagar", "Surat", "Tadepalligudem", "Tezpur", "Thiruvananthapuram", "Thiruvananthapura", "Tiruchirappalli", "Tirupati", "Trichy", "Tuljapur", "Udaipur", "Una", "Vadodara", "Varanasi", "Vijayawada", "Visakhapatnam", "Warangal", "Yupia", "RCB", "ICMR", "IMTECH", "CSIR", "SPPU", "IIT", "NIT", "IISC", "IIIT", "Central of University", "IGIB", "ICGEB", "CDRI", "SRM", "AIIMS" ]
    logger.info("Filtering for assemblies submitted from India.")
    if df.empty: return df
    pattern = re.compile(r"\b(" + "|".join(map(re.escape, indian_cities)) + r")\b", re.IGNORECASE)
    filtered_df = df[df["SubmitterOrganization"].str.contains(pattern, na=False)]
    logger.info(f"Found {len(filtered_df)} assemblies matching the India filter.")
    return filtered_df

def download_assemblies(filtered_df, output_dir, group, num_assemblies, logger):
    successful_accessions = []
    assembly_ids = filtered_df["AssemblyAccession"].tolist()
    total_to_download = len(assembly_ids) if num_assemblies is None else min(num_assemblies, len(assembly_ids))
    if not assembly_ids: return successful_accessions
    with tqdm(total=total_to_download, desc=f"Downloading {total_to_download} assemblies") as pbar:
        for accession in assembly_ids[:total_to_download]:
            source = "refseq" if accession.startswith("GCF") else "genbank"
            command = f"ncbi-genome-download -s {source} -F fasta -o {output_dir} -p {os.cpu_count()} all --assembly-accessions {accession}"
            try:
                run_command_logged(command, logger, error_message=f"Download failed for {accession}")
                successful_accessions.append(accession)
                pbar.update(1)
            except SystemExit:
                logger.error(f"Skipping accession {accession} due to download failure")
    return successful_accessions

def decompress_and_rename_assemblies(assembly_dir, accessions, logger):
    for accession in accessions:
        pattern = os.path.join(assembly_dir, f"**/{accession}*.fna.gz")
        files = glob.glob(pattern, recursive=True)
        if not files:
            logger.warning(f"No FASTA file found for {accession}")
            continue
        for file in files:
            run_command_logged(f"gunzip -f {file}", logger, error_message=f"Failed to decompress {file}")
            decompressed_file = file.removesuffix('.gz')
            new_file = os.path.join(assembly_dir, f"{accession}.fna")
            shutil.move(decompressed_file, new_file)

# --- NEW FUNCTION TO FIX FASTA HEADERS ---
def fix_fasta_headers(assembly_dir, accessions, logger):
    """
    Simplifies FASTA headers to only the first word (e.g., >CP12345.1)
    This prevents header mismatch errors in downstream tools like GeneMark.
    """
    logger.info("Simplifying FASTA headers to first word only.")
    for accession in accessions:
        fasta_file = os.path.join(assembly_dir, f"{accession}.fna")
        if not os.path.exists(fasta_file):
            logger.warning(f"Could not find FASTA file {fasta_file} to fix headers.")
            continue
            
        temp_file = f"{fasta_file}.tmp"
        
        try:
            with open(fasta_file, 'r') as infile, open(temp_file, 'w') as outfile:
                for line in infile:
                    if line.startswith(">"):
                        # Split the header line by space, take only the first part
                        outfile.write(line.split()[0] + "\n")
                    else:
                        outfile.write(line)
            
            # Replace the original file with the fixed file
            shutil.move(temp_file, fasta_file)
            
        except Exception as e:
            logger.error(f"Failed to fix headers for {fasta_file}: {e}")
            if os.path.exists(temp_file):
                os.remove(temp_file) # Clean up
            sys.exit(1)

# --- PIPELINE STAGE FUNCTIONS ---
def run_quast(accessions, assembly_dir, quast_dir, ref_genome, gff_genome, threads, quast_params, logger):
    assembly_files = [os.path.join(assembly_dir, f"{acc}.fna") for acc in accessions]
    quast_cmd = f"quast.py {' '.join(assembly_files)} -r {ref_genome} -g {gff_genome} -t {threads} {format_params(quast_params)} -o {quast_dir}"
    run_command_logged(quast_cmd, logger, error_message="QUAST failed.")

def run_prokka(accessions, assembly_dir, prokka_dir, kingdom, prokka_params, logger):
    for accession in accessions:
        input_fasta = os.path.join(assembly_dir, f"{accession}.fna")
        output_prefix = os.path.join(prokka_dir, accession)
        prokka_cmd = f"prokka --kingdom {kingdom} --outdir {output_prefix} --prefix {accession} {format_params(prokka_params)} {input_fasta}"
        run_command_logged(prokka_cmd, logger, error_message=f"Prokka failed for {accession}")

def run_repeat_modeling(accession_dir, genome_path, cpus, logger):
    logger.info("Stage 6.1: Building RepeatModeler database.")
    run_command_logged(f"BuildDatabase -name genome_db \"{genome_path}\"", logger, cwd=accession_dir, error_message="RepeatModeler BuildDatabase failed.")
    run_command_logged(f"RepeatModeler -database genome_db -threads {cpus}", logger, cwd=accession_dir, error_message="RepeatModeler failed.")
    rm_dirs = [d for d in os.listdir(accession_dir) if d.startswith("RM_") and os.path.isdir(os.path.join(accession_dir, d))]
    if not rm_dirs: logger.error("FATAL: RepeatModeler did not create an output directory."); sys.exit(1)
    latest_rm_dir = max(rm_dirs, key=lambda d: os.path.getmtime(os.path.join(accession_dir, d)))
    repeat_lib_path = os.path.abspath(os.path.join(accession_dir, latest_rm_dir, "consensi.fa.classified"))
    logger.info("Stage 6.2: Running RepeatMasker to mask the genome.")
    repeatmasker_cmd = f"RepeatMasker -engine ncbi -lib \"{repeat_lib_path}\" -pa {cpus} -gff -xsmall \"{genome_path}\""
    run_command_logged(repeatmasker_cmd, logger, cwd=accession_dir, error_message="RepeatMasker failed.")
    masked_genome_path = f"{genome_path}.masked"
    if not os.path.exists(masked_genome_path): logger.error(f"FATAL: Masked genome not created at {masked_genome_path}"); sys.exit(1)
    logger.info(f"Successfully created masked genome: {masked_genome_path}")
    return masked_genome_path


import os
import sys
import shutil
import subprocess
from termcolor import colored

# Helper functions assumed to be defined elsewhere in your script:
# - run_repeat_modeling(path, genome, cpus, logger)
# - format_params(dict)
# - run_command_logged(cmd, logger, error_message)

import os
import sys
import shutil
import subprocess

from termcolor import colored

# assumes you have these helpers somewhere in your codebase:
# - run_repeat_modeling(...)
# - format_params(...)
# - run_command_logged(...)


def run_braker_for_accession(accession, config, assembly_dir, braker_output_dir,
                             protein_evidence_path, logger):
    """
    Run BRAKER3 inside the official Docker image with robust handling for missing ProtHint outputs.

    - EP mode: protein only (ProtHint or precomputed hints)
    - ET mode: RNA-Seq only (SRA / FASTQ / BAM / hints)
    - ETP mode: protein + RNA-Seq (native BRAKER3 / GeneMark-ETP mode)

    IMPORTANT:
      * GeneMark-ETP is used from inside the container (/opt/ETP).
      * We do NOT override GENEMARK_PATH or mount a host GeneMark-ETP install.
      * gm_key is mounted from a host path resolved in this priority:
          1) BRAKER_GM_KEY env var (if set and file exists)
          2) project_root/.braker_keys/gm_key
          3) ~/.braker_keys/gm_key

    If BRAKER/GeneMark-ETP fails with messages about missing
    prothint/prothint.gff or evidence.gff, this function will:

      1. Search for any prothint.gff / evidence.gff under the working tree.
      2. Copy them into:
         - GeneMark-ETP/prothint/
         - GeneMark-ETP/proteins.fa/prothint/         (if proteins.fa dir exists)
         - GeneMark-ETP/proteins.fa/penalty/prothint/ (if penalty dir exists)
      3. Re-run BRAKER with the SAME arguments (no extra --hints in ETP mode).
    """
    logger.info(colored(f"--- Starting Eukaryotic Annotation with BRAKER3 (Docker) for {accession} ---", "green"))

    DOCKER_IMAGE = os.environ.get("BRAKER_DOCKER_IMAGE", "teambraker/braker3")
    DOCKER_TAG   = os.environ.get("BRAKER_DOCKER_TAG", "latest")

    # Paths & params
    original_genome_path = os.path.abspath(os.path.join(assembly_dir, f"{accession}.fna"))
    accession_out_dir = os.path.abspath(os.path.join(braker_output_dir, accession))
    os.makedirs(accession_out_dir, exist_ok=True)

    braker_params = config.get('braker_params', {})
    cpus = int(braker_params.get('cpus', 1))

    # --- Repeat modeling & masking on host ---
    masked_genome_path = run_repeat_modeling(
        os.path.dirname(original_genome_path),
        original_genome_path,
        cpus,
        logger
    )

    # --- Mode selection (interactive) ---
    print(colored("\n" + "=" * 50, "cyan"))
    print(colored(f"BRAKER3 MODE SELECTION for {accession}", "cyan", attrs=["bold"]))
    print("Please choose the BRAKER mode to run:")
    print("  [1] EP Mode (Protein evidence only)")
    print("  [2] ET Mode (Transcriptome evidence only)")
    print("  [3] ETP Mode (Protein AND Transcriptome evidence)")
    mode_choice = ""
    while mode_choice not in ["1", "2", "3"]:
        mode_choice = input(colored("Enter your choice (1, 2, or 3): ", "white", attrs=["bold"]))

    # Common BRAKER args (inside container)
    braker_args = [
        "braker.pl",
        f"--genome={masked_genome_path}",
        f"--species={accession}",
        f"--workingdir={accession_out_dir}",
        f"--threads={cpus}",
        "--softmasking",
    ]

    # Extra params from config (except ones we handle explicitly)
    extra_params_dict = {k: v for k, v in braker_params.items() if k not in ['cpus', 'protein_evidence']}
    if extra_params_dict:
        braker_args.append(format_params(extra_params_dict))

    # -----------------------
    #  EP MODE (protein only)
    # -----------------------
    if mode_choice == "1":
        logger.info("Configuring EP Mode (Protein only)...")

        default_protein = protein_evidence_path or braker_params.get('protein_evidence')
        print("\n" + "-" * 50)
        print(colored("Protein evidence for EP mode:", "cyan"))

        if default_protein:
            default_protein = os.path.abspath(default_protein)
            print(f"  [1] Use protein file from config: {default_protein}")
            print("  [2] Enter a different protein FASTA")
            choice = ""
            while choice not in ["1", "2"]:
                choice = input(colored("Choose protein input (1/2): ", "white", attrs=["bold"]))
            if choice == "1":
                protein_path_to_use = default_protein
            else:
                protein_path_to_use = ""
        else:
            print("  No protein_evidence in config; you must provide a protein FASTA.")
            protein_path_to_use = ""

        while not protein_path_to_use or not os.path.exists(protein_path_to_use):
            protein_path_to_use = input(colored("Enter FULL path to your protein.fa: ", "white"))
            if not os.path.exists(protein_path_to_use):
                logger.error(f"File not found: {protein_path_to_use}")

        protein_path_to_use = os.path.abspath(protein_path_to_use)
        logger.info(f"Using EP protein evidence file: {protein_path_to_use}")

        print(colored(f"Using protein FASTA: {protein_path_to_use}", "yellow"))
        print("-" * 50)

        print(colored("Protein Evidence Input:", "cyan"))
        print("  [1] Let ProtHint run inside BRAKER (recommended for EP).")
        print("  [2] Use a pre-computed ProtHint hints file (.gff) instead of running ProtHint.")
        prot_choice = ""
        while prot_choice not in ["1", "2"]:
            prot_choice = input(colored("Enter your choice (1 or 2): ", "white", attrs=["bold"]))

        if prot_choice == "1":
            braker_args.append(f"--prot_seq={protein_path_to_use}")
        else:
            hints_file = ""
            while not os.path.exists(hints_file):
                hints_file = input(colored("Enter the FULL path to your prothint_augustus.gff: ", "white"))
                if not os.path.exists(hints_file):
                    logger.error(f"File not found: {hints_file}")
            braker_args.append(f"--hints={os.path.abspath(hints_file)}")

    # ----------------------------
    #  ET MODE (transcript only)
    # ----------------------------
    elif mode_choice == "2":
        logger.info("Configuring ET Mode (Transcriptome only)...")
        print("\n" + "-" * 50)
        print(colored("RNA-Seq Data Input Method:", "cyan"))
        print("  [1] SRA IDs (auto-download)")
        print("  [2] Local FASTQ files")
        print("  [3] Local BAM files")
        print("  [4] Pre-computed hints file (.gff)")
        rna_choice = ""
        while rna_choice not in ["1", "2", "3", "4"]:
            rna_choice = input(colored("Enter your choice (1-4): ", "white", attrs=["bold"]))

        if rna_choice == "1":
            sra_ids = input(colored("Enter SRA IDs (comma-separated): ", "white"))
            braker_args.append(f"--rnaseq_sets_ids={sra_ids}")

        elif rna_choice == "2":
            sra_ids = input(colored("Enter library IDs (e.g., SRA_ID1,SRA_ID2): ", "white"))
            fastq_dir = ""
            while not os.path.isdir(fastq_dir):
                fastq_dir = input(colored("Enter path to directory containing FASTQs: ", "white"))
                if not os.path.isdir(fastq_dir):
                    logger.error(f"Directory not found: {fastq_dir}")
            braker_args.append(f"--rnaseq_sets_ids={sra_ids}")
            braker_args.append(f"--rnaseq_sets_dirs={os.path.abspath(fastq_dir)}")

        elif rna_choice == "3":
            bam_files_str = input(colored("Enter FULL paths to BAM files (comma-separated): ", "white"))
            bam_files = [os.path.abspath(f.strip()) for f in bam_files_str.split(',')]
            for f in bam_files:
                if not os.path.exists(f):
                    logger.error(f"BAM file not found: {f}")
                    sys.exit(1)
            braker_args.append(f"--bam={','.join(bam_files)}")

        elif rna_choice == "4":
            hints_file = ""
            while not os.path.exists(hints_file):
                hints_file = input(colored("Enter FULL path to RNA-Seq hints.gff: ", "white"))
                if not os.path.exists(hints_file):
                    logger.error(f"File not found: {hints_file}")
            braker_args.append(f"--hints={os.path.abspath(hints_file)}")

    # -----------------------------------
    #  ETP MODE (protein + transcript)
    # -----------------------------------
    else:
        logger.info("Configuring ETP Mode (Protein + Transcriptome)...")

        # ---- ETP protein selection (interactive confirm/override) ----
        default_protein = protein_evidence_path or braker_params.get('protein_evidence')
        print("\n" + "-" * 50)
        print(colored("Protein evidence for ETP mode:", "cyan"))

        if default_protein:
            default_protein = os.path.abspath(default_protein)
            print(f"  [1] Use protein file from config: {default_protein}")
            print("  [2] Enter a different protein FASTA")
            choice = ""
            while choice not in ["1", "2"]:
                choice = input(colored("Choose protein input (1/2): ", "white", attrs=["bold"]))
            if choice == "1":
                protein_path_to_use = default_protein
            else:
                protein_path_to_use = ""
        else:
            print("  No protein_evidence in config; you must provide a protein FASTA.")
            protein_path_to_use = ""

        while not protein_path_to_use or not os.path.exists(protein_path_to_use):
            protein_path_to_use = input(colored("Enter FULL path to your protein.fa: ", "white"))
            if not os.path.exists(protein_path_to_use):
                logger.error(f"File not found: {protein_path_to_use}")

        protein_path_to_use = os.path.abspath(protein_path_to_use)
        logger.info(f"Using ETP protein evidence file: {protein_path_to_use}")
        print(colored(f"Using protein FASTA for ETP: {protein_path_to_use}", "yellow"))
        print("-" * 50)

        # In ETP mode do NOT pass --hints; GeneMark-ETP manages ProtHint internally
        braker_args.append(f"--prot_seq={protein_path_to_use}")

        # ---- RNA input selection for ETP ----
        print(colored("RNA-Seq Data Input Method:", "cyan"))
        print("  [1] SRA IDs (auto-download)")
        print("  [2] Local FASTQ files")
        print("  [3] Local BAM files")
        rna_choice = ""
        while rna_choice not in ["1", "2", "3"]:
            rna_choice = input(colored("Enter your choice (1-3): ", "white", attrs=["bold"]))

        if rna_choice == "1":
            sra_ids = input(colored("Enter SRA IDs (comma-separated): ", "white"))
            braker_args.append(f"--rnaseq_sets_ids={sra_ids}")

        elif rna_choice == "2":
            sra_ids = input(colored("Enter library IDs (e.g., SRA_ID1,SRA_ID2): ", "white"))
            fastq_dir = ""
            while not os.path.isdir(fastq_dir):
                fastq_dir = input(colored("Enter path to directory containing FASTQs: ", "white"))
                if not os.path.isdir(fastq_dir):
                    logger.error(f"Directory not found: {fastq_dir}")
            braker_args.append(f"--rnaseq_sets_ids={sra_ids}")
            braker_args.append(f"--rnaseq_sets_dirs={os.path.abspath(fastq_dir)}")

        elif rna_choice == "3":
            bam_files_str = input(colored("Enter FULL paths to BAM files (comma-separated): ", "white"))
            bam_files = [os.path.abspath(f.strip()) for f in bam_files_str.split(',')]
            for f in bam_files:
                if not os.path.exists(f):
                    logger.error(f"BAM file not found: {f}")
                    sys.exit(1)
            braker_args.append(f"--bam={','.join(bam_files)}")

    # --- Host resources to mount into container ---
    augustus_cfg = os.environ.get("AUGUSTUS_CONFIG_PATH", os.path.expanduser("~/augustus_config"))
    if not os.path.isdir(augustus_cfg):
        logger.error(f"FATAL: AUGUSTUS config dir not found at {augustus_cfg}.")
        sys.exit(1)

    # --- GeneMark-ETP gm_key discovery ---
    # Priority:
    #   1) BRAKER_GM_KEY env var
    #   2) project-local .braker_keys/gm_key (created by setup script)
    #   3) ~/.braker_keys/gm_key (legacy / manual setup)
    gm_key_env = os.environ.get("BRAKER_GM_KEY", "").strip()

    # project_root = directory where this Python file lives
    project_root = os.path.dirname(os.path.abspath(__file__))
    project_gm_key = os.path.join(project_root, ".braker_keys", "gm_key")
    home_gm_key = os.path.expanduser("~/.braker_keys/gm_key")

    if gm_key_env and os.path.isfile(gm_key_env):
        gm_key_host = gm_key_env
    elif os.path.isfile(project_gm_key):
        gm_key_host = project_gm_key
    elif os.path.isfile(home_gm_key):
        gm_key_host = home_gm_key
    else:
        logger.error(
            "FATAL: GeneMark-ETP license key (gm_key) not found.\n"
            f"  Checked:\n"
            f"    BRAKER_GM_KEY={gm_key_env or '<not set>'}\n"
            f"    project-local: {project_gm_key}\n"
            f"    home fallback: {home_gm_key}\n"
            "  Please run the setup script again or set BRAKER_GM_KEY explicitly."
        )
        sys.exit(1)

    # Build docker volumes: mount everything we reference
    mounts = set()

    def add_mount_for(path: str):
        d = os.path.abspath(path if os.path.isdir(path) else os.path.dirname(path))
        if d:
            mounts.add(d)

    add_mount_for(masked_genome_path)
    add_mount_for(accession_out_dir)
    add_mount_for(augustus_cfg)

    # Also mount evidence paths if present in args
    for a in list(braker_args):
        if (
            a.startswith("--prot_seq=")
            or a.startswith("--bam=")
            or a.startswith("--hints=")
            or a.startswith("--rnaseq_sets_dirs=")
        ):
            val = a.split("=", 1)[1]
            for piece in val.split(","):
                add_mount_for(piece)

    # Compose docker -v flags
    vflags = " ".join([f'-v "{m}:{m}"' for m in sorted(mounts)])

    # Ensure expected ProtHint dirs exist inside GeneMark-ETP working dir
    expected_prothint_dir = os.path.join(accession_out_dir, "GeneMark-ETP", "prothint")
    gm_proteins_prothint_dir = os.path.join(accession_out_dir, "GeneMark-ETP", "proteins.fa", "prothint")
    os.makedirs(expected_prothint_dir, exist_ok=True)

    # Build base docker run command
    base_docker_cmd = (
        f'docker run --rm '
        f'-u {os.getuid()}:{os.getgid()} '
        f'{vflags} '
        f'-v "{gm_key_host}:/root/.gm_key:ro" '
        f'-e AUGUSTUS_CONFIG_PATH="{augustus_cfg}" '
        f'-w "{accession_out_dir}" '
        f'{DOCKER_IMAGE}:{DOCKER_TAG}'
    )

    docker_cmd = base_docker_cmd + " " + " ".join(braker_args)

    # --- Initial BRAKER run ---
    logger.info("Running BRAKER3 in Docker (initial run).")
    logger.info(f"Command: {docker_cmd}")
    proc = subprocess.run(
        docker_cmd,
        shell=True,
        text=True,
        capture_output=True,
        executable="/bin/bash"
    )

    if proc.returncode == 0:
        if proc.stdout:
            logger.info(f"STDOUT:\n{proc.stdout.strip()}")
        if proc.stderr:
            logger.warning(f"STDERR:\n{proc.stderr.strip()}")
        logger.info(colored(f"--- BRAKER3 (Docker) Annotation for {accession} COMPLETED ---", "green"))
        return

    # Log failure
    logger.error(f"Initial BRAKER Docker run failed for {accession} (exit {proc.returncode}).")
    if proc.stdout:
        logger.error(f"STDOUT:\n{proc.stdout.strip()}")
    if proc.stderr:
        logger.error(f"STDERR:\n{proc.stderr.strip()}")

    # --- Heuristic: check BRAKER stderr + GeneMark-ETP.stderr for ProtHint issues ---
    stderr_lower = (proc.stderr or "").lower()

    gm_etp_stderr_path = os.path.join(accession_out_dir, "errors", "GeneMark-ETP.stderr")
    gm_etp_stderr_lower = ""
    if os.path.exists(gm_etp_stderr_path):
        try:
            with open(gm_etp_stderr_path, "r") as fh:
                gm_etp_stderr_lower = fh.read().lower()
        except Exception as e:
            logger.warning(f"Could not read GeneMark-ETP.stderr for recovery heuristic: {e}")

    combined_err = stderr_lower + "\n" + gm_etp_stderr_lower

    if (
        "prothint/prothint.gff" in combined_err
        or "prothint/evidence.gff" in combined_err
        or "option --f1 prothint" in combined_err
    ):
        logger.info("Detected missing ProtHint outputs in BRAKER/GeneMark-ETP failure. Attempting automatic recovery...")

        # search for any prothint.gff / evidence.gff under the accession working tree
        candidates = []
        for root, dirs, files in os.walk(accession_out_dir):
            for f in files:
                if f.lower().endswith("prothint.gff") or f.lower().endswith("evidence.gff"):
                    candidates.append(os.path.join(root, f))

        # Also search some reasonable sibling output directories (older runs)
        extra_search_dirs = [
            os.path.join(os.path.dirname(accession_out_dir), ".."),
            os.path.join(os.path.dirname(accession_out_dir), "..", "..", "braker_output"),
            os.path.join(os.getcwd(), "braker"),
        ]
        for d in extra_search_dirs:
            if os.path.isdir(d):
                for root, dirs, files in os.walk(d):
                    for f in files:
                        if f.lower().endswith("prothint.gff") or f.lower().endswith("evidence.gff"):
                            candidates.append(os.path.join(root, f))

        candidates = sorted(set(candidates))
        if candidates:
            logger.info(
                f"Found {len(candidates)} candidate ProtHint/evidence files. "
                f"Using first match: {candidates[0]}"
            )
            found_prothint = None
            found_evidence = None

            # prefer exact prothint.gff if present
            for c in candidates:
                if os.path.basename(c).lower().endswith("prothint.gff"):
                    found_prothint = c
                    break
            if not found_prothint:
                found_prothint = candidates[0]

            # attempt to find evidence.gff in same dir
            maybe_dir = os.path.dirname(found_prothint)
            possible_evidence = os.path.join(maybe_dir, "evidence.gff")
            if os.path.exists(possible_evidence):
                found_evidence = possible_evidence
            else:
                for c in candidates:
                    if os.path.basename(c).lower().endswith("evidence.gff"):
                        found_evidence = c
                        break

            try:
                # 1) GeneMark-ETP/prothint
                os.makedirs(expected_prothint_dir, exist_ok=True)
                new_prothint = os.path.join(expected_prothint_dir, "prothint.gff")
                shutil.copy(found_prothint, new_prothint)
                logger.info(f"Copied ProtHint file to expected location: {new_prothint}")

                if found_evidence:
                    new_evidence = os.path.join(expected_prothint_dir, "evidence.gff")
                    shutil.copy(found_evidence, new_evidence)
                    logger.info(f"Copied evidence.gff to expected location: {new_evidence}")

                # 2) GeneMark-ETP/proteins.fa/prothint (if proteins.fa dir exists)
                gm_proteins_dir = os.path.join(accession_out_dir, "GeneMark-ETP", "proteins.fa")
                if os.path.isdir(gm_proteins_dir):
                    os.makedirs(gm_proteins_prothint_dir, exist_ok=True)
                    prot2 = os.path.join(gm_proteins_prothint_dir, "prothint.gff")
                    shutil.copy(found_prothint, prot2)
                    logger.info(f"Copied ProtHint file to proteins.fa/prothint: {prot2}")

                    if found_evidence:
                        evid2 = os.path.join(gm_proteins_prothint_dir, "evidence.gff")
                        shutil.copy(found_evidence, evid2)
                        logger.info(f"Copied evidence.gff to proteins.fa/prothint: {evid2}")
                else:
                    logger.info("GeneMark-ETP/proteins.fa directory not present yet; skipping mirror copy for now.")

                # 3) GeneMark-ETP/proteins.fa/penalty/prothint (for the gmes_petap.pl penalty runs)
                gm_penalty_dir = os.path.join(
                    accession_out_dir,
                    "GeneMark-ETP",
                    "proteins.fa",
                    "penalty"
                )
                if os.path.isdir(gm_penalty_dir):
                    penalty_prothint_dir = os.path.join(gm_penalty_dir, "prothint")
                    os.makedirs(penalty_prothint_dir, exist_ok=True)

                    prot3 = os.path.join(penalty_prothint_dir, "prothint.gff")
                    shutil.copy(found_prothint, prot3)
                    logger.info(f"Copied ProtHint file to penalty/prothint: {prot3}")

                    if found_evidence:
                        evid3 = os.path.join(penalty_prothint_dir, "evidence.gff")
                        shutil.copy(found_evidence, evid3)
                        logger.info(f"Copied evidence.gff to penalty/prothint: {evid3}")

                    logger.info(
                        "Contents of GeneMark-ETP/proteins.fa/penalty/prothint after recovery: "
                        f"{os.listdir(penalty_prothint_dir)}"
                    )
                else:
                    logger.info(
                        "GeneMark-ETP/proteins.fa/penalty does not exist yet; "
                        "cannot mirror ProtHint there."
                    )

                logger.info(
                    "Contents of GeneMark-ETP/prothint after recovery: "
                    f"{os.listdir(expected_prothint_dir)}"
                )
            except Exception as e:
                logger.error(f"Failed to copy ProtHint/evidence files during recovery: {e}")

            # Re-run BRAKER with the SAME arguments (no extra --hints in ETP mode).
            docker_cmd2 = base_docker_cmd + " " + " ".join(braker_args)
            logger.info("Re-running BRAKER with the same arguments after repairing ProtHint locations.")
            logger.info(f"Command: {docker_cmd2}")
            proc2 = subprocess.run(
                docker_cmd2,
                shell=True,
                text=True,
                capture_output=True,
                executable="/bin/bash"
            )

            if proc2.returncode == 0:
                if proc2.stdout:
                    logger.info(f"STDOUT:\n{proc2.stdout.strip()}")
                if proc2.stderr:
                    logger.warning(f"STDERR:\n{proc2.stderr.strip()}")
                logger.info(colored(
                    f"--- BRAKER3 (Docker) Annotation for {accession} COMPLETED (after recovery) ---",
                    "green"
                ))
                return
            else:
                logger.error(f"Recovery attempt failed (exit {proc2.returncode}).")
                if proc2.stdout:
                    logger.error(f"STDOUT:\n{proc2.stdout.strip()}")
                if proc2.stderr:
                    logger.error(f"STDERR:\n{proc2.stderr.strip()}")

    # If we get here, recovery did not succeed or it was not a ProtHint-related failure.
    logger.error(f"FATAL: BRAKER3 (Docker) failed for {accession} and automatic recovery did not succeed.")
    run_command_logged(docker_cmd, logger, error_message=f"BRAKER3 (Docker) failed for {accession}")

def run_busco(accessions, config, assembly_dir, prokka_dir, braker_dir, busco_dir, logger):
    logger.info(colored("[Step 7/7] Running BUSCO for annotation quality assessment...", "green"))
    os.makedirs(busco_dir, exist_ok=True)
    lineage = config['busco_lineage']
    busco_params = config.get('busco_params', {})
    org_type = config['organism_type']
    
    for accession in accessions:
        if org_type == 'PK':
            input_file = os.path.join(prokka_dir, accession, f"{accession}.faa")
        else:
            input_file = os.path.join(braker_dir, accession, "braker.aa")
            
        if not os.path.exists(input_file):
            logger.warning(f"No input file for BUSCO at {input_file}. Skipping {accession}.")
            continue
            
        output_prefix = os.path.join(busco_dir, accession)
        busco_cmd = f"busco -i {input_file} -o {output_prefix} -m protein -l {lineage} {format_params(busco_params)} --cpu {config.get('quast_threads', 1)} --force"
        run_command_logged(busco_cmd, logger, error_message=f"BUSCO failed for {accession}")
        
        summary_file = glob.glob(os.path.join(busco_dir, accession, f"short_summary.specific.{lineage}.*.txt"))
        if summary_file:
            with open(summary_file[0], 'r') as f:
                logger.info(f"BUSCO summary for {accession}:\n{f.read().strip()}")
        else:
            logger.warning(f"No BUSCO summary file generated for {accession}.")

# --- WORKFLOW FUNCTIONS ---
def get_species_for_kingdom(kingdom, num_species, filter_india, logger):
    logger.info(f"Searching for {num_species} species in kingdom '{kingdom}' with required assemblies...")
    species_to_process = []
    seen_species = set()
    retstart = 0
    batch_size = 200
    while len(species_to_process) < num_species:
        handle = Entrez.esearch(db="assembly", term=f'"{kingdom}"[Organism] AND "latest refseq"[Filter]', retmax=batch_size, retstart=retstart)
        record = Entrez.read(handle); handle.close()
        if not record["IdList"]:
            logger.warning("No more species found in the kingdom.")
            break
        retstart += batch_size
        summary_handle = Entrez.esummary(db="assembly", id=",".join(record["IdList"]), report="full")
        summaries = Entrez.read(summary_handle); summary_handle.close()
        potential_species = [s["SpeciesName"] for s in summaries["DocumentSummarySet"]["DocumentSummary"] if s["SpeciesName"] and s["SpeciesName"] not in seen_species]
        seen_species.update(potential_species)
        for species_name in potential_species:
            df = fetch_assembly_data(f'"{species_name}"[Organism]', logger)
            if df.empty: continue
            if filter_india:
                if not filter_for_india(df, logger).empty:
                    logger.info(colored(f"Found species with Indian assemblies: {species_name}", "cyan"))
                    species_to_process.append(species_name)
            else:
                species_to_process.append(species_name)
            if len(species_to_process) >= num_species:
                break
        if len(species_to_process) >= num_species:
            break
    if not species_to_process:
        logger.error(f"Could not find any species in '{kingdom}' that meet the specified criteria.")
    return species_to_process

def process_species(species_name, config, base_output_dir, num_assemblies_per_species, filter_india, is_species_mode, logger):
    species_slug = species_name.replace(" ", "_").replace("/", "_")
    species_output_dir = os.path.join(base_output_dir, species_slug)
    os.makedirs(species_output_dir, exist_ok=True)
    logger.info(colored(f"\n{'='*20} Processing species: {species_name} {'='*20}", "yellow"))

    dirs = { "ref": os.path.join(species_output_dir, "reference_data"), "asm": os.path.join(species_output_dir, "assemblies"), "quast": os.path.join(species_output_dir, "quast_output"), "prokka": os.path.join(species_output_dir, "prokka_output"), "braker": os.path.join(species_output_dir, "braker_output"), "busco": os.path.join(species_output_dir, "busco_output") }
    for d in dirs.values(): os.makedirs(d, exist_ok=True)

    logger.info(colored("[Step 1/7] Finding and downloading reference data...", "green"))
    fasta_url, gff_url, protein_url = fetch_reference_urls(species_name, logger)
    if not fasta_url or not gff_url:
        logger.error(f"Skipping species {species_name}: Could not find a suitable reference genome."); return False
    ref_genome = download_file_and_unzip(fasta_url, dirs["ref"], logger)
    gff_genome = download_file_and_unzip(gff_url, dirs["ref"], logger)
    if not ref_genome or not gff_genome:
        logger.error(f"Skipping species {species_name}: Failed to download reference data."); return False
    
    protein_evidence_path = None
    if config['organism_type'] == 'EK':
        protein_evidence_path = download_file_and_unzip(protein_url, dirs["ref"], logger)
        if not protein_evidence_path:
            logger.warning("Could not auto-download protein evidence. BRAKER may fall back to config file.")

    logger.info(colored("[Step 2/7] Preparing list of assemblies...", "green"))
    df = fetch_assembly_data(f'"{species_name}"[Organism]', logger)
    if df.empty:
        logger.warning(f"No other assemblies found for {species_name}."); return True

    if filter_india:
        filtered_df = filter_for_india(df, logger)
        if filtered_df.empty:
            if is_species_mode:
                logger.error(f"FATAL: No Indian-submitted assemblies found for '{species_name}'. Halting."); sys.exit(1)
            else:
                logger.warning(f"No Indian assemblies for {species_name}. Skipping this species."); return False
    else:
        filtered_df = df

    filtered_df['TotalLength'] = pd.to_numeric(filtered_df['TotalLength'], errors='coerce')
    min_genome_size = 100000 
    initial_count = len(filtered_df)
    filtered_df = filtered_df[filtered_df['TotalLength'] >= min_genome_size]
    logger.info(f"Filtered out {initial_count - len(filtered_df)} assemblies smaller than {min_genome_size} bp.")
    if filtered_df.empty:
        logger.warning(f"No assemblies remained for {species_name} after size filtering. Skipping.")
        return False

    metadata_path = os.path.join(species_output_dir, f"{species_slug}_metadata.tsv")
    filtered_df.to_csv(metadata_path, sep='\t', index=False)
    logger.info(f"Assembly metadata saved to {metadata_path}")
    
    logger.info(colored("[Step 3/7] Downloading assemblies...", "green"))
    accessions = download_assemblies(filtered_df, dirs["asm"], config.get('group', 'all'), num_assemblies_per_species, logger)
    if not accessions: return True

    logger.info(colored("[Step 4/7] Finalizing assembly files...", "green"))
    decompress_and_rename_assemblies(dirs["asm"], accessions, logger)
    
    # --- FIX: Simplify FASTA headers ---
    fix_fasta_headers(dirs["asm"], accessions, logger)
    # --- END FIX ---

    logger.info(colored("[Step 5/7] Running QUAST...", "green"))
    run_quast(accessions, dirs["asm"], dirs["quast"], ref_genome, gff_genome, config['quast_threads'], config.get('quast_params', {}), logger)

    if config['organism_type'] == 'EK':
        logger.info(colored("[Step 6/7] Running Eukaryotic Annotation (BRAKER3)...", "green"))
        for accession in accessions:
            run_braker_for_accession(accession, config, dirs["asm"], dirs["braker"], protein_evidence_path, logger)
    else:
        logger.info(colored("[Step 6/7] Running Prokaryotic Annotation (Prokka)...", "green"))
        run_prokka(accessions, dirs["asm"], dirs["prokka"], config.get('prokka_kingdom'), config.get('prokka_params', {}), logger)
    
    logger.info(colored("[Step 7/7] Running BUSCO...", "green"))
    run_busco(accessions, config, dirs["asm"], dirs["prokka"], dirs["braker"], dirs["busco"], logger)
    logger.info(colored(f"--- Finished processing species: {species_name} ---", "yellow"))
    return True

# --- MAIN WORKFLOW ---
def main():
    parser = argparse.ArgumentParser(description="Genomic Assembly and Annotation Pipeline (AquaG)", formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("-c", "--config", required=True, help="Path to the config file (YAML)")
    parser.add_argument("-o", "--output-dir", required=True, help="Main output directory")

    source_group = parser.add_argument_group('Assembly Source (choose one)')
    source_exclusive_group = source_group.add_mutually_exclusive_group(required=True)
    source_exclusive_group.add_argument("-s", "--species", type=str, help="Single species to search and process.")
    source_exclusive_group.add_argument("-k", "--kingdom", type=str, help="Kingdom to search for species to process.")
    source_exclusive_group.add_argument("--assembly-file", type=str, help="Path to a text file of assembly accession IDs.\nRequires a reference organism to be set in the config file.")
    
    filter_group = parser.add_argument_group('Filtering and Limiting')
    filter_group.add_argument("--num-assemblies", type=int, default=1, help="Max number of assemblies to process PER SPECIES. Default: 1")
    filter_group.add_argument("--num-species", type=int, default=1, help="In KINGDOM mode, the number of species to process. Default: 1")
    location_exclusive_group = filter_group.add_mutually_exclusive_group()
    location_exclusive_group.add_argument("-I", "--india-only", action='store_true', help="Process assemblies from Indian submitters only.")
    location_exclusive_group.add_argument("-a", "--any-source", action='store_true', help="Process assemblies from any submitter (default).")

    args = parser.parse_args()
    os.makedirs(args.output_dir, exist_ok=True)
    logger = setup_logging(args.output_dir)
    config = load_config(args.config)
    Entrez.email = config['email']

    error_count = 0
    if args.species:
        if not process_species(args.species, config, args.output_dir, args.num_assemblies, args.india_only, is_species_mode=True, logger=logger):
            error_count += 1
    
    elif args.kingdom:
        species_to_process = get_species_for_kingdom(args.kingdom, args.num_species, args.india_only, logger)
        for species_name in species_to_process:
            try:
                if not process_species(species_name, config, args.output_dir, args.num_assemblies, args.india_only, is_species_mode=False, logger=logger):
                    error_count += 1
            except Exception as e:
                logger.error(f"An unexpected error occurred while processing {species_name}. Skipping. Error: {e}")
                error_count += 1
                continue
    
    elif args.assembly_file:
        if 'organism' not in config or not config.get('organism'):
            logger.error("FATAL: For --assembly-file mode, you must specify a reference 'organism' in the config file."); sys.exit(1)
        
        reference_organism = config['organism']
        logger.info(colored(f"--- Running in Assembly File Mode for reference: {reference_organism} ---", "cyan"))
        
        dirs = { "ref": os.path.join(args.output_dir, "reference_data"), "asm": os.path.join(args.output_dir, "assemblies"), "quast": os.path.join(args.output_dir, "quast_output"), "prokka": os.path.join(args.output_dir, "prokka_output"), "braker": os.path.join(args.output_dir, "braker_output"), "busco": os.path.join(args.output_dir, "busco_output") }
        for d in dirs.values(): os.makedirs(d, exist_ok=True)
        
        fasta_url, gff_url, protein_url = fetch_reference_urls(reference_organism, logger)
        if not fasta_url or not gff_url: logger.error(f"Could not find a reference for '{reference_organism}'."); sys.exit(1)
        ref_genome = download_file_and_unzip(fasta_url, dirs["ref"], logger)
        gff_genome = download_file_and_unzip(gff_url, dirs["ref"], logger)
        if not ref_genome or not gff_genome: logger.error("Failed to download reference data."); sys.exit(1)
        
        protein_evidence_path = None
        if config['organism_type'] == 'EK':
            protein_evidence_path = download_file_and_unzip(protein_url, dirs["ref"], logger)

        try:
            with open(args.assembly_file, 'r') as f:
                accessions = [line.strip() for line in f if line.strip()]
            logger.info(f"Read {len(accessions)} accessions from {args.assembly_file}")
            df = pd.DataFrame({"AssemblyAccession": accessions})
            full_metadata_df = fetch_assembly_data(f'"{reference_organism}"[Organism]', logger)
            if not full_metadata_df.empty:
                subset_metadata = full_metadata_df[full_metadata_df['AssemblyAccession'].isin(accessions)]
                if not subset_metadata.empty:
                    metadata_path = os.path.join(args.output_dir, f"{reference_organism.replace(' ', '_')}_metadata.tsv")
                    subset_metadata.to_csv(metadata_path, sep='\t', index=False)
                    logger.info(f"Assembly metadata saved to {metadata_path}")
                else:
                    logger.warning(f"Could not find metadata for the specified accessions for organism '{reference_organism}'.")

            download_assemblies(df, dirs["asm"], config.get('group', 'all'), None, logger)
            decompress_and_rename_assemblies(dirs["asm"], accessions, logger)
            
            # --- FIX: Simplify FASTA headers ---
            fix_fasta_headers(dirs["asm"], accessions, logger)
            # --- END FIX ---
            
            run_quast(accessions, dirs["asm"], dirs["quast"], ref_genome, gff_genome, config['quast_threads'], config.get('quast_params', {}), logger)
            
            if config['organism_type'] == 'EK':
                for accession in accessions:
                    run_braker_for_accession(accession, config, dirs["asm"], dirs["braker"], protein_evidence_path, logger)
            else:
                run_prokka(accessions, dirs["asm"], dirs["prokka"], config.get('prokka_kingdom'), config.get('prokka_params', {}), logger)

            run_busco(accessions, config, dirs["asm"], dirs["prokka"], dirs["braker"], dirs["busco"], logger)
        except FileNotFoundError:
            logger.error(f"Assembly file not found: {args.assembly_file}"); sys.exit(1)

    if error_count > 0:
        logger.info(colored(f"\nPIPELINE COMPLETED WITH {error_count} ERRORS OR SKIPPED SPECIES. Please check the log for details. ⚠️", "red"))
    else:
        logger.info(colored("\nPIPELINE COMPLETED SUCCESSFULLY! ✅", "cyan"))

if __name__ == "__main__":
    main()

