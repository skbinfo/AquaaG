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
    """Custom filter for console to show only step messages."""
    def filter(self, record):
        return record.msg.startswith("[Step")

def setup_logging(output_dir):
    """Set up file and console logging."""
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.INFO)
    log_file = os.path.join(output_dir, "pipeline.log")
    file_handler = logging.FileHandler(log_file, mode='w')
    file_handler.setLevel(logging.INFO)
    file_formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    file_handler.setFormatter(file_formatter)
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(logging.INFO)
    console_handler.addFilter(StepFilter())
    logger.handlers.clear()
    logger.addHandler(file_handler)
    logger.addHandler(console_handler)
    return logger

# --- CONFIGURATION & PARAMETER HANDLING ---
def load_config(config_file):
    """Load and validate configuration from YAML file."""
    try:
        with open(config_file, 'r') as file:
            config = yaml.safe_load(file)

        required_fields = ['organism_query', 'email', 'organism', 'reference_genome_url', 'gff_url', 'group', 'organism_type', 'quast_threads']
        missing = [f for f in required_fields if f not in config]
        if missing:
            raise ValueError(f"Missing required fields in config: {', '.join(missing)}")

        org_type = config['organism_type']
        if org_type not in ['PK', 'EK']:
            raise ValueError("organism_type must be 'PK' or 'EK'")

        if org_type == 'PK' and 'prokka_kingdom' not in config:
            raise ValueError("Missing 'prokka_kingdom' for PK pipeline.")
        if org_type == 'EK':
            if 'maker_params' not in config:
                raise ValueError("Missing 'maker_params' section for EK pipeline.")
            # maker_exe_ctl is no longer required from the user as it's auto-generated
            required_maker_params = ['protein_evidence', 'est_evidence', 'maker_opts_ctl', 'maker_bopts_ctl', 'cpus']
            missing_maker = [f for f in required_maker_params if f not in config.get('maker_params', {})]
            if missing_maker:
                raise ValueError(f"Missing required fields in maker_params: {', '.join(missing_maker)}")

        config.setdefault('quast_params', {})
        config.setdefault('prokka_params', {})
        return config
    except Exception as e:
        print(f"FATAL ERROR loading configuration file '{config_file}': {e}", file=sys.stderr)
        raise

# --- UTILITY AND CORE FUNCTIONS ---

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

def download_reference_genome(url, output_dir, logger):
    filename = os.path.join(output_dir, os.path.basename(url))
    final_filename = filename.removesuffix('.gz')
    if not os.path.exists(final_filename):
        run_command_logged(f"wget -O {filename} {url}", logger, error_message=f"Failed to download {url}")
        run_command_logged(f"gunzip -f {filename}", logger, error_message=f"Failed to decompress {filename}")
    return final_filename

def read_assembly_file(assembly_file, logger):
    accession_pattern = re.compile(r"GC[AF]_\d{9}\.\d")
    try:
        with open(assembly_file, 'r') as f:
            valid_accessions = [line.strip() for line in f if accession_pattern.fullmatch(line.strip())]
        if not valid_accessions: raise ValueError(f"No valid accession IDs found in {assembly_file}")
        logger.info(f"Read {len(valid_accessions)} valid accession IDs from {assembly_file}")
        return valid_accessions
    except FileNotFoundError:
        logger.error(f"Assembly file {assembly_file} not found"); raise

def fetch_assembly_data(organism_query, logger):
    """Fetch assembly data from NCBI using Entrez."""
    command = f'esearch -db assembly -query "{organism_query}" | esummary | xtract -pattern DocumentSummary -element AssemblyAccession,AssemblyName,AssemblyStatus,SubmitterOrganization'
    logger.info(f"Fetching assembly data for '{organism_query}'")
    result = subprocess.run(command, shell=True, text=True, capture_output=True, check=True)
    
    if not result.stdout.strip():
        raise ValueError("No data returned from NCBI for the query.")
        
    data = [line.split("\t") for line in result.stdout.strip().split("\n")]
    df = pd.DataFrame(data, columns=["AssemblyAccession", "AssemblyName", "AssemblyStatus", "SubmitterOrganization"])
    logger.info(f"Fetched {len(df)} assembly records.")
    return df

def filter_rows_by_cities(df, location, logger):
    """Filter rows by location without limiting the number of assemblies."""
    # This function is retained for its specific India-related keywords
    indian_cities = [
        "IND", "Indian", "india", "India", "Agartala", "Ahmedabad", "Aizawl", "Ajmer", "Allahabad", "Amritsar",
        "Anand", "Avikanagar", "Aurangabad", "Amravati", "Bangalore", "Bareilly", "Batinda", "Belgavi", "Banglore", "Bengaluru",
        "Bhopal", "Berhampur", "Bhagalpur", "Bhilai", "Bhubaneswar", "Bhubaneshwar", "Bilaspur", "Bibinagar", "Calicut",
        "Chandigarh", "Chennai", "Cochin", "Dehradun", "Deoghar", "Delhi", "Dindigul", "Dhanbad",
        "Durgapur", "Faridabad", "Gangtok", "Gandhinagar", "Goa", "Greater Noida", "Gorakhpur", "Gurugram", "Guwahati",
        "Hyderabad", "Telangana", "Imphal", "Indore", "Itanagar", "Jaipur", "Jabalpur", "Jammu", "Jamshedpur",
        "Jalandhar", "Jodhpur", "Jorhat", "Kalaburagi", "Kalpakkam", "Kerala", "Kalyani", "Kanpur", "Karnal", "Kangra", "Kannur",
        "Kashipur", "Kochi", "Kolkata", "Kolkata, Hindupur", "Leh", "Lucknow", "Malda", "Manesar",
        "Manipur", "Mangalagiri", "Mandi", "Mumbai", "Mysuru", "Nagpur", "Nadiad", "Nellore", "Navi Mumbai", "New Delhi",
        "Noida", "Patiala", "Pilani", "Puducherry", "Pune", "Prayagraj", "Bareli", "Raichur", "Raipur", "Rajendranagar",
        "Rajkot", "Rishikesh", "Rohtak", "Roorkee", "Ropar", "Rourkela", "Sambalpur", "Secunderabad", "Sikkim", "Shibpur",
        "Shillong", "Shimla", "Silchar", "Sonepat", "Sri City", "Srinagar", "Surat", "Tadepalligudem", "Tezpur",
        "Thiruvananthapuram", "Thiruvananthapura", "Tiruchirappalli", "Tirupati", "Trichy", "Tuljapur", "Udaipur", "Una",
        "Vadodara", "Varanasi", "Vijayawada", "Visakhapatnam", "Warangal", "Yupia",
        "RCB", "ICMR", "IMTECH", "CSIR", "SPPU", "IIT", "NIT", "IISC", "IIIT", "Central of University", "IGIB", "ICGEB", "CDRI", "SRM", "AIIMS"
    ]
    
    logger.info(f"Filtering {len(df)} assemblies for location: {location}")
    if location.lower() == "india":
        pattern = re.compile(r"\b(" + "|".join(map(re.escape, indian_cities)) + r")\b", re.IGNORECASE)
        filtered_df = df[df["SubmitterOrganization"].str.contains(pattern, na=False)]
    else:
        # Fallback for other locations if ever needed, though default is India
        pattern = re.escape(location)
        filtered_df = df[df["SubmitterOrganization"].str.contains(pattern, case=False, na=False)]
    
    logger.info(f"Found {len(filtered_df)} assemblies matching the location '{location}'.")
    if filtered_df.empty:
        logger.warning("No assemblies matched the location filter.")
    return filtered_df

def download_assemblies(filtered_df, output_dir, group, num_assemblies, logger):
    successful_accessions = []
    assembly_ids = filtered_df["AssemblyAccession"].tolist()
    total_to_download = len(assembly_ids) if num_assemblies is None else min(num_assemblies, len(assembly_ids))
    with tqdm(total=total_to_download, desc="Downloading assemblies") as pbar:
        for accession in assembly_ids[:total_to_download]:
            source = "refseq" if accession.startswith("GCF") else "genbank"
            command = f"ncbi-genome-download -s {source} -F fasta -o {output_dir} -p {os.cpu_count()} {group} --assembly-accessions {accession}"
            try:
                run_command_logged(command, logger, error_message=f"Download failed for {accession}")
                successful_accessions.append(accession)
                pbar.update(1)
            except SystemExit: logger.error(f"Skipping accession {accession} due to download failure.")
    return successful_accessions

def decompress_and_rename_assemblies(output_dir, successful_accessions, logger):
    for accession in tqdm(successful_accessions, desc="Processing files"):
        find_command = f"find {output_dir} -name '*{accession}*.fna.gz'"
        result = subprocess.run(find_command, shell=True, check=True, capture_output=True, text=True)
        gz_file = result.stdout.strip().split('\n')[0]
        if gz_file:
            run_command_logged(f"gunzip -f '{gz_file}'", logger)
            fna_file = gz_file.removesuffix('.gz')
            shutil.move(fna_file, os.path.join(output_dir, f"{accession}.fna"))
    run_command_logged(f"find {output_dir} -mindepth 1 -type d -empty -delete", logger)

def run_quast(accessions, assembly_dir, quast_dir, ref_genome, ref_gff, threads, params, logger):
    for accession in tqdm(accessions, desc="Running QUAST"):
        cmd = f"quast.py '{os.path.join(assembly_dir, f'{accession}.fna')}' -r '{ref_genome}' -g '{ref_gff}' -o '{os.path.join(quast_dir, accession)}' --threads {threads} {format_params(params)}"
        run_command_logged(cmd, logger, error_message=f"QUAST failed for {accession}")

def run_prokka(accessions, assembly_dir, prokka_dir, kingdom, params, logger):
    for accession in tqdm(accessions, desc="Running Prokka"):
        cmd = f"prokka --outdir '{os.path.join(prokka_dir, accession)}' --prefix {accession} --kingdom {kingdom} {format_params(params)} '{os.path.join(assembly_dir, f'{accession}.fna')}'"
        run_command_logged(cmd, logger, error_message=f"Prokka failed for {accession}")

# --- EUKARYOTIC ANNOTATION WORKFLOW ---

def run_repeat_modeling(accession_dir, genome_path, cpus, logger):
    logger.info("Stage 6.1: Running RepeatModeler to build de novo repeat library.")
    run_command_logged(f"BuildDatabase -name genome_db '{genome_path}'", logger, cwd=accession_dir, error_message="RepeatModeler BuildDatabase failed.")
    run_command_logged(f"RepeatModeler -database genome_db -threads {cpus}", logger, cwd=accession_dir, error_message="RepeatModeler failed.")
    
    rm_dirs = [d for d in os.listdir(accession_dir) if d.startswith("RM_") and os.path.isdir(os.path.join(accession_dir, d))]
    if not rm_dirs: logger.error("FATAL: RepeatModeler did not create an output directory."); sys.exit(1)
    latest_rm_dir = max(rm_dirs, key=lambda d: os.path.getmtime(os.path.join(accession_dir, d)))
    repeat_lib_path = os.path.abspath(os.path.join(accession_dir, latest_rm_dir, "consensi.fa.classified"))

    logger.info("Stage 6.2: Running RepeatMasker to mask the genome.")
    repeatmasker_cmd = f"RepeatMasker -engine ncbi -lib '{repeat_lib_path}' -pa {cpus} -gff -xsmall '{genome_path}'"
    run_command_logged(repeatmasker_cmd, logger, cwd=accession_dir, error_message="RepeatMasker failed.")
    
    masked_genome_path = f"{genome_path}.masked"
    if not os.path.exists(masked_genome_path): logger.error(f"FATAL: Masked genome not created."); sys.exit(1)
    logger.info(f"Successfully created masked genome: {masked_genome_path}")
    return repeat_lib_path, masked_genome_path

def train_snap_model(round1_output_dir, logger):
    """Train a SNAP HMM model from MAKER round 1."""
    logger.info("Stage 6.4: Training SNAP HMM model.")
    abs_r1_dir = os.path.abspath(round1_output_dir)
    snap_dir = os.path.join(abs_r1_dir, "snap_training")
    os.makedirs(snap_dir, exist_ok=True)
    
    # THE FIX IS HERE: Correct the search path to look inside the *.maker.output directory
    datastore_glob = os.path.join(abs_r1_dir, "*.maker.output", "*_master_datastore_index.log")

    logger.info(f"Searching for MAKER datastore log using pattern: {datastore_glob}")
    datastore_files = glob.glob(datastore_glob)
    if not datastore_files:
        logger.error(f"FATAL: MAKER Round 1 output not found. No datastore log file matching the pattern was found. Check MAKER Round 1 logs for errors."); sys.exit(1)
        
    # We only need to merge the one datastore, not use the glob pattern in the command
    datastore_to_merge = datastore_files[0]
    run_command_logged(f"gff3_merge -s -d '{datastore_to_merge}' > round1.all.gff", logger, cwd=snap_dir, error_message="gff3_merge failed.")
    
    for cmd in ["maker2zff -l 50 -x 0.5 round1.all.gff", "fathom genome.ann genome.dna -categorize 1000", "fathom uni.ann uni.dna -export 1000 -plus", "forge export.ann export.dna", "hmm-assembler.pl snap_model . > snap_model.hmm"]:
        run_command_logged(cmd, logger, cwd=snap_dir, error_message=f"SNAP training failed during: {cmd.split()[0]}")
        
    return os.path.abspath(os.path.join(snap_dir, "snap_model.hmm"))

def create_maker_ctl_files(run_dir, maker_params, logger):
    """Programmatically creates maker_exe.ctl and copies other control files."""
    # Create maker_exe.ctl from scratch to guarantee it's correct
    maker_exe_path = os.path.join(run_dir, 'maker_exe.ctl')
    logger.info(f"Programmatically creating {maker_exe_path}")
    repeatmasker_exe = shutil.which('RepeatMasker')
    if not repeatmasker_exe:
        logger.error("FATAL: Cannot find 'RepeatMasker' executable in the environment's PATH via shutil.which()."); sys.exit(1)
    with open(maker_exe_path, 'w') as f:
        f.write(f"RepeatMasker={repeatmasker_exe}\n")

    # Copy bopts file from template
    shutil.copy(maker_params['maker_bopts_ctl'], run_dir)

def run_maker_for_accession(accession, config, assembly_dir, maker_output_dir, logger):
    """Run the full eukaryotic annotation pipeline for a single accession."""
    logger.info(f"--- Starting Eukaryotic Annotation for {accession} ---")
    original_genome_path = os.path.abspath(os.path.join(assembly_dir, f"{accession}.fna"))
    accession_out_dir = os.path.join(maker_output_dir, accession)
    
    # THE FIX: Create the output directory for this accession before using it.
    os.makedirs(accession_out_dir, exist_ok=True)
    
    maker_params = config['maker_params']
    cpus = maker_params['cpus']

    repeat_lib_path, masked_genome_path = run_repeat_modeling(accession_out_dir, original_genome_path, cpus, logger)
    
    # --- Stage 3: MAKER Round 1 ---
    logger.info(f"Stage 6.3: Running MAKER Round 1 for {accession}.")
    r1_out_dir = os.path.join(accession_out_dir, "round1")
    os.makedirs(r1_out_dir, exist_ok=True)
    create_maker_ctl_files(r1_out_dir, maker_params, logger)
    
    r1_opts_path = os.path.join(r1_out_dir, 'maker_opts.ctl')
    shutil.copy(maker_params['maker_opts_ctl'], r1_opts_path)
    with open(r1_opts_path, 'a') as f:
        f.write(f"\ngenome={os.path.abspath(masked_genome_path)}\nprotein={os.path.abspath(maker_params['protein_evidence'])}\nest={os.path.abspath(maker_params['est_evidence'])}\nrmlib={os.path.abspath(repeat_lib_path)}\ncpus={cpus}\n")

    run_command_logged(f"maker -base {accession}_r1 maker_opts.ctl maker_bopts.ctl maker_exe.ctl", logger, cwd=r1_out_dir, error_message="MAKER Round 1 failed.")
    
    # --- Stage 4: Train SNAP ---
    snaphmm_path = train_snap_model(r1_out_dir, logger)
    
    # --- Stage 5: MAKER Round 2 ---
    logger.info(f"Stage 6.5: Running MAKER Round 2 for {accession}.")
    r2_out_dir = os.path.join(accession_out_dir, "round2")
    os.makedirs(r2_out_dir, exist_ok=True)
    create_maker_ctl_files(r2_out_dir, maker_params, logger)
    
    r2_opts_path = os.path.join(r2_out_dir, 'maker_opts.ctl')
    shutil.copy(maker_params['maker_opts_ctl'], r2_opts_path)
    with open(r2_opts_path, 'a') as f:
        f.write(f"\ngenome={os.path.abspath(masked_genome_path)}\nsnaphmm={os.path.abspath(snaphmm_path)}\nest2genome=0\nprotein2genome=0\ncpus={cpus}\n")

    run_command_logged(f"maker -base {accession}_r2 maker_opts.ctl maker_bopts.ctl maker_exe.ctl", logger, cwd=r2_out_dir, error_message="MAKER Round 2 failed.")
    logger.info(f"--- Eukaryotic Annotation for {accession} COMPLETED ---")

# --- MAIN WORKFLOW ---
def main():
    parser = argparse.ArgumentParser(description="Genomic Assembly and Annotation Pipeline")
    parser.add_argument("-c", "--config", required=True, help="Path to the config file (YAML)")
    parser.add_argument("-o", "--output-dir", required=True, help="Main output directory")
    parser.add_argument("--num-assemblies", type=int, help="Max number of assemblies to process")
    parser.add_argument("--assembly-file", type=str, help="Path to a text file of assembly accession IDs. Overrides default search.")
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    logger = setup_logging(args.output_dir)
    config = load_config(args.config)
    Entrez.email = config['email']

    dirs = { "ref": os.path.join(args.output_dir, "reference_data"), "asm": os.path.join(args.output_dir, "assemblies"), "quast": os.path.join(args.output_dir, "quast_output"), "prokka": os.path.join(args.output_dir, "prokka_output"), "maker": os.path.join(args.output_dir, "maker_output") }
    for d in dirs.values(): os.makedirs(d, exist_ok=True)

    logger.info(colored("[Step 1/6] Downloading reference data...", "green"))
    ref_genome = download_reference_genome(config['reference_genome_url'], dirs["ref"], logger)
    gff_genome = download_reference_genome(config['gff_url'], dirs["ref"], logger)

    logger.info(colored("[Step 2/6] Preparing list of assemblies...", "green"))
    if args.assembly_file:
        assembly_ids = read_assembly_file(args.assembly_file, logger)
        filtered_df = pd.DataFrame({"AssemblyAccession": assembly_ids})
    else:
        logger.info("No assembly file provided. Searching for assemblies from India by default.")
        df = fetch_assembly_data(config['organism_query'], logger)
        # The location is now hardcoded to "India"
        filtered_df = filter_rows_by_cities(df, "India", logger)
        if filtered_df.empty:
            logger.error("No assemblies from India matched the query. To process other assemblies, please provide a list via the --assembly-file flag.")
            return
            
    logger.info(colored("[Step 3/6] Downloading assemblies...", "green"))
    accessions = download_assemblies(filtered_df, dirs["asm"], config['group'], args.num_assemblies, logger)
    if not accessions:
        logger.error("No assemblies downloaded.")
        return

    logger.info(colored("[Step 4/6] Finalizing assembly files...", "green"))
    decompress_and_rename_assemblies(dirs["asm"], accessions, logger)

    logger.info(colored("[Step 5/6] Running QUAST...", "green"))
    run_quast(accessions, dirs["asm"], dirs["quast"], ref_genome, gff_genome, config['quast_threads'], config.get('quast_params', {}), logger)

    if config['organism_type'] == 'PK':
        logger.info(colored("[Step 6/6] Running Prokaryotic Annotation (Prokka)...", "green"))
        run_prokka(accessions, dirs["asm"], dirs["prokka"], config.get('prokka_kingdom', 'Bacteria'), config.get('prokka_params', {}), logger)
    elif config['organism_type'] == 'EK':
        logger.info(colored("[Step 6/6] Running Eukaryotic Annotation (MAKER)...", "green"))
        for accession in accessions:
            run_maker_for_accession(accession, config, dirs["asm"], dirs["maker"], logger)
    
    logger.info(colored("\nPIPELINE COMPLETED SUCCESSFULLY!", "cyan"))

if __name__ == "__main__":
    main()
