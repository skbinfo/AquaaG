#!/usr/bin/env python3
"""
Download reference genome + India-specific assemblies for species in a kingdom from NCBI.
Optimized to process species one at a time (faster for small --num values).
"""

import os
import re
import sys
import gzip
import shutil
import argparse
import subprocess
import logging
from concurrent.futures import ThreadPoolExecutor, as_completed
from Bio import Entrez
import pandas as pd

# ---------- CONFIG ----------
Entrez.email = "your@email.com"  # Change this
OUTPUT_DIR = "ncbi_kingdom_data"
THREADS = 5  # Parallel assembly downloads
MAX_SPECIES_FETCH = 200  # max records per esearch batch (small to avoid overload)

# ---------- INDIA FILTER ----------
INDIAN_CITIES = [
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
INDIA_PATTERN = re.compile(r"\b(" + "|".join(map(re.escape, INDIAN_CITIES)) + r")\b", re.IGNORECASE)

# ---------- LOGGING ----------
def setup_logger(output_dir):
    log_file = os.path.join(output_dir, "india_download.log")
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        handlers=[
            logging.FileHandler(log_file, mode='w'),
            logging.StreamHandler(sys.stdout)
        ]
    )
    return logging.getLogger(__name__)

# ---------- UTILS ----------
def download_and_extract(url, out_path, logger):
    try:
        # Check if URL exists
        check = subprocess.run(f"wget --spider -q {url}", shell=True)
        if check.returncode != 0:
            logger.warning(f"URL does not exist: {url}")
            return False

        subprocess.run(f"wget -q -O temp.gz {url}", shell=True, check=True)
        with gzip.open("temp.gz", 'rb') as f_in, open(out_path, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
        os.remove("temp.gz")
        logger.info(f"Downloaded: {url}")
        return True
    except Exception as e:
        logger.warning(f"Failed to download {url}: {e}")
        if os.path.exists("temp.gz"):
            os.remove("temp.gz")
        return False


def fetch_reference_for_species(species_name, output_dir, logger):
    """Fetch the latest valid reference/representative genome for a species."""
    logger.info(f"Fetching reference genome for {species_name}")

    # Broader query: reference OR representative OR latest
    search_terms = [
        f'"{species_name}"[Organism] AND "reference genome"[Filter] AND latest[filter]',
        f'"{species_name}"[Organism] AND "representative genome"[Filter] AND latest[filter]',
        f'"{species_name}"[Organism] AND latest[filter]'
    ]

    for term in search_terms:
        handle = Entrez.esearch(db="assembly", term=term, retmax=5, sort="relevance")
        record = Entrez.read(handle)
        handle.close()

        if not record["IdList"]:
            continue

        summary_handle = Entrez.esummary(db="assembly", id=",".join(record["IdList"]), report="full")
        summaries = Entrez.read(summary_handle)
        summary_handle.close()

        for doc in summaries["DocumentSummarySet"]["DocumentSummary"]:
            ftp_path = doc["FtpPath_RefSeq"] or doc["FtpPath_GenBank"]
            accession = doc["AssemblyAccession"]

            if not ftp_path:
                continue

            ftp_path = ftp_path.replace("ftp://", "https://")
            base_name = os.path.basename(ftp_path)

            fasta_url = f"{ftp_path}/{base_name}_genomic.fna.gz"
            gff_url_gz = f"{ftp_path}/{base_name}_genomic.gff.gz"
            gff_url_plain = f"{ftp_path}/{base_name}_genomic.gff"

            fasta_ok = download_and_extract(fasta_url, os.path.join(output_dir, f"{accession}.fna"), logger)
            if not fasta_ok:
                continue

            gff_ok = download_and_extract(gff_url_gz, os.path.join(output_dir, f"{accession}.gff"), logger) \
                     or download_and_extract(gff_url_plain, os.path.join(output_dir, f"{accession}.gff"), logger)

            if not gff_ok:
                continue

            logger.info(f"Reference genome downloaded for {species_name} ({accession})")
            return True

    logger.warning(f"No valid downloadable reference genome found for {species_name}")
    return False

def fetch_india_specific_assemblies(species_name, out_dir, logger, meta_list):
    """Fetch India-specific assemblies for a species."""
    logger.info(f"Fetching India-specific assemblies for {species_name}")

    cmd = (
        f'esearch -db assembly -query "{species_name}[Organism] AND latest[filter]" | '
        f'esummary | '
        f'xtract -pattern DocumentSummary '
        f'-element AssemblyAccession,SubmitterOrganization,FtpPath_RefSeq,FtpPath_GenBank'
    )
    result = subprocess.run(cmd, shell=True, text=True, capture_output=True)
    if result.returncode != 0 or not result.stdout.strip():
        return

    rows = []
    for line in result.stdout.strip().split("\n"):
        parts = line.split("\t")
        while len(parts) < 4:
            parts.append("-")
        rows.append(parts)

    df = pd.DataFrame(rows, columns=["Accession", "Submitter", "RefSeq", "GenBank"])
    df = df[df["Submitter"].str.contains(INDIA_PATTERN, na=False)]
    if df.empty:
        return

    for _, row in df.iterrows():
        ftp_path = row["RefSeq"] if row["RefSeq"] != "-" else row["GenBank"]
        if ftp_path and ftp_path != "-":
            ftp_path = ftp_path.replace("ftp://", "https://")
            base_name = os.path.basename(ftp_path)
            fasta_url = f"{ftp_path}/{base_name}_genomic.fna.gz"
            out_file = os.path.join(out_dir, f"{row['Accession']}.fna")

            if download_and_extract(fasta_url, out_file, logger):
                meta_list.append({
                    "Species": species_name,
                    "Accession": row["Accession"],
                    "Submitter": row["Submitter"],
                    "FTP_Path": ftp_path,
                    "Saved_File": out_file
                })

def get_species_names_from_kingdom(kingdom, start, batch_size):
    """Get a batch of species names from a kingdom."""
    handle = Entrez.esearch(
        db="assembly",
        term=f'"{kingdom}"[Organism]',
        retstart=start,
        retmax=batch_size
    )
    record = Entrez.read(handle)
    handle.close()
    if not record["IdList"]:
        return []
    summaries_handle = Entrez.esummary(db="assembly", id=",".join(record["IdList"]), report="full")
    summaries = Entrez.read(summaries_handle)
    summaries_handle.close()
    return list(dict.fromkeys([doc["SpeciesName"] for doc in summaries["DocumentSummarySet"]["DocumentSummary"]]))


def main(args):
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    logger = setup_logger(OUTPUT_DIR)

    count = 0
    start = 0
    while args.num is None or count < args.num:
        species_batch = get_species_names_from_kingdom(args.kingdom, start, MAX_SPECIES_FETCH)
        if not species_batch:
            break

        for species in species_batch:
            if count >= args.num:
                break

            logger.info(f"Checking India-specific assemblies for: {species}")
            meta_list = []
            temp_dir = os.path.join(OUTPUT_DIR, "temp_check")
            os.makedirs(temp_dir, exist_ok=True)

            fetch_india_specific_assemblies(species, temp_dir, logger, meta_list)
            shutil.rmtree(temp_dir, ignore_errors=True)

            if not meta_list:
                logger.info(f"No India-specific assemblies for {species}, skipping...")
                continue

            logger.info(f"India-specific assemblies found for {species}, downloading...")
            species_dir = os.path.join(OUTPUT_DIR, species.replace(" ", "_"))
            os.makedirs(species_dir, exist_ok=True)

            ref_dir = os.path.join(species_dir, "reference")
            os.makedirs(ref_dir, exist_ok=True)
            fetch_reference_for_species(species, ref_dir, logger)

            asm_dir = os.path.join(species_dir, "assemblies")
            os.makedirs(asm_dir, exist_ok=True)
            meta_list = []
            fetch_india_specific_assemblies(species, asm_dir, logger, meta_list)

            if meta_list:
                pd.DataFrame(meta_list).to_csv(
                    os.path.join(species_dir, "india_assemblies_metadata.tsv"),
                    sep="\t", index=False
                )

            count += 1

        start += len(species_batch)  # move forward only by actual batch size


    logger.info("Completed.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Download NCBI reference + India-specific assemblies for a kingdom"
    )

    parser.add_argument(
        "-s", "--kingdom", required=True,
        help="Kingdom name (e.g., fungi, bacteria, viruses)"
    )
    parser.add_argument(
    "--num", type=int, default=None,
    help="Number of species to process (if not given, process all)"
    ) 

    args = parser.parse_args()

    main(args)

