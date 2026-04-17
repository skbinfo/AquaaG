# AquaaG: Automated Quality Assessment and Annotation of Genomes

> A comprehensive, fully containerized Python-based pipeline for downloading, quality-assessing, annotating, and functionally mapping genomic assemblies from NCBI.

---

## 📋 Table of Contents

- [Overview](#overview)
- [Features](#features)
- [Requirements](#requirements)
- [Installation](#installation)
- [Running the Pipeline](#running-the-pipeline)
- [Output Structure](#output-structure)
- [Regional Filtering](#regional-filtering)
- [Troubleshooting](#troubleshooting)
- [download.py Companion Tool](#downloadpy-companion-tool)
- [Contributions & Acknowledgments](#contributions--acknowledgments)

---

## Overview

AquaaG automates the full genome analysis workflow:

- **Download** genomic assemblies from NCBI
- **Quality assessment** with QUAST
- **Structural annotation** — Prokaryotes via Prokka, Eukaryotes via BRAKER3 (with automated repeat masking)
- **Completeness evaluation** with BUSCO
- **Functional annotation mapping** with EggNOG-mapper
- **Interactive reporting** via a Chart.js Dashboard and MultiQC report

Supports **regional metadata filtering** (e.g., country or institution-specific submissions), user-provided assembly ID lists, and direct local `.fna` files. Designed for production-grade bioinformatics with **zero host-system dependencies**.

---

## ✨ Features

| Feature | Description |
|---|---|
| 🐳 **100% Dockerized** | Eliminates dependency conflicts. Runs entirely via a master Docker container, including heavy tools like BRAKER3 and EggNOG-mapper. |
| ⏯️ **Zero-API Resumption** | Caches NCBI metadata locally so interrupted runs resume instantly without re-querying NCBI. |
| 🌍 **Regional Filtering** | Optional `-R` / `--region-only` flag to process submissions from specific geographic locations or institutions. |
| 📊 **Quality Assessment** | QUAST evaluation against a reference genome + GFF. |
| 🧬 **End-to-End Annotation** | PK: Prokka. EK: RepeatModeler → RepeatMasker → BRAKER3. |
| 🔬 **Functional Mapping** | COG, KEGG, and GO term mapping via EggNOG-mapper. |
| 📈 **Automated Reporting** | Unified `AquaG_Dashboard.html` (Chart.js) and `multiqc_report.html`. |
| 🗂️ **Flexible Input Modes** | Species name, Kingdom, accession list file, or local FASTA directory. |

---

## 🧩 Requirements

Because AquaaG is fully containerized, **no Conda, Python, or bioinformatics tools** need to be installed on your host machine.

- **OS:** Linux or macOS (Ubuntu-like environments and HPC clusters with Docker access supported)
- **Docker Engine:** Must be installed and running
  - Install: https://docs.docker.com/engine/install/
  - Verify: `docker run hello-world`
- **Internet Access:** Required for NCBI downloads, BUSCO lineages, and Docker image pulls
- **Disk Space:** 10–20 GB recommended (~50 GB additional if downloading the EggNOG functional database)

---

## 📦 Installation

### 1. Clone the Repository

```bash
git clone https://github.com/skbinfo/AquaaG.git
cd AquaaG
```

### 2. Place Your GeneMark Key *(Required for Eukaryotes only)*

To run BRAKER3 on eukaryotic genomes, obtain a free academic GeneMark-ETP license key:

1. Download `gm_key` from: https://exon.gatech.edu/GeneMark/license_download.cgi
2. Place the `gm_key` file directly in your `AquaaG/` directory. The pipeline will auto-detect and mount it.

### 3. Run the docker build

```bash
docker build -t aquaag-pipeline:latest .
```

This will:

- Build the `aquaag-pipeline:latest` Docker image
- Detect your `gm_key`
- Optionally pre-download the 50 GB EggNOG database into `./eggnog_db/` (required for `--run-func`)

---

## 🚀 Running the Pipeline

AquaaG is executed entirely through `docker run`, using environment variables (`-e`) to configure the pipeline on the fly.

### Standard Docker Boilerplate

Every command must begin with this block:

```bash
docker run --rm -it \
  -u $(id -u):$(id -g) \
  --group-add $(stat -c '%g' /var/run/docker.sock) \
  -v /var/run/docker.sock:/var/run/docker.sock \
  -v "$PWD":"$PWD" -w "$PWD" \
  -e HOME=/tmp
```

---

### Mode 1: Species Mode (`-s`)

Process one specific species by scientific name.

**Example — Eukaryote, 1 assembly, functional annotation enabled:**

```bash
docker run --rm -it \
  -u $(id -u):$(id -g) --group-add $(stat -c '%g' /var/run/docker.sock) \
  -v /var/run/docker.sock:/var/run/docker.sock \
  -v "$PWD":"$PWD" -w "$PWD" -e HOME=/tmp \
  -e EMAIL="your@email.com" \
  -e TYPE="EK" \
  -e GROUP="plant" \
  -e THREADS=40 \
  -e BUSCO="embryophyta_odb10" \
  aquaag-pipeline -o arabidopsis_results -s "Arabidopsis thaliana" --num-assemblies 1 --run-func
```

> **Note:** For eukaryotes, BRAKER3 will interactively ask for RNA-Seq or protein evidence paths on the first assembly, then autonomously cache those choices for all subsequent assemblies.

---

### Mode 2: Kingdom Mode (`-k`)

Automatically discover species within a kingdom and process a subset.

**Example — Prokaryote, 3 bacteria species:**

```bash
docker run --rm -it \
  -u $(id -u):$(id -g) --group-add $(stat -c '%g' /var/run/docker.sock) \
  -v /var/run/docker.sock:/var/run/docker.sock \
  -v "$PWD":"$PWD" -w "$PWD" -e HOME=/tmp \
  -e EMAIL="your@email.com" \
  -e TYPE="PK" \
  -e GROUP="bacteria" \
  -e THREADS=20 \
  -e BUSCO="bacteria_odb10" \
  aquaag-pipeline -o bulk_bacteria_results -k "Bacteria" --num-species 3 --num-assemblies 1 -a --run-func
```

---

### Mode 3: Local Assembly Directory (`--assembly-dir`)

Annotate your own unpublished genomes without downloading from NCBI.

**Example — Prokaryote, custom FASTA files in `./my_genomes/`:**

```bash
docker run --rm -it \
  -u $(id -u):$(id -g) --group-add $(stat -c '%g' /var/run/docker.sock) \
  -v /var/run/docker.sock:/var/run/docker.sock \
  -v "$PWD":"$PWD" -w "$PWD" -e HOME=/tmp \
  -e EMAIL="your@email.com" \
  -e TYPE="PK" \
  -e ORGANISM="Vibrio cholerae" \
  -e THREADS=20 \
  -e BUSCO="vibrionales_odb10" \
  aquaag-pipeline -o custom_vibrio_results --assembly-dir ./my_genomes/ --run-func
```

> **Note:** `-e ORGANISM` is still required so QUAST knows which reference genome to download for quality grading.

---

## 📁 Output Structure

```text
<OUTPUT_DIR>/
├── AquaG_Dashboard.html          # Interactive Chart.js dashboard (QUAST/BUSCO)
├── MultiQC_Report/               # Aggregated HTML report
├── pipeline_main.log             # Master execution log
└── <Species_Slug>/               # e.g., Vibrio_cholerae/
    ├── reference_data/           # Cached NCBI reference genomes (.fna, .gff, .faa)
    ├── assemblies/               # Decompressed target FASTA files
    ├── quast_output/             # Contiguity metrics
    ├── prokka_output/            # (PK mode) Structural annotations (.gff, .faa)
    ├── braker_output/            # (EK mode) Structural annotations (braker.aa, .gtf)
    ├── busco_output/             # Completeness grading
    ├── functional_annotation/    # EggNOG-mapper outputs (.emapper.annotations)
    └── <Species>_metadata.tsv    # Cached NCBI metadata
```
---

## 🌍 Regional Filtering

AquaaG supports metadata filtering (`-R` / `--region-only`) to exclusively process submissions from specific regions or institutions.

By default, the pipeline includes a built-in list of Indian cities and institutes (e.g., `"NIPGR"`, `"Delhi"`, `"CSIR"`).

**To customize for your own region:**

1. Open `AquaaG.py`
2. Locate the `filter_by_region` (or `filter_for_india`) function
3. Update the `indian_cities` array with your desired keywords (e.g., `["USA", "NIH", "Oxford"]`)
4. Rebuild the Docker image:
```bash
   docker build -t aquaag-pipeline:latest .
```
5. Run the pipeline with the `-R` flag

---
## Quick Start: Verifying Your Installation

To verify that AquaaG and Docker are correctly installed and functioning, we provide a minimal toy dataset. This test will download a very small bacterial genome, assess its quality, annotate it, and evaluate its completeness in just a few minutes.

**1. Run the Prokaryotic Test:**
```bash
docker run --rm -it \
  -u $(id -u):$(id -g) --group-add $(stat -c '%g' /var/run/docker.sock) -v /var/run/docker.sock:/var/run/docker.sock -v "$PWD":"$PWD" -w "$PWD" -e HOME=/tmp \
  -e EMAIL="test@example.com" \
  -e TYPE="PK" \
  -e THREADS=4 \
  -e BUSCO="bacteria_odb10" \
  aquaag-pipeline -o test_output --assembly-file test_data/test_pk_assembly.txt -a
---
## 🛠️ Troubleshooting

### ❌ `FATAL: 'gm_key' not found!`

You are running Eukaryotic mode (`TYPE="EK"`) but the `gm_key` file is missing from the directory where you execute the `docker run` command. Move the key to your current directory and rerun.

### ⚠️ EggNOG Database Warnings
If you pass `--run-func` but see a warning that functional annotation was skipped, it means the required 50GB database is missing. As intended, the EggNOG database download is processed securely via Docker to avoid local environment conflicts. 

You can resolve this by running `bash setup.sh` and accepting the prompt to download `eggnog_db`. Alternatively, you can trigger the Docker download manually in your current directory:
```bash
mkdir -p "$PWD/eggnog_db"
docker run --rm -it \
  -v "$PWD/eggnog_db":/eggnog_db \
  --entrypoint /bin/bash \
  quay.io/biocontainers/eggnog-mapper:2.1.12--pyhdfd78af_2 \
  -c "sed -i 's/eggnogdb.embl.de/eggnog5.embl.de/g' \$(which download_eggnog_data.py) && download_eggnog_data.py -P -y --data_dir /eggnog_db && chown -R \$(id -u):\$(id -g) /eggnog_db"
```
### 🔄 Pipeline Crashes or NCBI Timeouts

Do not panic. AquaaG has **Zero-API Resumption**. Simply re-run the exact same `docker run` command. The pipeline reads its local cache, skips all completed steps, and resumes from where it stopped.

### 💬 BRAKER3 Interactive Prompts

AquaaG caches your BRAKER3 answers (mode, SRA IDs, protein paths) after the first assembly of a species and autonomously applies them to all subsequent assemblies, enabling hands-free bulk processing.

---

## 🌐 `download.py` Companion Tool

`download.py` focuses exclusively on downloading and cataloging assemblies from NCBI, with special support for regional submissions.

### Features

- Detects regional assemblies using submitter info (cities, institutes, keywords)
- Downloads reference or representative genomes (FASTA + GFF)
- Produces a metadata table (e.g., `assemblies_metadata.tsv`)
- Parallelized downloads with logging and error handling

### Usage

```bash
# Process all species in a kingdom (e.g., fungi)
python download.py -s fungi

# Process only a fixed number of species
python download.py -s fungi --num 5
```

Outputs from `download.py` can be used to curate custom `assembly.txt` lists for AquaaG's `--assembly-file` workflow.

---

## 🤝 Contributions & Acknowledgments

Contributions are welcome via pull requests or issues for bug reports, feature requests, or documentation improvements.

**Built on open-source tools including:**
[BRAKER3](https://github.com/Gaius-Augustus/BRAKER) · [GeneMark-ETP](http://exon.gatech.edu/GeneMark/) · [AUGUSTUS](https://github.com/Gaius-Augustus/Augustus) · [Prokka](https://github.com/tseemann/prokka) · [QUAST](https://github.com/ablab/quast) · [BUSCO](https://busco.ezlab.org/) · [EggNOG-mapper](https://github.com/eggnogdb/eggnog-mapper) · [MultiQC](https://multiqc.info/) · [RepeatModeler](https://github.com/Dfam-consortium/RepeatModeler)
