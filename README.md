# AquaaG: Automated Quality Assessment and Annotation of Genomes

AquaaG is a Python-based pipeline for:

* Downloading genomic assemblies from NCBI
* Performing quality assessment with **QUAST**
* Annotating:

  * **Prokaryotic genomes** with **Prokka**
  * **Eukaryotic genomes** with **BRAKER3** (inside Docker, with repeat masking)
* Evaluating annotation completeness with **BUSCO**

It supports **India-specific filtering** of assemblies (based on submitter metadata) and can also process **user-provided assembly IDs** via a file (e.g. `assembly.txt`). The pipeline is designed for production-grade bioinformatics workflows and supports both **prokaryotic (PK)** and **eukaryotic (EK)** organisms.

---

## ✨ Features

* **Assembly Fetching**

  * Downloads assemblies from the NCBI Assembly database
  * Optional filtering for **India-specific submissions** (by submitter organization)

* **Quality Assessment**

  * Uses **QUAST** to evaluate assembly quality against a reference genome + GFF

* **Annotation**

  * **Prokaryotic (PK)**: Annotation with **Prokka**
  * **Eukaryotic (EK)**:

    * Repeat library construction using **RepeatModeler**
    * Genome masking with **RepeatMasker**
    * **BRAKER3** gene prediction inside a **Docker** container (GeneMark-ETP, AUGUSTUS)

* **Completeness Evaluation**

  * Runs **BUSCO** on predicted proteins (Prokka or BRAKER3 output)

* **Flexible Modes**

  * **Species mode**: `--species "Genus species"`
  * **Kingdom mode**: `--kingdom Bacteria/Fungi/...` with automatic species discovery
  * **Assembly list mode**: `--assembly-file assembly.txt` for explicit accessions
  * **India-only filtering**: `-I / --india-only`
  * **Any source (default)**: `-a / --any-source`

* **Parallel & Configurable**

  * Multi-threaded components (QUAST, Prokka, BRAKER3, BUSCO)
  * All major parameters configured via **YAML** files

---

## 🧩 Requirements

* **Operating System**
  Linux (tested on Ubuntu-like environments; HPC clusters should also work)

* **Conda**
  Miniconda or Anaconda installed and available in `PATH`

* **Docker**
  AquaaG uses **BRAKER3 via Docker** for eukaryotic annotation. **Docker Engine must be installed and running before you use AquaaG.**

  * Installation instructions are available from the official Docker documentation:
    [https://docs.docker.com/engine/install/](https://docs.docker.com/engine/install/)
  * After installation, verify Docker with:

    ```bash
    docker run hello-world
    ```

* **Internet Access**

  * Required for:

    * NCBI downloads (assemblies, references, proteins)
    * BUSCO lineage datasets (first-time download)
    * Docker image pulls (BRAKER3 container)

* **Disk Space**

  * Recommended: **10–20 GB** free space (more for many or large genomes)

---

## 📦 Installation and Environment Setup

### 1. Clone the Repository

```bash
git clone https://github.com/skbinfo/AquaaG.git
cd AquaaG
```

### 2. Ensure Docker Is Installed

Install Docker Engine following the official instructions:
[https://docs.docker.com/engine/install/](https://docs.docker.com/engine/install/)

Then test that Docker works:

```bash
docker run hello-world
```

If this succeeds, Docker is correctly installed and the current user can run containers.

### 3. Run the Setup Script

AquaaG provides a setup helper script that:

* Creates two Conda environments:

  * `assembly_tool_env` (main pipeline tools)
  * `braker_env` (Docker client / helper tools)
* Pulls the **BRAKER3** Docker image
* Copies **AUGUSTUS config** from inside the Docker image to your host
* Configures `AUGUSTUS_CONFIG_PATH`
* Registers your **GeneMark-ETP `gm_key`** for BRAKER3

Run:

```bash
bash setup_new.sh
```

During setup, you will be prompted for:

* A target directory for `AUGUSTUS_CONFIG_PATH`
* The local path to your **GeneMark-ETP `gm_key`** file (downloaded from the GeneMark website)

After the script completes, run:

```bash
source ~/.bashrc
```

or open a **new shell** so that `AUGUSTUS_CONFIG_PATH` and `BRAKER_GM_KEY` are available.

### 4. Activate the Main Environment

Most of the pipeline is run from the **main environment**:

```bash
conda activate assembly_tool_env
```

### 5. Verify Installation

Check that core tools are visible:

```bash
quast.py --version
prokka --version
busco --version
docker --version
```

BRAKER3 itself is run **inside Docker**, so you do not need a host-level `braker.pl` installation.

---

## ⚙️ Configuration

AquaaG uses separate YAML configuration files for **eukaryotic** and **prokaryotic** pipelines:

* **Eukaryotic config:** `Eu_config.yaml`
* **Prokaryotic config:** `Pr_config.yaml`

Edit these files **before** running the pipeline.

---

### 🧬 Eukaryotic Config (`Eu_config.yaml`)

Example:

```yaml
#### Eukaryotic Annotation Config (BRAKER3) ####

email: "your@email.com"

organism: "Arabidopsis thaliana"
organism_type: "EK"

quast_threads: 50

busco_lineage: "fungi_odb10"
busco_params:
  cpu: 50

braker_params:
  cpus: 50
  # Optional extras:
  # fungus: true
  # UTR: "on"
  # protein_evidence: "/absolute/path/to/custom_proteins.fasta"

quast_params: {}
```

**Key fields:**

* `email` – used for NCBI Entrez queries (required).
* `organism` – scientific name used to fetch **reference genome + GFF + proteins** from NCBI.
* `organism_type` – must be `"EK"` for the eukaryotic pipeline.
* `quast_threads` – number of CPU threads for QUAST.
* `busco_lineage` – BUSCO lineage database; e.g.:

  * `embryophyta_odb10` for plants
  * `fungi_odb10` for fungi
  * `metazoa_odb10` for animals
* `busco_params` – extra BUSCO parameters; `cpu` is recommended.
* `braker_params` – configuration for BRAKER3:

  * `cpus` (required) – number of threads for BRAKER
  * optional flags such as `fungus: true`, `UTR: "on"`, or a custom `protein_evidence` path.
* `quast_params` – additional QUAST options (may be left empty `{}`).

---

### 🦠 Prokaryotic Config (`Pr_config.yaml`)

Example:

```yaml
#### Prokaryotic Annotation Config ####

email: "your@email.com"

organism: "Mycobacterium tuberculosis"
group: "bacteria"

organism_type: "PK"

quast_threads: 8

busco_lineage: "bacteria_odb10"
busco_params:
  cpu: 8

prokka_kingdom: "Bacteria"
prokka_params:
  cpus: 8

quast_params: {}
```

**Key fields:**

* `organism` – reference organism used to fetch **reference genome + GFF + proteins**.
* `group` – NCBI download group for `ncbi-genome-download` (e.g. `bacteria`, `fungi`, `viral`).
* `organism_type` – must be `"PK"` for the prokaryotic pipeline.
* `quast_threads` – number of CPU threads for QUAST.
* `busco_lineage` – BUSCO lineage for prokaryotes (`bacteria_odb10`, `archaea_odb10`, etc.).
* `busco_params` – extra BUSCO parameters.
* `prokka_kingdom` – Prokka kingdom (e.g. `"Bacteria"`, `"Archaea"`).
* `prokka_params` – additional Prokka options (e.g. `cpus`, `genus`, `species`, `strain`).

---

### 🔖 Assembly List File (`assembly.txt`)

For **assembly-file mode**, create a text file with one NCBI Assembly Accession per line:

```text
GCF_000000000.1
GCA_123456789.1
```

Use this file with `--assembly-file assembly.txt`.

---

## 🚀 Running the Pipeline

The main entry point is **`AquaaG.py`**.

General pattern:

```bash
python AquaaG.py -c <CONFIG.yaml> -o <OUTPUT_DIR> [MODE OPTIONS] [FILTER OPTIONS]
```

---

### 1. Species Mode (`--species`)

Process one specific species by scientific name.

Example (eukaryote):

```bash
python AquaaG.py \
  -c Eu_config.yaml \
  -o output_eu \
  --species "Arabidopsis thaliana" \
  -I \
  --num-assemblies 1
```

* `--species` – species name for assembly search.
* `-I / --india-only` – only process assemblies whose submitter organization matches India-related keywords.
* `--num-assemblies` – max number of assemblies to download for this species.

If you **omit `-I`**, the pipeline processes assemblies from **any submitter**.

---

### 2. Kingdom Mode (`--kingdom`)

Automatically discover species within a kingdom and process a subset of them.

Example (bacteria, any submitter):

```bash
python AquaaG.py \
  -c Pr_config.yaml \
  -o output_pk_kingdom \
  --kingdom Bacteria \
  --num-species 3 \
  --num-assemblies 1 \
  -a
```

* `--kingdom` – kingdom-level search term (e.g. `Bacteria`, `Fungi`).
* `--num-species` – maximum number of species to process.
* `--num-assemblies` – max assemblies per species.
* `-a / --any-source` – process assemblies from any submitter (no India-only filter).

If you use `-I` instead of `-a`, only species with at least one **India-submitted** assembly are selected.

---

### 3. Assembly List Mode (`--assembly-file`)

Process **specific accessions** defined in a text file, using `organism` in the config as the reference.

Example (prokaryote):

```bash
python AquaaG.py \
  -c Pr_config.yaml \
  -o output_pk_assemblies \
  --assembly-file assembly.txt \
  -a
```

* Uses `organism` from the config as the reference species.
* Downloads only the listed accessions.
* India filtering is not applied in this mode (you control the accessions explicitly).

---

## 📁 Output Structure

Depending on mode, AquaaG creates an output directory structure like:

### Species / Kingdom Modes

```text
<OUTPUT_DIR>/
  <Species_Slug>/                # e.g. Arabidopsis_thaliana
    reference_data/
      *.fna                      # reference genome
      *.gff                      # reference annotation
      *.faa                      # reference proteins (if available)
    assemblies/
      GCF_XXXXXXX.fna
      ...
    quast_output/
      report.txt
      report.html
      ...
    prokka_output/               # PK mode
      GCF_XXXXXXX/
        *.gff
        *.faa
        ...
    braker_output/               # EK mode
      GCF_XXXXXXX/
        braker.aa                # predicted proteins
        braker.gtf               # gene models
        GeneMark-ETP/            # GeneMark internals
        ...
    busco_output/
      GCF_XXXXXXX/
        short_summary.specific.<lineage>.*.txt
    pipeline_main.log
```

### Assembly-File Mode

In `--assembly-file` mode (single reference organism), the structure is:

```text
<OUTPUT_DIR>/
  reference_data/
  assemblies/
  quast_output/
  prokka_output/      # PK mode
  braker_output/      # EK mode
  busco_output/
  pipeline_main.log
```

---

## 🇮🇳 Customizing India-Specific Filtering

India-specific filtering is implemented by the `filter_for_india` function in `AquaaG.py`. It uses a list of keywords (`indian_cities`) matched against the `SubmitterOrganization` field from NCBI.

Snippet:

```python
indian_cities = [
    "IND", "Indian", "india", "India",
    "Agartala", "Ahmedabad", "Aizawl", "Ajmer",
    "Allahabad", "Amritsar", "Anand", "Avikanagar",
    "Aurangabad", "Amravati", "Bangalore", "Bareilly",
    # ... many more cities/institutes ...
    "ICMR", "IMTECH", "CSIR", "IIT", "NIT", "IISC",
    "AIIMS", "SRM", "CDRI", "ICGEB"
]
```

### Extend the List

1. Open `AquaaG.py` and locate `filter_for_india`.
2. Add new city/institution keywords to `indian_cities`, for example:

```python
indian_cities = [
    ...
    "Tamil Nadu",
    "NCCS"
]
```

3. Save `AquaaG.py` and rerun AquaaG with `-I` / `--india-only`.

**Best practices:**

* Avoid overly generic terms like `"Institute"` to reduce false positives.
* Inspect submitter organizations manually via:

  ```bash
  esearch -db assembly -query "Mycobacterium tuberculosis[Organism]" \
    | esummary \
    | xtract -pattern DocumentSummary -element AssemblyAccession,SubmitterOrganization
  ```

---

## 🛠 Troubleshooting

### 1. No Assemblies Found

* Check `pipeline_main.log` under your output directory.
* Verify that `organism` or `kingdom` in the config file is correct.
* In India-only mode (`-I`), there may simply be no India-submitted assemblies.
* Use `--any-source` or `--assembly-file` as a fallback.

### 2. Docker / BRAKER3 Issues

If BRAKER fails or logs mention missing `prothint.gff`, `evidence.gff`, or GeneMark license issues:

* Confirm Docker is running:

  ```bash
  docker run hello-world
  ```

* Check that the BRAKER image is available:

  ```bash
  docker images | grep braker
  ```

* Verify environment variables:

  * `AUGUSTUS_CONFIG_PATH` points to the directory created by `setup_new.sh`.
  * `BRAKER_GM_KEY` points to your `gm_key` file under `.braker_keys/`.

AquaaG attempts an **automatic recovery** if GeneMark-ETP complains about missing ProtHint outputs by copying found `prothint.gff` / `evidence.gff` into expected directories and re-running BRAKER.

### 3. BUSCO Failures

* Check that `busco_lineage` in the config matches your organism type (`bacteria_odb10`, `archaea_odb10`, `fungi_odb10`, `embryophyta_odb10`, etc.).
* Ensure BUSCO lineages are downloaded (BUSCO will log if databases are missing).
* Confirm the BUSCO input file exists:

  * PK: Prokka `.faa` file under `prokka_output/`
  * EK: `braker.aa` under `braker_output/`

### 4. RepeatModeler / RepeatMasker Issues (EK Only)

* Confirm tools are installed in `assembly_tool_env` and in `PATH`:

  ```bash
  BuildDatabase -h
  RepeatModeler -h
  RepeatMasker
  ```

* Check the log for missing dependencies or incompatible options.

### 5. Conda / Environment Issues

If tools are missing or versions look wrong:

```bash
conda env remove -n assembly_tool_env
conda env remove -n braker_env

conda env create -f environment_main.yml
conda env create -f environment_braker.yml

bash setup_new.sh
```

---

## 🤝 Contributions & Acknowledgments

* **Contributions** are welcome via pull requests or issues:

  * Bug reports
  * Feature requests
  * Documentation improvements

* **Tools and Libraries**

  * **BRAKER3**, **GeneMark-ETP**, **AUGUSTUS**
  * **Prokka**
  * **QUAST**
  * **BUSCO**
  * **RepeatModeler**
  * **RepeatMasker**
  * **pandas**, **Biopython**, and others

* **Community**

  * Inspired by open-source efforts in genome assembly and annotation
  * Thanks to all users and collaborators who tested and suggested improvements

---

# 🇮🇳 `download.py`

**Download Reference Genomes + India-Specific Assemblies from NCBI**

`download.py` is a companion tool in this repository that focuses on **downloading and cataloging assemblies** from NCBI, with special support for **India-specific submissions**.

---

## 🔍 Features

* Detects **India-specific assemblies** using submitter info (cities, institutes, keywords)
* Downloads **reference or representative genomes** (FASTA + GFF)
* Fetches **India-specific assemblies** for multiple species in a kingdom
* Produces a **metadata table** (e.g. `india_assemblies_metadata.tsv`)
* Parallelized downloads with logging and error handling

---

## ▶️ Usage

```bash
# Process all species in a kingdom (e.g. fungi)
python download.py -s fungi

# Process only a fixed number of species
python download.py -s fungi --num 5

# Example: Cyanobacteria, first 3 species
python download.py -s Cyanobacteriota --num 3
```

You can use outputs from `download.py` to:

* Select species to feed into `AquaaG.py` (species or assembly-file modes)
* Curate `assembly.txt` lists for `--assembly-file` workflows

---

## 🙏 Acknowledgments for `download.py`

* Uses the same India-filtering strategy as AquaaG
* Built on top of:

  * NCBI **Entrez**, **esearch/esummary/xtract**
  * **pandas** for metadata management

For questions or issues related to `AquaaG` or `download.py`, please **open an issue** on the GitHub repository.
