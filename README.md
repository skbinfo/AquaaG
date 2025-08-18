# AquaaG: Automated Quality Assessment and Annotation of Genomes


AquaaG is a Python-based pipeline for downloading genomic assemblies from NCBI, performing quality assessment with QUAST, annotating prokaryotic genomes with Prokka or eukaryotic genomes with MAKER, and evaluating completeness with BUSCO. It supports filtering for India-specific assemblies by default and allows processing user-provided assembly IDs via a file (e.g., genome.txt). The pipeline is designed for bioinformatics workflows and handles both prokaryotic (PK) and eukaryotic (EK) organisms.

---

## Features

- **Assembly Fetching**: Downloads assemblies from NCBI, with optional filtering for India-specific submissions.
- **Quality Assessment**: Uses QUAST to evaluate assembly quality.
- **Annotation**:
  - **Prokaryotic**: Prokka for rapid annotation.
  - **Eukaryotic**: MAKER with repeat modeling, SNAP training, and Augustus.
- **Completeness Evaluation**: BUSCO for annotation quality assessment.
- **Flexible Modes**:
  - Default: All India-specific assemblies.
  - Limited: User-defined number via `--num-assemblies`.
  - Specific: Assemblies from file via `--assembly-file`.
- **Multi-threading**: Configurable CPU usage.
- **Custom Filtering**: Extendable India-specific keyword list.

---

## Requirements

- **OS**: Linux (tested on Ubuntu)
- **Conda**: Miniconda or Anaconda
- **Internet Access**: For NCBI downloads & BUSCO datasets
- **Disk Space**: 10–20 GB

---

## Installation and Environment Setup
1. **Clone the Repository**:
```bash
git clone https://github.com/skbinfo/AquaaG.git
cd AquaaG
```
2. **Create Conda Environment**:
The pipeline uses a Conda environment defined in `environment.yml`. Create and activate it:
```bash
conda env create -f environment.yml
conda activate assembly_tool_env
```
3. **Verify Installation**:
Run test commands to ensure tools are installed:
```bash
maker --version
busco --version
prokka --version
```
## Configuration
The pipeline uses YAML configuration files for prokaryotic and eukaryotic modes. Edit these files before running.

- **Eukaryotic Config (`config.yaml`)**:
- Update `organism_query` (e.g., `"Saccharomyces cerevisiae[Organism]"`), `email`, `reference_genome_url`, `gff_url`.
- In `maker_params`, provide absolute paths to `protein_evidence` (FASTA file of homologous proteins), `est_evidence` (optional ESTs), and control files (`maker_opts.ctl`, etc.).
- Set `augustus_species` if using a pre-trained Augustus model (e.g., `"saccharomyces_cerevisiae_S288C"`).
- Example snippet:
 ```yaml
 organism_query: "Saccharomyces cerevisiae[Organism]"
 maker_params:
   protein_evidence: "/path/to/input/sc_protein.fasta"
   augustus_species: "saccharomyces_cerevisiae_S288C"
   cpus: 16
 ```

- **Prokaryotic Config (`prokaryotic_config.yaml`)**:
- Update `organism_query` (e.g., `"Salmonella enterica[Organism]"`), `email`, `reference_genome_url`, `gff_url`.
- Set `prokka_kingdom` (e.g., `"Bacteria"`).
- Example snippet:
 ```yaml
 organism_query: "Salmonella enterica[Organism]"
 prokka_kingdom: "Bacteria"
 cpus: 8
 ```

- **User File Updates**:
- **Evidence Files for MAKER (Eukaryotic Mode)**:
 - `protein_evidence`: Provide an absolute path to a FASTA file of protein sequences from related organisms (required). Example: Download Saccharomyces cerevisiae proteins from UniProt.
 - `est_evidence`: Provide an absolute path to a FASTA file of EST/cDNA sequences (optional; comment out with `est_evidence: ""` if not used).
 - Place files in the `input/` directory (e.g., `input/sc_protein.fasta`).
- **MAKER Control Files**:
 - Edit `maker_opts.ctl`: Set `organism_type=eukaryotic`, `est2genome=0`, `protein2genome=1`, `keep_preds=1`.
 - Edit `maker_bopts.ctl`: Adjust BLAST parameters if needed (e.g., `blastp_evalue=1e-06`).
 - Edit `maker_evm.ctl`: Set weights for evidence integration (e.g., `evmtrans=10`, `evmab:snap=10`).
 - Provide absolute paths to these files in `config.yaml`.
- **Augustus Species Model**:
 - If using Augustus, set `augustus_species` in `config.yaml` (e.g., `"saccharomyces_cerevisiae_S288C"`).
 - Download or train a model using `augustus --species=help` and ensure it’s accessible to MAKER.
- **BUSCO Lineage**:
 - Set `busco_lineage` to match your organism (e.g., `fungi_odb10` for fungi, `bacteria_odb10` for bacteria).
- **Assembly File (e.g., `genome.txt`)**:
 - Create a text file with NCBI assembly accession IDs (e.g., `GCA_965239625.1`), one per line.
 - Example `genome.txt`:
   ```
   GCA_965239625.1
   GCA_123456789.1
   ```
 - Use with the `--assembly-file` flag to process specific assemblies.

## Running the Pipeline
The pipeline supports three modes for selecting assemblies to process:

1. **Default Mode (All India-Specific Assemblies)**:
- Downloads and processes all assemblies from India for the organism specified in `config.yaml`.
- Command:
  ```
  python AquaaG.py -c config.yaml -o output_directory
  ```
- Behavior: Queries NCBI for assemblies matching `organism_query` (e.g., `"Saccharomyces cerevisiae[Organism]"`), filters for India-specific submitters (based on `SubmitterOrganization`), and processes all matching assemblies.

2. **Limited India-Specific Assemblies (`--num-assemblies`)**:
- Limits the number of India-specific assemblies to process (e.g., 1 or 2).
- Command:
  ```
  python AquaaG.py -c config.yaml -o output_directory --num-assemblies 1
  ```
  or
  ```
  python AquaaG.py -c config.yaml -o output_directory --num-assemblies 2
  ```
- Behavior: Filters for India-specific assemblies but processes only the first N assemblies (sorted by NCBI default order) as specified by `--num-assemblies`.

3. **Specific Assemblies (`--assembly-file`)**:
- Processes assemblies listed in a user-provided file (e.g., `genome.txt`), bypassing the India-specific filter.
- Command:
  ```
  python AquaaG.py -c config.yaml -o output_directory --assembly-file genome.txt
  ```
- Behavior: Reads accession IDs from `genome.txt` (e.g., `GCA_965239625.1`) and processes only those assemblies, regardless of their submitter location.

4. **Prokaryotic Mode**:
- Use `prokaryotic_config.yaml` for prokaryotic organisms (e.g., Salmonella enterica):
  ```
  python AquaaG.py -c prokaryotic_config.yaml -o output_directory --num-assemblies 1
  ```
- Supports `--assembly-file` or `--num-assemblies` as above.

5. **Outputs**:
- `output_directory/reference_data/`: Reference genome and GFF files.
- `output_directory/assemblies/`: Downloaded and decompressed FASTA files (e.g., `GCA_965239625.1.fna`).
- `output_directory/quast_output/`: QUAST reports (e.g., `report.txt`, `report.html`).
- `output_directory/prokka_output/` (PK mode): Prokka annotations per accession.
- `output_directory/maker_output/` (EK mode): MAKER results, including masked genomes, repeat libraries, SNAP models, and annotations (e.g., `GCA_965239625.1_r2.all.maker.proteins.fasta`).
- `output_directory/busco_output/`: BUSCO summaries (e.g., `short_summary.specific.fungi_odb10.GCA_965239625.1.txt`).
- `output_directory/pipeline.log`: Detailed log file for debugging.

## Customizing India-Specific Keywords
The pipeline filters for India-specific assemblies using a predefined list of keywords (e.g., cities like "Delhi", "Mumbai", and institutions like "IIT", "ICMR") in the `filter_rows_by_cities` function in `AquaaG.py`. To add more India-specific keywords:

1. **Locate the Keyword List**:
- Open `AquaaG.py` and find the `filter_rows_by_cities` function.
- The `indian_cities` list contains keywords:
  ```python
  indian_cities = [
      "IND", "Indian", "india", "India", "Agartala", "Ahmedabad", "Aizawl", ..., "IIT", "NIT", "IISC", ...
  ]
  ```

2. **Add New Keywords**:
- Append new keywords to the `indian_cities` list. For example, to add "Tamil Nadu" and "NCCS":
  ```python
  indian_cities = [
      "IND", "Indian", "india", "India", "Agartala", "Ahmedabad", "Aizawl", ...,
      "Tamil Nadu", "NCCS"  # Add new keywords here
  ]
  ```
- Ensure keywords are unique and relevant to `SubmitterOrganization` fields in NCBI Assembly records.

3. **Save and Test**:
- Save `AquaaG.py`.
- Test the updated filter:
  ```
  python AquaaG.py -c config.yaml -o test_output
  ```
- Check `test_output/pipeline.log` for the filtered assemblies:
  ```
  Found X assemblies matching the location 'India': [list of accession IDs]
  ```

4. **Best Practices**:
- Avoid overly broad keywords (e.g., "Institute") to prevent false positives.
- Verify new keywords by querying NCBI manually:
  ```
  esearch -db assembly -query "Saccharomyces cerevisiae[Organism]" | esummary | xtract -pattern DocumentSummary -element AssemblyAccession,SubmitterOrganization
  ```
- Ensure keywords match `SubmitterOrganization` exactly (case-insensitive).

## Troubleshooting
- **No India-Specific Assemblies**:
- If the pipeline exits with "No assemblies from India matched the query," check `pipeline.log` for available assemblies.
- Use `--assembly-file genome.txt` to process specific assemblies or verify `organism_query`.
- **MAKER/BUSCO Failures**:
- Ensure `protein_evidence` and `est_evidence` (if used) are valid FASTA files.
- Check `maker_output/[accession]/round2/[accession]_r2.maker.output/maker_opts.log` for MAKER errors.
- Verify BUSCO lineage matches the organism (e.g., `fungi_odb10` for Saccharomyces cerevisiae).
- **RepeatModeler Errors**:
- Run `BuildDatabase -version` to confirm RepeatModeler compatibility (should be 2.0.4).
- **Conda Issues**:
- If tools are missing, recreate the environment: `conda env create -f environment.yml`.
- **Low CPU Usage**:
- Update `cpus` in `config.yaml` or `prokaryotic_config.yaml` to match your system (e.g., 16 for a 16-core machine).

## AND OTHERS
- **Contributors**: Special thanks to all contributors who have provided feedback, tested the pipeline, or suggested improvements. Contributions are welcome via pull requests!
- **Tools and Libraries**: Gratitude to the developers of MAKER, Prokka, QUAST, BUSCO, RepeatModeler, RepeatMasker, pandas, BioPython, and other dependencies that power AquaaG.
- **Community**: Inspired by the bioinformatics community’s open-source efforts to streamline genomic analysis.



# 🇮🇳 download.py  
**Download Reference Genomes + India-Specific Assemblies from NCBI**

`download.py` is a Python tool that automates downloading **reference/representative genomes** and **India-specific assemblies** from the **NCBI Assembly database** for any given kingdom (e.g., fungi, bacteria, viruses, cyanobacteria).  

---

##  Features
- ✅ Detects **India-specific assemblies** using submitter info (cities, institutes, keywords).  
- ✅ Downloads **reference or representative genomes** (FASTA + GFF).  
- ✅ Fetches **India-specific assemblies** for species.  
- ✅ Creates a **metadata table** (`india_assemblies_metadata.tsv`) with accession, submitter, FTP path, and saved filename.  
- ✅ Parallelized downloads with logging and error handling.  

---
##  Usage

```bash
#  Process all species in a kingdom
python download.py -s fungi

#  Process only a fixed number of species
python download.py -s fungi --num 5

#  Example: Cyanobacteria, first 3 species
python download.py -s Cyanobacteriota --num 3
```

## Contributions
Contributions are welcome! Submit pull requests for bug fixes or new features. For major changes, open an issue first.

## Acknowledgments
- Built with tools like MAKER, Prokka, QUAST, BUSCO, and RepeatModeler.
- Inspired by bioinformatics pipelines for automated annotation.

For questions, open an issue on GitHub.
