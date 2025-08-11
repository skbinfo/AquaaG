# AquaaG: Automated Genomic Assembly and Annotation Pipeline

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

## Installation

### 1. Clone the Repository
```bash
git clone https://github.com/yourusername/AquaaG.git
cd AquaaG

2. Create Conda Environment

conda env create -f environment.yml
conda activate assembly_tool_env

3. Verify Installation

maker --version
busco --version
prokka --version

Configuration
Eukaryotic Config (config.yaml)

Example:

organism_query: "Saccharomyces cerevisiae[Organism]"
maker_params:
  protein_evidence: "/path/to/input/sc_protein.fasta"
  augustus_species: "saccharomyces_cerevisiae_S288C"
  cpus: 16

Prokaryotic Config (prokaryotic_config.yaml)

Example:

organism_query: "Salmonella enterica[Organism]"
prokka_kingdom: "Bacteria"
cpus: 8

Running the Pipeline
Default Mode

python assembly_tool.py -c config.yaml -o output_directory

Limited Assemblies

python assembly_tool.py -c config.yaml -o output_directory --num-assemblies 2

Specific Assemblies

python assembly_tool.py -c config.yaml -o output_directory --assembly-file genome.txt

Prokaryotic Mode

python assembly_tool.py -c prokaryotic_config.yaml -o output_directory --num-assemblies 1

Output Structure

    reference_data/: Reference genome and GFF

    assemblies/: Downloaded assemblies

    quast_output/: QUAST reports

    prokka_output/: Prokka results (PK mode)

    maker_output/: MAKER results (EK mode)

    busco_output/: BUSCO summaries

    pipeline.log: Run logs

Customizing India-Specific Keywords

To edit location filtering:

    Open assembly_tool.py

    Find filter_rows_by_cities() and indian_cities list

    Add your keywords, e.g.:

indian_cities = [
    "IND", "Indian", "India", "Delhi", "Mumbai",
    "Tamil Nadu", "NCCS"
]

    Save and rerun.

Troubleshooting

    No Assemblies Found: Check pipeline.log or use --assembly-file.

    MAKER/BUSCO Failures: Verify FASTA inputs & lineage settings.

    RepeatModeler Errors: Check version compatibility.

    Low CPU Usage: Increase cpus in config file.

Contributions

PRs welcome! Open issues for bugs/features.
Acknowledgements

Built with MAKER, Prokka, QUAST, BUSCO, RepeatModeler.
