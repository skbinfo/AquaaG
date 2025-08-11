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
