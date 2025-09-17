
# NCBI Assembly Metadata Harvester

Fetch assembly-level metadata from NCBI using the **NCBI Datasets CLI**, enrich it with **Taxonomy Lineage** (ETE3) and **BioSample** fields (Entrez), append to a **single TSV**, and rebuild a **filtered TSV** with high-quality geography/biome entries.

> **Script**: `get_metadata_update.py`  
> **Language**: Python 3.8+  
> **OS**: Linux/macOS/Windows (via Conda)

---

## Features

- **Resume-safe ingestion**  
  Skips assemblies already present in the output (and in an optional `--resume` file).

- **Auto-enrichment**  
  Adds `Lineage` (ETE3), plus BioSample‐derived `Location`, `BiomeDistribution`, `Latitude`, `Longitude` (Entrez).

- **Unified placeholders**  
  Any missing value is **normalized** to `"Unknown"` (including `Latitude` and `Longitude`).

- **Filtered view**  
  Rebuilds `filtered_<output>.tsv` at the end using:
  - `BiomeDistribution != "Unknown"`
  - `Latitude != "Unknown"`
  - `Longitude != "Unknown"`

- **Retry modes** (healing):
  - `--retry-missing` (heal after a normal run)
  - `--retry-only` (heal the existing file and exit)
  - `--retry-columns biome,latlon,location,lineage`
  - `--retry-changed-only` (log only updated rows)

---

## What This Script Produces

- **Main output TSV** (you choose the filename), containing a header and one line per assembly with the following columns (fixed order):

```
Assembly Accession, Organism Name, Organism Common Name, Organism Tax ID, Lineage,
Assembly Level, BioProject Accession, BioSample Accession, GC Percent,
Total Sequence Length, Sequencing Technology, Release Date, Collection Date,
BioSample Description, Location, BiomeDistribution, Latitude, Longitude
```

- **Filtered output TSV**: `filtered_<output>.tsv`  
  Rebuilt from scratch at the end of each run with only high-quality rows (see [How the Filtered File Works](#how-the-filtered-file-works)).

---

## Installation

### Install Conda (or Mamba)

Choose one of the following:

**Miniconda (recommended)**
- Download: <https://docs.conda.io/en/latest/miniconda.html>  
- Install (Linux/macOS, bash example):
  ```bash
  # Run the installer you downloaded, e.g.:
  bash Miniconda3-latest-Linux-x86_64.sh
  # Follow the prompts and allow it to initialize your shell
  ```

**Mambaforge (Conda + Mamba)**
- Download: <https://github.com/conda-forge/miniforge#miniforge3-and-mambaforge>  
- Install similarly to Miniconda.

**Windows**
- Use the official Miniconda installer (.exe) and open the “Anaconda Prompt” for commands below.

> Verify installation:
```bash
conda --version
# (optional if installed) 
mamba --version
```

### Create the Project Environment

```bash
# Create an isolated environment (Python version is flexible; 3.10 works well)
conda create -n ncbi-meta python=3.10 -y

# Activate it
conda activate ncbi-meta
```

Install Python dependencies:

```bash
# Using pip inside the env
pip install ete3 biopython
```

> **Optional (Conda alternative)**  
> You can also install via Conda:  
> `conda install -c conda-forge ete3 biopython -y`

### Install NCBI Datasets CLI

You need **both** `datasets` and `dataformat` on your PATH. The **easiest** way is via Conda:

```bash
conda install -c conda-forge ncbi-datasets-cli -y
```

Verify:

```bash
datasets --version
dataformat --help
```

macOS users can alternatively use Homebrew:

```bash
brew install ncbi-datasets
```

---

## Environment Variables

NCBI requires you to identify yourself via email, and an API key improves rate limits.

**Linux/macOS (bash/zsh):**
```bash
export ENTREZ_EMAIL="you@example.org"
export NCBI_API_KEY="YOUR_KEY"   # optional
```

**Windows (PowerShell):**
```powershell
$env:ENTREZ_EMAIL="you@example.org"
$env:NCBI_API_KEY="YOUR_KEY"   # optional
```

You can also pass `--email you@example.org` on the command line; it overrides `ENTREZ_EMAIL`.

---

## Input & Output Formats

### Input (`input_file`)
- A **TSV (or whitespace-separated)** file with **Assembly Accessions in the 1st column**.
- By default, we assume the file has a **header**.  
  If your file does **not** have a header, pass `--no-input-header`.

Example:
```
Assembly Accession
GCF_000001405.40
GCA_000146045.2
...
```

### Output (`output_file`)
- A **single TSV** to which the script **appends** enriched rows.  
- It will create the header if the file is missing/empty.

### Filtered Output
- Rebuilt at the end of each run as:
```
filtered_<output_file_name>
```
- Contains **only** rows where `BiomeDistribution`, `Latitude`, and `Longitude` are **not** `"Unknown"`.

---

## Getting a Full List of Accessions (“Download the whole dataset”)

This script **does not** crawl NCBI; it processes whatever accessions you provide.  
To **generate** a full list of accessions (e.g., for a taxon or bioproject), use **NCBI Datasets**:

### A) Accessions by **Taxon** (replace `TAXON_ID`)
```bash
datasets summary genome taxon TAXON_ID --as-json-lines | dataformat tsv genome --fields accession > ids_genomes.tsv

# Add a header for clarity
{ echo "Assembly Accession"; cat ids_genomes.tsv; } > ids_with_header.tsv
```

### B) Accessions by **BioProject** (replace `PRJNAxxxxxx`)
```bash
datasets summary genome bioproject PRJNA123456 --as-json-lines | dataformat tsv genome --fields accession > ids_genomes.tsv

{ echo "Assembly Accession"; cat ids_genomes.tsv; } > ids_with_header.tsv
```

You can merge/clean multiple lists before running the script.

---

## Usage

> Always activate your environment first:
```bash
conda activate ncbi-meta
```

### Fresh Full Run

```bash
nohup python get_metadata_update.py   ids_with_header.tsv   metadata_all.tsv   --email you@example.org   --delay 0.35   > run.log 2>&1 &
```

- Reads IDs from `ids_with_header.tsv`
- Skips IDs already present in `metadata_all.tsv`
- Appends enriched rows to `metadata_all.tsv`
- Rebuilds `filtered_metadata_all.tsv` at the end

### Resume a Partial Run

```bash
nohup python get_metadata_update.py   ids_with_header.tsv   metadata_all.tsv   --resume metadata_all.tsv   --email you@example.org   --delay 0.35   > run_resume.log 2>&1 &
```

Using `--resume metadata_all.tsv` simply tells the script to also consider those accessions as “already processed”.

### Retry Modes (Healing Incomplete Rows)

**Why?**  
If the network dropped, Entrez rate-limited, or BioSample data was temporarily unavailable, some fields may be `"Unknown"`.  
Use **retry** to fill them **later** without reprocessing everything.

#### A) Heal **after** a normal run
```bash
python get_metadata_update.py   ids_with_header.tsv   metadata_all.tsv   --retry-missing   --retry-columns biome,latlon,location,lineage   --retry-changed-only   --email you@example.org   --delay 0.35
```

#### B) Heal **only** the existing file and exit
```bash
python get_metadata_update.py   ids_with_header.tsv   metadata_all.tsv   --retry-only   --retry-columns biome,latlon   --retry-changed-only   --email you@example.org   --delay 0.35
```

**What retry tries to fix** (based on `--retry-columns`):
- `lineage` → `Lineage == "Unknown"`  
- `location` → `Location == "Unknown"`  
- `biome` → `BiomeDistribution == "Unknown"`  
- `latlon` → `Latitude == "Unknown"` or `Longitude == "Unknown"`

If the lookup succeeds, the row is updated **in place** and reported (when `--retry-changed-only` is set).

### Command Reference

```text
usage: get_metadata_update.py input_file output_file [options]

positional arguments:
  input_file            TSV with Assembly Accessions in the 1st column.
  output_file           Single metadata TSV (will be appended).

optional arguments:
  --resume TSV          Also skip accessions already present in this TSV.
  --no-input-header     Use if input_file has NO header row.
  --email EMAIL         Entrez email (overrides env ENTREZ_EMAIL).
  --delay FLOAT         Delay (seconds) between Entrez requests; default 0.35.

  --retry-missing       After normal processing, heal missing fields in output_file.
  --retry-only          Do not process new IDs; only heal output_file and exit.
  --retry-columns CSV   Which fields to heal: biome,latlon,location,lineage. Default: all.
  --retry-changed-only  During retry, log only rows that actually improved.
```

---

## How the Filtered File Works

The script rebuilds `filtered_<output>.tsv` **from scratch** at the end. A row is included **only if**:
- `BiomeDistribution != "Unknown"`
- `Latitude != "Unknown"`
- `Longitude != "Unknown"`

This gives you a clean subset suitable for geospatial/biome analyses.

---

## Placeholders and “Missing” Values

To keep column counts stable and TSVs tidy, **every** missing field is stored as `"Unknown"`:

- `Lineage`: `"Unknown"` if ETE3 lineage lookup fails.  
- `Location`: `"Unknown"` if BioSample lacks geographic info.  
- `BiomeDistribution`: `"Unknown"` if cannot be inferred from BioSample.  
- `Latitude`: `"Unknown"` if missing/unparsed.  
- `Longitude`: `"Unknown"` if missing/unparsed.

Retry modes try to replace `"Unknown"` with real values later.

---

## Performance & Operational Tips

- **Be gentle with NCBI**: keep `--delay` around `0.3–0.5s`.  
  Adding `NCBI_API_KEY` boosts your quotas.

- **Long runs**: use `nohup ... &` and direct logs to file.  
- **Selective healing**: `--retry-columns` (e.g., `latlon` only) reduces API calls.  
- **ETE3 first-run**: it will download taxonomy data on first use; allow network access.  
  You can pre-warm via Python:
  ```python
  from ete3 import NCBITaxa
  NCBITaxa().update_taxonomy_database()
  ```

---

## Troubleshooting

**`datasets` or `dataformat` not found**  
- Ensure NCBI Datasets CLI is installed and on PATH.  
  `conda install -c conda-forge ncbi-datasets-cli -y`

**ETE3 taxonomy errors / lineage = `"Unknown"`**  
- First run downloads taxonomy (SQLite DB) to your home (commonly `~/.etetoolkit/`).  
- If it fails, check your network or pre-warm as shown above.

**Many `"Unknown"` values**  
- Some BioSamples genuinely lack metadata.  
- Try a healing pass later: `--retry-only` with `--retry-columns biome,latlon,location,lineage`.

**HTTP 429 / Rate limits**  
- Increase `--delay`, set `NCBI_API_KEY`, and/or pause between runs.

**Header handling**  
- The script writes **its own header** if the output is empty.  
- Use `--no-input-header` if your input has **no** header row.

---

## Citations, Compliance & Fair Use

- Use of NCBI resources must comply with their policies:
  - Identify yourself: set **ENTREZ_EMAIL** (or `--email`).
  - Respect rate limits: use `--delay` and consider **NCBI_API_KEY**.
- If you publish results, cite:
  - **NCBI Datasets CLI**
  - **NCBI Entrez / E-utilities**
  - **ETE3 Toolkit** (for taxonomy lineage)

---

## License

Add your preferred license here (e.g., MIT). Ensure compliance with third-party tools’ licenses (NCBI Datasets CLI, ETE3).

---

### Quick Copy-Paste: End-to-End Example

```bash
# 0) Activate env
conda activate ncbi-meta

# 1) Build an IDs list (taxon example: 1117 is Pseudomonas)
datasets summary genome taxon 1117 --as-json-lines | dataformat tsv genome --fields accession > ids_genomes.tsv
{ echo "Assembly Accession"; cat ids_genomes.tsv; } > ids_with_header.tsv

# 2) Full run
nohup python get_metadata_update.py   ids_with_header.tsv   metadata_all.tsv   --email you@example.org   --delay 0.35   > run.log 2>&1 &

# 3) Heal only lat/lon & biome later
python get_metadata_update.py   ids_with_header.tsv   metadata_all.tsv   --retry-only   --retry-columns biome,latlon   --retry-changed-only   --email you@example.org   --delay 0.35
```
