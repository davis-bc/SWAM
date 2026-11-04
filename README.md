# SWAM — Surface Water AMR Monitoring

Snakemake workflow for end-to-end metagenomic analysis of antibiotic resistance in environmental bacterial communities.

## What it produces

SWAM starts from paired-end FASTQs and generates:

- **`fastp_summary.csv`** - per-sample read QC and Nonpareil coverage
- **`short_reads_output.csv`** - per-gene short-read AMR abundances
- **`ARGnet_short-reads_summary.tsv`** - discovery-mode ARGnet calls on host-filtered short reads
- **`assembly_qa.tsv`** - per-sample assembly size, N50, and read-mapping stats
- **`contig_summary.tsv`** - per-contig abundance, taxonomy, AMR, MGE, and molecule type
- **`AMR_abundance_summary.csv`** - per-sample AMR abundance plus additive and multiplicative risk scores
- **`mag_summary.tsv`** - per-MAG abundance and annotations
- **`ARGnet_contig_summary.tsv`** - discovery-mode ARGnet calls on contig proteins with contig context
- **`ARGnet_mag_summary.tsv`** - discovery-mode ARGnet calls on MAG proteins

It also keeps lower-level outputs under `out_dir/data/`, including assemblies, geNomad results, MobMess output, MAG files, and the marker support table `data/QAQC/markers_cpg.csv`.

---

## Pipeline at a glance

1. **Short reads** - fastp QC, host filtering, KMA AMR alignment, optional ARGnet discovery, SCG alignment, anthropogenic marker alignment, Nonpareil
2. **Contigs** - MEGAHIT assembly, geNomad classification, MobMess, Prodigal, AMRFinderPlus, optional ARGnet discovery, MobileElementFinder, MMseqs2 taxonomy, contig abundance
3. **MAGs** - MetaBAT2 binning, per-bin AMR/MGE annotation, optional ARGnet discovery, abundance, and optional GTDB-tk, CheckM2, and METABOLIC
4. **Summary** - sample-level AMR abundance and risk scoring from AMR, mobility, host context, and marker-derived exposure

[![SWAM rulegraph](docs/rulegraph.png)](docs/rulegraph.png)

---

## Installation

```bash
git clone https://github.com/davis-bc/SWAM.git
cd SWAM
conda install -c conda-forge -c bioconda "snakemake>=8" snakemake-executor-plugin-slurm
python -m pip install -e .        # installs the `swam` CLI (textual, ruamel.yaml)
```

A conda recipe is in `recipe/meta.yaml`. All workflow tools are installed automatically into per-rule conda environments on first run.

The supported entrypoint is the `swam` CLI:

| Command | What it does |
|---|---|
| `swam config [--profile slurm\|local]` | Terminal UI to toggle stages, set paths, and configure Slurm + Snakemake |
| `swam run [--profile slurm\|local]` | Stages the packaged workflow and runs Snakemake, then writes the run report |

Profiles are user-owned copies seeded on first use at `~/.config/swam/profiles/<profile>/config.yaml` (every TUI save keeps a backup in `~/.config/swam/backups/`). The workflow is staged into `~/.cache/swam/runtime/<version>/`; databases default to `~/.local/share/swam/dbs/` (set **Database dir** in the TUI to reuse an existing `dbs/`).

---

## Quick start

### 1. Configure

```bash
swam config --profile local     # or --profile slurm
```

| Side tab | Settings |
|---|---|
| **Stages** | Short reads / Contigs / MAGs toggles, AMR mode, geNomad splits, min read length, GTDB-tk/METABOLIC/CheckM2 skips |
| **Paths** | input FASTQ dir, output dir, database dir, GTDB-tk data, clean-reads and contigs dirs for partial runs |
| **Slurm** | account, partition, default mem/walltime, per-rule threads / mem / walltime / group size |
| **Snakemake** | jobs, cores, restarts, latency wait, conda frontend, keep-going / rerun-incomplete etc. |

Keys: `1`–`4` switch tabs, `Ctrl+S` save, `Ctrl+R` reload, `Q` quit. The orb bar shows how many samples were found in the input dir, and the Stages tab flags missing inputs for partial runs.

`amr_mode: operational` keeps the AMRFinderPlus-based outputs only. `amr_mode: discovery` keeps those outputs and adds ARGnet-based discovery summaries for short reads, contig proteins, and MAG proteins.

### 2. Run

```bash
swam run --profile local -n                       # dry run using the profile's in/out dirs
swam run --profile local --input /path/to/fastq --output /path/to/output
```

`--input/--output` override the profile; anything else is passed to Snakemake (`-n`, `--until`, `--forcerun <rule>`, `--config key=value`, ...). `--scheduler greedy` is always added (the default MILP scheduler needs the missing `cbc` solver).

### Test mode

Runs the full workflow on two bundled mock samples and mini databases (needs this source checkout):

```bash
swam run --profile local --test -n
swam run --profile local --test
```

Expected runtime is about 30-60 minutes on a laptop. Outputs go to `test/output/`.

### Production mode on SLURM

Set the Slurm account and partition with `swam config --profile slurm`, then submit the driver job:

```bash
#!/bin/bash
#SBATCH --job-name=swam
#SBATCH --time=7-00:00:00
#SBATCH --mem=8G
swam run --profile slurm --input /path/to/fastq_files --output /path/to/output
```

`bash run_swam.sh /path/to/fastq_files /path/to/output` still works and forwards to `swam run --profile slurm`. Each run writes:

- `{out_dir}/logs/run_status/run_report.txt`
- `{out_dir}/logs/run_status/sample_rule_status.tsv`

---

## Databases

Most databases are downloaded and prepared automatically in the database dir (`dbs_dir`, default `~/.local/share/swam/dbs/`) on first use, including:

- AMRFinderPlus databases
- ARGNet models (discovery mode only)
- human reference for host filtering
- anthropogenic markers (`pBI143`, `crAss001`)
- geNomad database
- UniRef50 MMseqs2 taxonomy database
- CheckM2 database
- METABOLIC

The database dir is organized by tool or reference bundle:

```text
dbs/
  amrfinderplus/       AMRFinderPlus protein database for contig/MAG annotation
  amrfinderplus_kma/   KMA-formatted AMR CDS database for short reads
  checkm2/             CheckM2 reference bundle
  genomad/             geNomad database bundle
  human/               host-filter reference genome
  markers/             pBI143 / crAss001 FASTAs and KMA index
  metabolic/           METABOLIC checkout
  scg/                 DIAMOND SCG database
  uniref50/            MMseqs2 UniRef50 taxonomy database
```

**GTDB-tk is the only manual setup step.**

```bash
conda run -n gtdbtk download-db.sh /path/to/gtdbtk_db
```

Then set **GTDB-tk data dir** on the TUI Paths tab (`gtdbtk_db`).

Optional stages can be disabled on the TUI Stages tab (profile `config:` block):

```yaml
skip_gtdbtk: False
skip_metabolic: False
skip_checkm2: False
genomad_splits: 1
```

Increase `genomad_splits` if geNomad needs less peak memory.

---

## Running on SLURM

SWAM ships with a single default Slurm profile (`src/swam/assets/config/slurm/config.yaml`), seeded into `~/.config/swam/profiles/slurm/` on first use. It enables `use-conda`, the greedy scheduler, retries, `keep-going`, and grouped execution for lightweight sample-parallel rules. `contigs`, `genomad`, and `bin` scale memory and walltime per retry attempt inside the rules, so only their threads are set in the profile.

---

## Running only part of the workflow

Toggle stages on the TUI Stages tab, or override per run with `--config`:

```yaml
run_short_reads: True
run_contigs: True
run_mags: True
```

| `run_short_reads` | `run_contigs` | `run_mags` | Result | Extra input needed |
|---|---|---|---|---|
| `True` | `True` | `True` | Full pipeline | None |
| `True` | `True` | `False` | Short reads + contigs | None |
| `True` | `False` | `False` | Short reads only | None |
| `False` | `True` | `True` | Contigs + MAGs | `clean_reads_dir` |
| `False` | `True` | `False` | Contigs only | `clean_reads_dir` |
| `False` | `False` | `True` | MAGs only | `clean_reads_dir`, `contigs_dir` |

When skipping upstream stages:

- `clean_reads_dir` must contain host-filtered paired FASTQs
- `contigs_dir` must contain `{sample}.contigs.fa`

If `run_short_reads: False`, contig abundance falls back to mean depth rather than fully normalized cpg because the SCG-based genome estimate is unavailable.

---

## Main outputs

| Path | Description |
|---|---|
| `{out_dir}/fastp_summary.csv` | Per-sample fastp metrics and Nonpareil coverage |
| `{out_dir}/short_reads_output.csv` | Short-read AMR gene abundances |
| `{out_dir}/assembly_qa.tsv` | Per-sample assembly QC metrics |
| `{out_dir}/contig_summary.tsv` | Per-contig annotations and abundance |
| `{out_dir}/AMR_abundance_summary.csv` | Per-sample AMR abundance and risk scores |
| `{out_dir}/mag_summary.tsv` | Per-MAG summary table |
| `{out_dir}/data/QAQC/markers_cpg.csv` | `pBI143` and `crAss001` copies per genome |

`AMR_abundance_summary.csv` contains:

| Column | Meaning |
|---|---|
| `sample` | Sample ID |
| `AMR_total_cpg` | Total AMR abundance across detected genes |
| `pBI143_cpg`, `crAss001_cpg` | Anthropogenic marker abundance |
| `E_exposure` | Exposure score from marker abundance |
| `R_mean` | Mean resistance-hazard score |
| `M_mean` | Mean mobility score |
| `H_mean` | Mean host/pathogenicity score |
| `amr_risk_additive_raw`, `amr_risk_multiplicative_raw` | Raw study-level risk scores |
| `amr_risk_additive`, `amr_risk_multiplicative` | Min-max normalized 0-100 scores |

---

## Abundance and risk scoring

All AMR and contig abundances are reported as **copies per genome (cpg)**:

```text
n_genomes = Σ(alignment_length / gene_length) / 40
cpg       = mean_depth / n_genomes
```

The 40-gene denominator comes from the bundled single-copy gene reference in `src/swam/assets/workflow/resources/SCGs_40_All.fasta`.

Risk scoring combines four components:

- **R** - resistance hazard from AMRFinderPlus class/subclass
- **M** - mobility from contig type and nearby MGE evidence
- **H** - host/pathogenicity from MMseqs2 taxonomy
- **E** - exposure from `pBI143` and `crAss001`

Both additive and multiplicative sample-level scores are reported.

---

## Repository layout

```text
pyproject.toml          `swam` package (CLI + TUI)
recipe/                 conda recipe
run_swam.sh        shim for `swam run --profile slurm`
src/swam/
  cli.py                `swam config` / `swam run`
  config.py             user profile paths, seeding, backups
  runtime.py            stages the workflow and runs Snakemake
  tui.py, tui.tcss      Textual configuration UI
  assets/
    config/             pipeline defaults (config.yaml) + slurm/ and local/ profile templates
    workflow/
      Snakefile         workflow entry point
      envs/             per-rule conda environments
      resources/        bundled reference files
      rules/            short-read, contig, MAG, and summary rules
      scripts/          Python/R scripts used by rules
docs/                   rule graph and session log
test/                   mock data, mini databases, reproducible test scripts
```
