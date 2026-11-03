# SWAM — Claude instructions

## Session continuity

**At the start of every session:** read `docs/session-log.md` to restore context on project history, current pipeline state, and any outstanding issues.

**At the end of every session** (or after completing significant work): append a new dated entry to `docs/session-log.md` summarising what was done, the current pipeline state, and any known issues or next steps. Commit the updated log.

## Overview

SWAM (Surface Water AMR Monitoring) is a Snakemake workflow for end-to-end metagenomic analysis of antibiotic resistance in environmental bacterial communities. Full pipeline: read QC + host filtering → short-read AMR/marker alignment → de novo assembly → contig classification (chromosome/plasmid/phage) → taxonomy (MMseqs2/UniRef50) → AMR annotation (AMRFinderPlus) → MGE annotation (MobileElementFinder) → contig abundance → unified AMR table → MAG binning → MAG characterisation (AMR, MGE, taxonomy, metabolism, quality, abundance).

## Running the Workflow

The supported entrypoint is the installed `swam` CLI (`python -m pip install -e .` from the repo). It stages the packaged workflow (`src/swam/assets/`) into `~/.cache/swam/runtime/<version>/` and runs Snakemake from there.

```bash
swam config --profile local|slurm      # Textual TUI: stages, paths, Slurm, Snakemake
swam run --profile local -n            # dry run (always adds --scheduler greedy)
swam run --profile local --input <fastq_dir> --output <out_dir>
swam run --profile slurm               # HPC driver job; dispatches via Slurm
swam run --profile local --forcerun <rule_name>     # extra args pass through to Snakemake
swam run --profile local --config run_mags=False -n # --config is merged, not replaced
bash run_swam.sh <in> <out> [args]     # shim -> swam run --profile slurm
```

> **`--scheduler greedy` is required** (the default MILP scheduler needs the missing `cbc` solver); `swam run` always adds it.

```bash
# Test mode — two synthetic mock samples (≤16 GB RAM, ≤4 cores); needs this checkout's test/
swam run --profile local --test -n     # dry run
swam run --profile local --test        # full run (~30–60 min), outputs in test/output/
```

Direct Snakemake (dev only): `snakemake -s src/swam/assets/workflow/Snakefile --scheduler greedy --use-conda --cores 4 --config test=True test_dir=test dbs_dir=dbs -n`.

## Configuration

- `src/swam/assets/config/config.yaml` — pipeline defaults, loaded by the Snakefile via `workflow.basedir`.
- `src/swam/assets/config/{slurm,local}/config.yaml` — Snakemake profile templates. Seeded on first use into `~/.config/swam/profiles/<profile>/config.yaml` (user-owned; TUI saves back up to `~/.config/swam/backups/`). Each profile has a `config:` block that overrides the defaults.
- Pipeline keys: `in_dir`, `out_dir`, `dbs_dir`, `gtdbtk_db`, `clean_reads_dir`, `contigs_dir`, `run_short_reads`, `run_contigs`, `run_mags`, `amr_mode` (`operational`/`discovery`), `skip_gtdbtk`, `skip_metabolic`, `skip_checkm2`, `genomad_splits`, `min_read_length`. Internal: `test`, `test_dir`.
- `dbs_dir` empty → `swam run` uses `~/.local/share/swam/dbs` (outside the runtime cache, which is wiped when packaged assets change). Raw Snakemake falls back to `<runtime root>/dbs`.
- Run path precedence: `--input/--output` > profile `config.in_dir/out_dir` > `./input`, `./output`.

## Architecture

### CLI package (`src/swam/`)

| File | Responsibility |
|------|----------------|
| `cli.py` | argparse entrypoint: `config`, `run` |
| `config.py` | XDG profile/cache/data paths, profile seeding, empty-section pruning, backups |
| `runtime.py` | fingerprinted staging of `assets/workflow` + `assets/config/config.yaml`, `--config` merging, driver log, `build_run_report.py` |
| `tui.py` + `tui.tcss` | Textual editor; OSRS-inspired `swam-osrs` theme (+ `classic`); edits only the profile YAML via ruamel round-trip |

### Rule files (`src/swam/assets/workflow/rules/`, loaded in `Snakefile` in this order)

| File | Responsibility |
|------|----------------|
| `rules/common.smk` | Sample discovery, test-mode path overrides, `res()` resource helper |
| `rules/short_reads.smk` | DB init, fastp QC, host filtering, KMA+DIAMOND alignments, Nonpareil, short-reads summary |
| `rules/contigs.smk` | MEGAHIT assembly, geNomad, MobMess, Prodigal, AMRFinderPlus (contigs), MobileElementFinder, MMseqs2 taxonomy, contig abundance, contig summary |
| `rules/mags.smk` | MetaBAT2 binning (checkpoint), per-bin Prodigal/AMR/MGE/taxonomy/metabolism, CheckM2 QC, CoverM abundance |
| `rules/summary.smk` | `amr_unified` — cross-stage: merges short-reads + contig AMR into `AMR_unified.csv` + `AMR_abundance_summary.csv` |

`localrules` (head-node only, not submitted to cluster): `initiate_dbs`, `short_reads_summary`, `init_genomad`, `init_mmseqs_db`, `init_mge_tool`, `contig_summary`, `amr_unified`.

### Pipeline stages (dependency order)

1. **`initiate_dbs`** — one-time: downloads AMRFinderPlus DB + human reference, indexes with KMA/DIAMOND. Sentinel `.done.txt` outputs.
2. **`short_reads`** — per-sample: fastp → minimap2 host removal → KMA (AMR + markers) → DIAMOND (SCGs) → Nonpareil.
3. **`short_reads_summary`** — all-sample aggregation via R script (`short_reads_processing.R`): outputs `fastp_summary.csv`, top-level `short_reads_output.csv`, and `data/QAQC/markers_cpg.csv`.
4. **`contigs`** — MEGAHIT assembly (`--presets meta-large`, min contig 1000 bp); contig headers renamed `{sample}-<original>` in-place.
5. **`init_genomad`** + **`genomad`** — downloads geNomad DB on first run; classifies contigs into chromosome/plasmid/phage (`--relaxed --cleanup`).
6. **`mobmess`** — maps reads to plasmid contigs, calls `detect_circular_contigs.py`, runs MobMess to infer plasmid systems. Has `|| true` fallback when all plasmids are singleton components.
7. **`init_mmseqs_db`** + **`mmseqs_taxonomy`** — builds MMseqs2 index from pre-built UniRef50 DB; runs `easy-taxonomy` per sample (`--lca-ranks superkingdom,...,species`).
8. **`prodigal`** — shared ORF prediction per sample (`.faa` + `.gff` consumed by both `contig_amr` and `contig_summary`).
9. **`contig_amr`** (`init_mge_tool` + **`mge_annotation`**) — AMRFinderPlus protein mode on Prodigal ORFs; MobileElementFinder on contigs.
10. **`contig_abundance`** — minimap2 (`-ax sr`) → samtools sort/index/depth → per-contig cpg via `contig_abundance.py`.
11. **`contig_summary`** — all-sample aggregation via `contig_summary.py`: joins geNomad + MMseqs2 + AMR + MGE + abundance into `contig_summary.tsv`.
12. **`amr_unified`** (`summary.smk`) — `amr_unified.py`: outer-joins `short_reads_output.csv` + `contig_summary.tsv` (AMR rows) on `sample × gene_symbol`. Outputs `AMR_unified.csv` + `AMR_abundance_summary.csv`.
13. **`bin`** — MetaBAT2 (`jgi_summarize_bam_contig_depths` → `metabat2`); reuses contig BAM from step 10.
14. **`mag_prodigal`** (checkpoint) — Prodigal per bin; dynamically discovers bins after `bin` completes.
15. **`mag_amr`** / **`mag_mge`** — AMRFinderPlus + MobileElementFinder per bin; `2>/dev/null || true` handles empty bins.
16. **`mag_qc`** — CheckM2 per sample (skipped when `skip_checkm2: True`).
17. **`mag_taxonomy`** — GTDB-tk `classify_wf` per sample (skipped when `skip_gtdbtk: True`).
18. **`mag_metabolism`** — METABOLIC-G.pl per sample (skipped when `skip_metabolic: True`).
19. **`mag_abundance`** — CoverM `genome` mode, trimmed-mean, `--min-covered-fraction 0.1`.

### Sample discovery (`rules/common.smk`)

Globs `in_dir` for `*R1*.fastq*` / `*_1.fastq*` patterns and pairs with R2 counterparts. The `sample` wildcard is constrained to `[^/]+` in `Snakefile`.

### Conda environments (`src/swam/assets/workflow/envs/`)

| File | Tools |
|------|-------|
| `shortreads.yaml` | fastp, kma, nonpareil, diamond, minimap2, samtools, pigz |
| `contigs.yaml` | megahit, minimap2, samtools, prodigal, mmseqs2, blast, pysam, MobMess + PlasX (pip from GitHub) |
| `genomad.yaml` | genomad |
| `Renv.yaml` | R, tidyverse, Nonpareil |
| `mge.yaml` | mobilefinder (MobileElementFinder) |
| `mags.yaml` | metabat2, coverm, prodigal, ncbi-amrfinderplus, mobilefinder, gtdbtk |
| `metabolic.yaml` | METABOLIC (separate env — Perl dependency) |
| `checkm2.yaml` | checkm2 |
| `gtdbtk.yaml` | gtdbtk |

**MobMess/PlasX** are installed via pip from GitHub inside `contigs.yaml` — internet access required on first `--use-conda` setup.

### Key scripts (`src/swam/assets/workflow/scripts/`)

| Script | Interface | Description |
|--------|-----------|-------------|
| `short_reads_processing.R` | Snakemake R (`snakemake@input`) | Normalises AMR + marker depth to cpg; outputs top-level `short_reads_output.csv` plus backend `markers_cpg.csv` |
| `contig_abundance.py` | Snakemake Python | Computes per-contig cpg from BAM depth + SCG-derived `n_genomes` |
| `contig_summary.py` | Snakemake Python | Joins geNomad, MMseqs2, AMR, MGE, abundance per sample into long-format TSV |
| `amr_unified.py` | Snakemake Python | Outer-joins short-reads + contig AMR on `sample × gene_symbol`; cpg priority: contig > short-reads; produces `evidence` column (`both`/`short_reads_only`/`contigs_only`) |
| `detect_circular_contigs.py` | Standalone (pysam) | Yu et al. 2024 circularity: RF-oriented reads spanning ≥ `contig_length − 3×median_insert` |

## Testing

```bash
swam run --profile local --test -n     # dry run
swam run --profile local --test        # full end-to-end (~30–60 min, ≤16 GB RAM)
```

- Mock FASTQ pairs in `test/data/` (mock1, mock2): biologically realistic reads simulated from 6 real NCBI reference genomes (E. coli, E. faecalis, Salmonella, Campylobacter, Carjivirus phage, Influenza A). ~280k–395k read pairs each.
- Mini databases in `test/dbs/`: AMRFinderPlus subset (49 real genes + 11 synthetic), UniRef50 subset, 40 SCGs, anthropogenic markers, mini human genome.
- `test/` is **not** packaged; `--test` auto-detects it for editable installs (else pass `--test-dir`). Rules read it through `_TEST_DIR` in `common.smk`.
- `test/scripts/generate_mock_data.py` documents reproducible mock data generation.
- Test mode sets `skip_gtdbtk=True`, `skip_metabolic=True`, `skip_checkm2=True`, and `genomad_splits=4` automatically.
- `test/output/` is gitignored.
- TUI check: drive `ConfigEditorApp` headlessly with `App.run_test()` (Pilot), save, and diff the profile YAML; `app.save_screenshot()` for visuals.

## Key Conventions

- **`--scheduler greedy`**: always pass this — MILP scheduler requires `cbc` solver which is not installed.
- **`res()` helper** (`common.smk`): `res(production_value, test_value=4000)` — all resource values use this pattern so test runs are lightweight. Both `threads:` (top-level directive) and `resources: threads` are declared; shell commands reference `{resources.threads}`.
- **`threads:` must be top-level**: `threads: lambda wc: res(32, 4)` — NOT inside `resources:`. Snakemake only uses the top-level `threads:` directive for core-based scheduling.
- **Resource declarations**: all compute-intensive rules declare `mem_mb`, `threads`, and `time` under both `threads:` and `resources:`.
- **Temporary files**: `mktemp -d` inside the output directory keeps temp files on the same filesystem (atomic moves). Temp dirs are cleaned up in shell blocks.
- **No `protected()` outputs**: previously used on clean reads and contigs but removed — causes `ProtectedOutputException` on re-runs after failures.
- **`localrules`**: these rules run on the head node, not submitted to the cluster. Currently: `initiate_dbs`, `short_reads_summary`, `init_genomad`, `init_mmseqs_db`, `init_mge_tool`, `contig_summary`, `amr_unified`.
- **cpg normalisation**: `n_genomes = Σ(alignment_length / gene_length) / 40` across the 40 SCGs (DIAMOND). `cpg = mean_depth / n_genomes`. The same `n_genomes` is shared between short-read AMR and contig-level abundance.
- **AMR_unified cpg priority**: contig cpg preferred when assembled (higher specificity); short-reads cpg used otherwise (captures low-depth genes that don't assemble).
- **MobMess fallback**: `mobmess systems ... || true` + `touch` sentinel — tool crashes when all plasmid contigs are singleton components (known edge case).
- **AMRFinder on small bins**: `amrfinder ... 2>/dev/null || true` — tool fails GFF validation on small bins; graceful skip is intentional.
- **contig BAM reuse**: the BAM from `contig_abundance` (minimap2 all-contigs) is reused by `bin` (MetaBAT2) to avoid a redundant mapping.
- **Prodigal shared output**: the `prodigal` rule in `contigs.smk` outputs `{sample}.faa` + `{sample}.gff`; both `contig_amr` and `contig_summary` consume these — do not run Prodigal redundantly per tool.
- **MGE output format**: MobileElementFinder outputs CSV (comma-separated) with comment lines (`#date:`, etc.) — parse with `pd.read_csv(sep=",", comment="#")`. Contig IDs have extra metadata appended; strip to first token.
- **Profiles**: one Slurm and one local profile template under `src/swam/assets/config/`; the old `small-batch` / `large-batch` split has been retired.
- **`_REPO` / `_TEST_DIR` / `_DBS_DIR`** (`common.smk`): `_REPO = dirname(workflow.basedir)` (runtime root); test assets and DBs come from `test_dir` / `dbs_dir` config when set. Never derive paths from cwd.
- **Repeated lists**: the TUI's `R1_GLOBS`, `LOCAL_RULES`, and `RETRY_SCALED_RULES` (`tui.py`) mirror sample discovery in `common.smk`, `localrules:` in the Snakefile, and the attempt-scaled rules in `contigs.smk`/`mags.smk` — update together.
- **Logging layout**: per-rule logs live under `{output_dir}/logs/` (`logs/<rule>/<sample>.log` for sample rules, `logs/<rule>.log` for global rules). Run-status summaries go to `{output_dir}/logs/run_status/`.
- **`.gitignore`**: `input/`, `output/`, `slurm/`, `dbs/`, `.snakemake/`, `run_snakemake.sh` are excluded.
- **New rule files**: MAG rules in `rules/mags.smk`, cross-stage summary rules in `rules/summary.smk` (under `src/swam/assets/workflow/`). New files must match `package-data` globs in `pyproject.toml`. Add terminal outputs to `all_targets()` in `Snakefile`.

## Claude Code notes
- Session log: `docs/session-log.md` (read the latest entry at session start)
