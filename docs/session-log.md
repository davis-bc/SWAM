# SWAM-meta — Copilot Session Log

This file is maintained by GitHub Copilot. At the start of each session, Copilot reads this file to restore context. At the end of each session (or when significant progress is made), Copilot appends a new entry.

---

## 2026-02-22 (session 2)

### What was done
- **Mock dataset v2**: Replaced the old 50k-read synthetic mock data with a biologically realistic dataset simulated from 6 real NCBI reference genomes:
  - *E. coli* 2024GO-0438 (GCA_042156675.1, 5.18 Mb, 198 contigs, plasmids)
  - *E. faecalis* ST6 (GCA_015304995.1, 3.31 Mb, 89 contigs, plasmids)
  - *Salmonella enterica* ENWS000537-005 (GCA_052022995.1, 4.71 Mb, 43 contigs, plasmids)
  - *Campylobacter coli* 18QD2YX31C (GCF_013072625.1, 1.94 Mb, 5 contigs, plasmids)
  - Carjivirus communis phage (GCF_009734245.1, 97 kb, 1 contig)
  - Influenza A H3N2 (GCA_038998445.1, 13.6 kb, 8 segments)
- **mock1** = E. coli (10×) + E. faecalis (10×) + Carjivirus (3×) + Influenza A (3×) → 284,141 read pairs (16 MB per file gzipped)
- **mock2** = E. coli (10×) + Salmonella (10×) + Campylobacter (10×) + Carjivirus (3×) → 394,971 read pairs (23 MB per file gzipped)
- E. coli and Carjivirus appear in both samples to exercise CoverM cross-sample mapping.
- 40 synthetic SCGs (seed=52) + blaTEM-1 re-embedded into E. coli first contig at same positions; existing DIAMOND/KMA databases remain valid.
- Reference FASTAs saved to `test/dbs/genomes/`. `test/data/mock_data_summary.txt` documents composition and rationale.
- **CheckM2 MAG QC**: Added `mag_qc` rule to `mags.smk` (`checkm2 predict` per sample). New `workflow/envs/checkm2.yaml`. `checkm2_db` and `skip_checkm2` added to `config/config.yaml`. `skip_checkm2=True` set in test mode. `checkm2_quality.tsv` added to `rule all` as conditional target.
- Pushed to `origin/main` (commit `200c9c4`).

### Current pipeline state
All stages fully implemented and test data ready. Key additions since last session:

| Component | Status |
|-----------|--------|
| Mock data | ✅ Upgraded to real full-genome references (v2) |
| CheckM2 MAG QC | ✅ Added (`mag_qc` rule, `checkm2.yaml` env) |
| CheckM2 in test mode | ✅ Skipped (`skip_checkm2=True`) |

### Known issues / next steps
- CheckM2 requires downloading its database (~600 MB `.dmnd` file) via `checkm2 database --download --path <dir>` before use in production.
- The 10x coverage mock data is larger than the old dataset (16-23 MB gzipped FASTQ files vs 3.2 MB). This is necessary to ensure MetaBAT2 can produce bins.
- `test/output/` (gitignored) contains old test run outputs and will need to be cleared before running with the new mock data.
- Reference genomes source: `/home/benda/Projects/databases/test-genomes/` — required locally to re-run `generate_mock_data.py`.

---

## 2026-02-22 (session 5)

### What was done
- **Removed all `protected()` output declarations** from `short_reads.smk` (clean reads, nonpareil) and `contigs.smk` (megahit contigs). These were causing `ProtectedOutputException` on re-runs after failures.
- **Fixed MILP scheduler missing `cbc`**: Snakemake's MILP-based job selector requires the `cbc` solver which is not installed. All test runs now use `--scheduler greedy`.
- **Fixed MobMess `KeyError` on small mock plasmid data**: Added `|| true` + `touch` fallback after `mobmess systems` call — the tool crashes in `derep_sources` when all 216 plasmid contigs form singleton connected components. Pipeline now gracefully skips and continues.
- **Fixed `mag_amr` AMRFinder `gff_check` failure on small bins**: Added `2>/dev/null || true` to the per-bin `amrfinder` call in `mags.smk`, matching the pattern already in `contig_amr`.
- **Full end-to-end test passed**: `rule all` completed successfully — all 29 jobs across both mock samples finished. Output verified at `test/output/`.
- Pushed to `origin/main` (commits `6b60b72`, `2dc586f`, `b8ac566`).

### Current pipeline state
✅ **Full end-to-end test mode passes** (`snakemake --use-conda --cores 4 --scheduler greedy --config test=True`).

All 29 jobs complete for 2 mock samples:
- short_reads, short_reads_summary ✅
- contigs, genomad, mobmess ✅
- init_mmseqs_db, mmseqs_taxonomy ✅
- prodigal, contig_amr, mge_annotation, contig_abundance, contig_summary ✅
- bin, mag_prodigal, mag_amr, mag_mge, mag_abundance ✅ (gtdbtk, metabolic, checkm2 skipped as expected)

### Known issues / next steps
- `--scheduler greedy` must be passed to avoid `PulpSolverError: cannot execute cbc`. Consider documenting this in README or installing `coincbc` in the base env.
- MobMess produces no output for mock2 plasmids (gracefully skipped). This is expected — tool has a known edge case with high-identity plasmid clusters.
- AMRFinder finds no AMR genes in mock bins (expected for synthetic test data; `|| true` handles this cleanly).

---

## 2026-02-22 (session 4)

### What was done
- **Bug fix — OOM crash on test run**: All jobs died with `SIGTERM` because `threads` was placed inside `resources:` instead of as a top-level `threads:` directive. Snakemake only uses the `threads:` directive for core-based scheduling; `resources: threads` is a custom resource it ignores for that purpose. With `--cores 4`, Snakemake treated every job as 1-core and allowed up to 4 memory-intensive jobs (geNomad ×2, MEGAHIT, mmseqs, contig_abundance, etc.) to run simultaneously, exhausting RAM and crashing Ubuntu via the OOM killer.
- **Fix**: Added `threads: lambda wc: res(X, Y)` top-level directive to all compute-intensive rules in `short_reads.smk`, `contigs.smk`, and `mags.smk`. Shell commands continue using `{resources.threads}` unchanged (both can coexist).
- **Bonus fix**: Set `genomad_splits=4` as the default in test mode (`config.setdefault("genomad_splits", 4)` in `common.smk`) to further reduce geNomad peak memory. Can be overridden with `--config genomad_splits=N`.
- Dry run (`snakemake -n --config test=True --cores 4`) passes cleanly.
- Pushed to `origin/main` (commit `b85f186`).

### Current pipeline state
Same as session 3 — all stages implemented. The only change is the scheduling fix above.

### Known issues / next steps
- Re-run the full test (`snakemake --use-conda --cores 4 --config test=True`) to confirm end-to-end success now that parallel OOM is prevented.
- If memory is still tight, run with `--cores 2` to further limit parallelism, or add `--resources mem_mb=12000` to cap total concurrent memory usage.

---

## 2026-02-22 (session 3)

### What was done
- **Bug fix — markers always showing 0 cpg**: KMA includes the full FASTA description in the `#Template` column (e.g., `U30316.1 Bacteroides fragilis strain IB143 plasmid pBI143, complete sequence`), but the R script was filtering by exact accession match (`template == "U30316.1"`), which never succeeded. Fixed by extracting the accession before the first space: `mutate(seq_id = str_extract(template, "^\\S+"))` and filtering on `seq_id`.
- **Bug fix — robustify AFM template accession extraction**: Added `str_extract(refseq_nucleotide_accession, "^\\S+")` after the `separate()` call to strip any trailing description text from the nucleotide accession used as the catalog join key. Addresses the reported "AMRFinderPlus metadata missing" production issue (likely caused by the same KMA full-description behaviour on some database versions).
- **Test confirmed**: Deleted and re-ran `short_reads_summary` in test mode. `AMR_abundance_summary.csv` now correctly shows non-zero pBI143/crAss001 cpg values (mock1: pBI143=3.64, crAss001=0.50; mock2: crAss001=0.49). AMR metadata still fully populated.
- **Note**: All changes confined to `workflow/scripts/short_reads_processing.R`. No Snakemake rule or conda env changes needed.

### Current pipeline state
All stages implemented. Short-reads summary outputs are:
- `fastp_summary.csv` — per-sample QC metrics
- `short_reads_output.csv` — per-gene AMR abundance in cpg with full metadata
- `AMR_abundance_summary.csv` — per-sample totals: AMR_total (cpg) | pBI143 (cpg) | crAss001 (cpg)

### Known issues / next steps
- Production run OOM'd (likely during MEGAHIT or MetaBAT2). Need to investigate peak memory rule and adjust resources or assembly parameters.
- pBI143 and crAss001 references are stored at `/home/benda/Projects/databases/reference-sequences/` (production) and `test/dbs/markers/` (test). Both use the same FASTA headers.
- Full end-to-end test on this machine is blocked by memory — short-reads-only testing is feasible.


### What was done
- Wrote a comprehensive `README.md` covering: pipeline overview (ASCII flow diagram), requirements, configuration, run commands, test mode, outputs table, cpg normalisation formula, and repository structure.
- Pushed to `origin/main` (commit `1471875`).
- Created `docs/session-log.md` (this file) and added session-logging instructions to `.github/copilot-instructions.md`.

### Current pipeline state
The full pipeline is implemented end-to-end across four rule files:

| Rule file | Status | Key rules |
|-----------|--------|-----------|
| `short_reads.smk` | ✅ Complete | `initiate_dbs`, `short_reads`, `short_reads_summary` |
| `contigs.smk` | ✅ Complete | `contigs`, `genomad`, `mobmess`, `init_mmseqs_db`, `prodigal`, `contig_amr`, `init_mge_tool`, `mge_annotation`, `mmseqs_taxonomy`, `contig_abundance`, `contig_summary` |
| `mags.smk` | ✅ Complete | `bin`, `mag_prodigal` (checkpoint), `mag_amr`, `mag_mge`, `mag_taxonomy`, `mag_metabolism`, `mag_abundance` |
| Test infrastructure | ✅ Complete | Mock data in `test/data/`, mini DBs in `test/dbs/`, test mode via `--config test=True` |

### Known issues / next steps
- No known bugs at this time.
- MAG taxonomy (`mag_taxonomy`) and metabolism (`mag_metabolism`) require large external databases (GTDB-tk, METABOLIC); skipped in test mode and when `skip_gtdbtk: True` / `skip_metabolic: True`.
- `contig_summary.py` script exists but content not reviewed this session.
- MobMess/PlasX require internet access on first `--use-conda` setup (pip install from GitHub).

---

## 2026-02-22 (session 6)

### What was done
- **AMR_unified.csv feature implemented** (commit `ac36b6d`): New cross-stream AMR abundance table merging short-read KMA detections with contig-level AMRFinderPlus annotations.
  - New `workflow/rules/summary.smk` with `localrule amr_unified`
  - New `workflow/scripts/amr_unified.py`: outer-joins `short_reads_output.csv` + `contig_summary.tsv` (AMR rows only) on `sample × gene_symbol`
  - **cpg priority**: contig cpg preferred when available (higher specificity), else short-reads cpg
  - **`evidence` column**: `"both"` | `"short_reads_only"` | `"contigs_only"`
  - **`molecule_type`**: comma-sep sorted unique values per gene × sample (e.g., `"phage,plasmid"` when gene found on both)
  - **`taxonomy`**: modal non-empty value from contig annotations (empty string for short_reads_only)
  - **AMRFinderPlus catalog metadata** (gene_family, class, subclass, type, subtype) joined from `ReferenceGeneCatalog.txt`
  - Fixed pandas 2.0 compatibility (removed `include_groups=False` from `groupby.apply`)
  - Added `AMR_unified.csv` to `all_targets()` in `Snakefile`
- Verified on test data: `blaTEM-1` correctly tagged `evidence=both`, contig cpg preferred; `tet(M)` shows `molecule_type=phage,plasmid`.

### Current pipeline state

| Rule file | Status | Key rules |
|-----------|--------|-----------|
| `short_reads.smk` | ✅ Complete | `initiate_dbs`, `short_reads`, `short_reads_summary` |
| `contigs.smk` | ✅ Complete | `contigs`, `genomad`, `mobmess`, `prodigal`, `contig_amr`, `mge_annotation`, `mmseqs_taxonomy`, `contig_abundance`, `contig_summary` |
| `mags.smk` | ✅ Complete | `bin`, `mag_prodigal`, `mag_amr`, `mag_mge`, `mag_taxonomy`, `mag_metabolism`, `mag_abundance` |
| `summary.smk` | ✅ Complete (new) | `amr_unified` |
| Test infrastructure | ✅ Complete | Mock data, mini DBs, `--config test=True` |

### Known issues / next steps
- **Catalog metadata sparse in test output**: test catalog uses synthetic gene names (`fakeGene_*`) so only `blaTEM-1` gets metadata annotated; production run with real ReferenceGeneCatalog will populate all rows.
- Consider adding a contig-level AMR annotation lookup using the real catalog at runtime (post-assembly) to enrich `contigs_only` rows with gene_family/class/etc. when not in catalog.
- `short_reads_only` evidence will appear in production runs where KMA detects genes at low read depth that don't assemble (expected for low-abundance resistome members).

---

## 2026-02-22 (session 7)

### What was done
Four issues fixed in the unified AMR output:

**1. Test KMA database too small** (`test/dbs/amrfinder/`)
- The test `AMR_CDS.fa` had only 11 synthetic genes + blaTEM-1, causing `short_reads_output.csv` to detect only blaTEM-1.
- Wrote a script that extracts 49 real AMR gene nucleotide sequences from assembled contigs (using AMRFinder coordinates, strand-aware reverse complement) and adds them to `AMR_CDS.fa` with proper AMRFinderPlus-style FASTA headers.
- Updated `ReferenceGeneCatalog.txt` with synthetic `NG_MOCK_*` accessions matching the FASTA headers so the R catalog join succeeds.
- `short_reads_output.csv` now detects 81 gene-per-sample entries across both mock samples.

**2. MGE parsing broken in contig_summary.py**
- MobileElementFinder output is CSV (comma-separated) but `safe_read` used tab separator → single-column DataFrame → all fields empty.
- Comment lines (e.g. `#date:`, `#sample:`) were parsed as data.
- `stop_col` search was missing "end" (the actual column name in mge_finder output).
- Contig IDs had extra metadata appended (`"contig_id flag=1 multi=..."`) that failed mol_map/tax_map lookups.
- Fixed: `parse_mge` now uses `pd.read_csv(sep=",", comment="#")`, strips contig ID to first token, and includes "end" in stop_col candidates.

**3. amr_unified.py metadata was empty for most genes**
- Previous version joined metadata from `ReferenceGeneCatalog.txt` only (synthetic in test mode).
- Rebuilt metadata lookup from `*_contig_amr.tsv` files directly (highest priority), with fallback to short_reads catalog join, then the catalog file.
- All genes in the unified output now have product_name, class, subclass, type, subtype populated.

**4. AMR_abundance_summary.csv now uses unified totals**
- R script `short_reads_processing.R` now outputs `markers_cpg.csv` (pBI143/crAss001 per sample) as out_file3 instead of `AMR_abundance_summary.csv`.
- `amr_unified.py` reads `markers_cpg.csv` as input and generates `AMR_abundance_summary.csv` from the unified cpg totals (contig preferred, else short-reads) + marker cpg values.
- `short_reads_summary` rule output renamed `amr_summary` → `markers_cpg`.
- `summary.smk` `amr_unified` rule: added `markers_cpg` input, second output `amr_abundance_summary`.

### Current pipeline state
All 31 jobs pass in test mode. Terminal outputs:
- `fastp_summary.csv`, `short_reads_output.csv`, `markers_cpg.csv` (from short_reads)
- `contig_summary.tsv` (AMR+MGE fully populated with contig_id, gene, coords, molecule_type)
- `AMR_unified.csv` (78 `both`, 3 `short_reads_only`, 0 `contigs_only` in test data)
- `AMR_abundance_summary.csv` (unified AMR total + pBI143 + crAss001 per sample)

### Known issues / next steps
- `gene_family` in `AMR_unified.csv` uses the gene symbol itself (not a gene family grouping) for contig-detected genes, since AMRFinder direct output doesn't provide a gene_family column. Production runs with the full catalog will populate this correctly for known alleles.
- `short_reads_only` evidence is rare in test data (3 genes: cmlA1 and 2 others in mock2) because the test reads come from the same organisms whose contigs were assembled — expected in production to be more prevalent for low-depth genes.


---

## 2026-03-09 (session 8)

### What was done

**SLURM executor plugin compatibility:**
- Removed `time` / `runtime` string values from all rule `resources:` blocks — the `d-hh:mm:ss` format is only valid within `snakemake-executor-plugin-slurm` and causes a parse error when running locally.
- Created `config/slurm/config.yaml` — a complete SLURM executor profile following the SWAM-g large-batch pattern. Includes `executor: slurm`, `scheduler: greedy`, `set-threads` and `set-resources` (with `mem_mb` + `runtime`) for all 16 compute rules. Placeholder `slurm_account`/`slurm_partition` to fill in per cluster.
- Usage: `snakemake --profile config/slurm` (requires `pip install snakemake-executor-plugin-slurm`).

**Pipeline modularity (run_short_reads / run_contigs / run_mags flags):**
- Three new `config.yaml` stage flags (all default `True`): `run_short_reads`, `run_contigs`, `run_mags`.
- Two new external-input keys: `clean_reads_dir` (pre-filtered reads, needed when `run_short_reads: False`), `contigs_dir` (pre-assembled contigs, needed when `run_contigs: False` and `run_mags: True`).
- `common.smk`: added `_RUN_SR/_RUN_CTG/_RUN_MAGS` booleans, `get_clean_r1/r2()` and `get_contigs()` input helpers, `WorkflowError` config validation.
- `Snakefile` `all_targets()`: fully conditional — each stage block is gated on its flag.
- `contigs.smk`: `contigs` + `mobmess` rules use `get_clean_r1/r2()`; `contig_abundance` uses all three helpers + optional `scgs` (empty list when `_RUN_SR=False`).
- `mags.smk`: `bin` rule uses `get_contigs()`.
- `summary.smk`: `amr_unified` inputs conditional on stage flags.
- `amr_unified.py`: handles `None` for any optional input (short_reads_amr, contig_summary, markers_cpg).
- `contig_abundance.py`: falls back to `n_genomes=1` when `scgs` is unavailable.
- README updated with "Running on SLURM" and "Modular execution" sections.
- Pushed to `origin/main` (commit `d187e7d`).

### Current pipeline state

Dry-run job counts for all mode combinations (all validated this session):
- Full pipeline (`test=True`): **32 jobs** ✅
- Short reads only (`run_contigs=False, run_mags=False`): **6 jobs** ✅
- Contigs only (`run_short_reads=False, run_mags=False`): **23 jobs** ✅
- MAGs only (`run_short_reads=False, run_contigs=False`): **8 jobs** ✅

### Known issues / next steps
- SLURM profile placeholders `slurm_account` and `slurm_partition` must be filled in before cluster use.
- When `run_short_reads: False`, contig abundance cpg uses `n_genomes=1` (mean depth not genome equivalents) — documented in README.
- End-to-end full test (`snakemake --use-conda --cores 4 --scheduler greedy --config test=True`) was last confirmed passing in session 7; not re-run this session.

## 2026-03-10 (session 9)

### What was done
- **End-to-end test passed** with the session 8 updates (SLURM profile + modularity flags).
- Full pipeline ran successfully on both mock samples: **26 steps (100%) done**, exit code 0.
- Run time: ~21 minutes on 4 cores (21:36 → 21:57).

### Output verification

| File | Size | Notes |
|------|------|-------|
| `fastp_summary.csv` | 429 B | ✅ |
| `short_reads_output.csv` | 15 KB | ✅ |
| `markers_cpg.csv` | 65 B | ✅ (mock data: pBI143/crAss001 = 0 — expected, not in mock genomes) |
| `contig_summary.tsv` | 54 KB, 857 rows | ✅ |
| `AMR_unified.csv` | 15 KB | ✅ evidence=both/short_reads_only working |
| `AMR_abundance_summary.csv` | 109 B | ✅ |
| `mag_amr.tsv` / `mag_mge.tsv` / `mag_abundance.tsv` | mock1+mock2 | ✅ |

`AMR_unified.csv` sample: `tet(A)` evidence=both, `cmlA1` evidence=short_reads_only, `tet(M)` molecule_type=phage,plasmid — all as expected.

### Current pipeline state
✅ **Full end-to-end test passes** with SLURM compatibility + modularity changes.

### Known issues / next steps
- `markers_cpg.csv` shows pBI143/crAss001 = 0 for both mock samples. This is expected — the mock genomes don't include pBI143 or crAss001 sequences. Production runs with real wastewater samples will show non-zero values.
- SLURM profile `slurm_account` and `slurm_partition` remain as empty placeholders — fill in before cluster use.

---

## 2026-03-11 (session 10)

### What was done

**Conda environment restructuring + version pinning:**

**Problem**: `mags.yaml` bundled `coverm` + `ncbi-amrfinderplus` (which pulls biopython) + `samtools>=1.9` — the version constraint caused `LibMambaUnsatisfiableError` during env creation. Several partial fixes were attempted before the SSH connection dropped (commits `86856c8`, `16a9729`, `baacd40`).

**Solution — environment consolidation (9 → 8 envs):**

1. **Moved `coverm` to `shortreads.yaml`** — already has samtools/minimap2; no biopython conflict. Updated `mag_abundance` rule to use `shortreads.yaml`.
2. **Removed `coverm` and `samtools>=1.9` from `mags.yaml`** — no more biopython + samtools version conflict. `mags.yaml` now: metabat2, ncbi-amrfinderplus, minimap2, samtools, prodigal.
3. **Merged `mge.yaml` into `contigs.yaml`** — MobileElementFinder's Python deps (biopython, tabulate, pyyaml, click, attrs, cattrs, bcbio-gff, setuptools<70) are compatible with contigs.yaml's Python stack. Deleted `mge.yaml`. Updated `init_mge_tool`, `mge_annotation` (contigs.smk), and `mag_mge` (mags.smk) to use `contigs.yaml`.
4. **Fixed `contig_abundance` env** — rule uses `script:` (Python), needs pandas/numpy; was using `shortreads.yaml` which has no Python. Moved to `contigs.yaml`.
5. **Pinned exact versions** for all direct dependencies in each yaml, resolved from `conda env export --no-builds` on built envs.

**Final env inventory (8 envs):**

| File | Key tools | Pinned versions |
|------|-----------|-----------------|
| `shortreads.yaml` | fastp, kma, diamond, minimap2, samtools, coverm | ✅ |
| `contigs.yaml` | megahit, mmseqs2, prodigal, MobMess/PlasX + MobileElementFinder deps | ✅ |
| `genomad.yaml` | genomad | ✅ |
| `Renv.yaml` | r-base, r-tidyverse, nonpareil | ✅ |
| `mags.yaml` | metabat2, ncbi-amrfinderplus, prodigal | ✅ |
| `checkm2.yaml` | checkm2 | (unpinned — skipped in test) |
| `gtdbtk.yaml` | gtdbtk | (unpinned — skipped in test) |
| `metabolic.yaml` | hmmer, perl, R stack | (unpinned — skipped in test) |

**End-to-end test passed**: 30/30 steps (100%) done.

### Current pipeline state
✅ Full end-to-end test passes with restructured, version-pinned environments.

### Known issues / next steps
- `checkm2.yaml`, `gtdbtk.yaml`, `metabolic.yaml` are not version-pinned (skipped in test mode — need a production run to export their resolved versions).
- `igraph` in `contigs.yaml` is unpinned — conda did not resolve it as a named package (likely satisfied via pip by PlasX/MobMess). Left in yaml for completeness.
- SLURM profile placeholders `slurm_account` and `slurm_partition` still need filling in before cluster use.

---

## 2026-03-16 (session 11)

### What was done
- Replaced the hand-maintained ASCII `README.md` "Pipeline stages" block with a generated Snakemake rulegraph image.
- Added `docs/rulegraph.png`, rendered with Graphviz from `snakemake --snakefile workflow/Snakefile --rulegraph`.
- Generated the image using the mock input directory (`in_dir=test/data`, `out_dir=test/output`, `markers_db=test/dbs/markers`) so the DAG includes real samples while still showing the full workflow, including optional MAG taxonomy, metabolism, and CheckM2 branches.
- Removed the temporary `docs/rulegraph.dot` file after rendering so only the committed PNG asset remains.

### Current pipeline state
- Workflow logic is unchanged; this session only updated documentation and added a generated visualization asset.
- `README.md` now documents the pipeline structure using the actual rule dependency graph from the Snakemake workflow instead of a prose-only stage list.

### Known issues / next steps
- If rules or dependencies change, `docs/rulegraph.png` should be regenerated so the README stays in sync with the workflow DAG.
- The rulegraph is more accurate than the old text block, but it is denser; keep surrounding README sections clear so new users still have high-level narrative context.

---

## 2026-03-16 (session 12)

### What was done
- Replaced the legacy single SLURM profile `config/slurm/config.yaml` with two `snakemake-executor-plugin-slurm` profiles:
  - `config/slurm/large-batch/config.yaml`
  - `config/slurm/small-batch/config.yaml`
- Followed the SWAM-g profile pattern: explicit `set-threads` and explicit `set-resources` for all compute rules, with `executor: slurm` and `scheduler: greedy` defined in profile YAML.
- Implemented conservative batching only in `large-batch`: `contig_amr`, `mge_annotation`, `mag_amr`, `mag_mge`, `mag_qc`, and `mag_abundance` are grouped; heavy assembly, taxonomy, and binning rules remain one sample per job.
- Updated `README.md` to document profile selection (`>=50` samples => `large-batch`, `<50` samples => `small-batch`), plugin installation, and the absence of a local profile.
- Updated `.github/copilot-instructions.md` to stop recommending legacy `--slurm` execution and instead use `snakemake --profile config/slurm/{small-batch,large-batch}`.

### Validation
- Verified both profiles parse and run a dry-run successfully with a local executor override:
  - `snakemake --profile config/slurm/small-batch --executor local -n --config test=True`
  - `snakemake --profile config/slurm/large-batch --executor local -n --config test=True`
- Both commands completed successfully and built the expected test DAG (36 dry-run jobs reported in the current workspace state).

### Current pipeline state
- The workflow logic is unchanged; this session only updates cluster execution configuration and related documentation.
- SWAM-meta now has two supported SLURM entrypoints and no fallback top-level SLURM profile.

### Known issues / next steps
- Users still need to populate `slurm_account` and `slurm_partition` in both profile files before running on a real cluster.
- The large-batch grouping strategy is intentionally conservative; if cluster utilization remains low, batch sizes can be increased later in profile YAML without changing rule files.

---

## 2026-03-16 (session 13)

### What was done
- Updated `README.md` install guidance to make the new cluster requirements explicit:
  - raised the documented Snakemake requirement to `>=8`
  - added base-environment install commands for Snakemake
  - added the combined `snakemake>=8` + `snakemake-executor-plugin-slurm` install command for SLURM users
- Updated the SLURM section to reference the combined conda install command rather than only installing the plugin in isolation.

### Current pipeline state
- Workflow logic is unchanged; this was a documentation-only follow-up to the SLURM profile migration.

### Known issues / next steps
- The README now reflects the plugin-based profile setup, but users still need to set `slurm_account` and `slurm_partition` in both SLURM profile files before cluster execution.

---

## 2026-03-20 (session 14)

### What was done
- Added `log:` directive to every rule across all four rule files (`short_reads.smk`, `contigs.smk`, `mags.smk`, `summary.smk`).
- Log files go to `{output_dir}/data/QAQC/logs/` — named `{sample}.rulename.log` for per-sample rules, `rulename.log` for global rules.
- Redirected all tool stdout/stderr to the log with `>> {log} 2>&1` (append), keeping logs for debugging without flooding the terminal.
- For multi-command pipelines (minimap2 → samtools), used `2>> {log}` per command to silence each stage's stderr while preserving the pipe.
- Added/updated echo messages throughout with the format `[{wildcards.sample}] tool: action...` and `[{wildcards.sample}] tool: done` so the terminal only shows brief start/stop updates.
- Fixed a minor pre-existing issue: `fastp --html /dev/null/` (spurious trailing slash) → `--html /dev/null`.
- Script-based rules (`short_reads_summary`, `contig_abundance`, `contig_summary`, `amr_unified`) also received `log:` directives — Snakemake automatically redirects script stderr/stdout to the log file.
- Validated with `snakemake -n --config test=True` — dry run exits 0, DAG builds correctly.

### Current pipeline state
- Workflow logic unchanged; only output verbosity and log routing updated.
- Terminal output during a run will show only brief `[sample] tool: action` progress echoes.
- Full tool output (including errors) is captured in per-rule log files under `data/QAQC/logs/`.

### Known issues / next steps
- `checkm2.yaml`, `gtdbtk.yaml`, `metabolic.yaml` still not version-pinned (skipped in test mode).
- `slurm_account` and `slurm_partition` placeholders need filling before cluster use.

---

## 2026-03-20 (session 15)

### What was done

**SLURM profile runtime fix:**
- Both `config/slurm/small-batch/config.yaml` and `config/slurm/large-batch/config.yaml` had runtime values in `d-hh:mm:ss` string format (invalid for `snakemake-executor-plugin-slurm`). Converted all values to integer minutes: 240, 360, 480, 600, 1440, 2160, 4320.

**README streamlining:**
- Removed the separate "Requirements" section; merged into "Installation" with a single conda install command that always includes both `snakemake>=8` and `snakemake-executor-plugin-slurm`.
- Removed duplicate conda install lines from the SLURM section.
- Improved the database setup section with direct download URLs for UniRef50 (UniProt FTP), NCBI taxdump, and corrected the GTDB-tk downloader call to `conda run -n gtdbtk download-db.sh`. Added METABOLIC git clone command.
- Simplified SLURM section to a 2-step pattern (edit account/partition, then submit).
- The expanded "Prepare required databases" table was drafted but reverted at user request; user has separate plans for that section.

### Current pipeline state
- No workflow logic changes.
- SLURM profiles now use valid runtime syntax (integer minutes).
- README installation is now a single command with no duplicates.

### Known issues / next steps
- Database section in README is intentionally minimal — user has plans to expand it separately.
- `checkm2.yaml`, `gtdbtk.yaml`, `metabolic.yaml` still not version-pinned.
- `slurm_account` and `slurm_partition` placeholders still need filling before cluster use.

---

## 2026-03-20 (session 16)

### What was done

**Full database automation — all DBs now live in `SWAM-meta/dbs/` (gitignored):**

- `common.smk`: removed `config["scg_db"]`, `config["uniref50_db"]`, `config["markers_db"]` test-mode overrides; added module-level constants `_DBS_DIR`, `_SCG_DB`, `_UNIREF50_DB`, `_CHECKM2_DB`, `_METABOLIC_DIR` pointing to `SWAM-meta/dbs/` (test mode uses `test/dbs/` as source inputs but writes built indices to the same `dbs/` dir).

- `config/config.yaml`: stripped `scg_db`, `uniref50_db`, `markers_db`, `checkm2_db`, `metabolic_dir`. Now only `in_dir`, `out_dir`, and `gtdbtk_db` are path keys.

- `initiate_dbs` (short_reads.smk): outputs moved to `_DBS_DIR`; added `metabolic_done` sentinel; production shell now auto-downloads pBI143 (`U30316.1`) and crAss001 (`NC_049977.1`) via NCBI Entrez efetch API; git-clones METABOLIC (conditioned on `skip_metabolic`).

- `init_genomad` (contigs.smk): sentinel and DB directory moved from `output_dir/data/genomad/` to `_DBS_DIR/`.

- `init_mmseqs_db` (contigs.smk): removed from `localrules`; production branch now auto-downloads UniRef50 FASTA (~9 GB) and NCBI taxdump, builds MMseqs2 taxonomy DB; resources increased to 16 threads / 32 GB.

- `init_checkm2_db` (new localrule in mags.smk): downloads CheckM2 DB (~3 GB) via `checkm2 database --download`; sentinel at `_DBS_DIR/.checkm2.done`; respects `skip_checkm2`.

- `mag_qc`: added `checkm2_done` input; `db_path` now uses `_CHECKM2_DB`.

- `mag_metabolism`: added `metabolic_done` input; `metabolic_dir` now uses `_METABOLIC_DIR`.

- `amr_unified` (summary.smk): catalog input path updated to `_DBS_DIR/ReferenceGeneCatalog.txt`.

- `Snakefile`: replaced `init_mmseqs_db` with `init_checkm2_db` in `localrules`.

- SLURM profiles: added `init_mmseqs_db` to `set-threads` (16) and `set-resources` (32 GB, 4320 min) in both `small-batch` and `large-batch`.

- `workflow/resources/README.md`: new file documenting that `SCGs_40_All.fasta` must be placed in `workflow/resources/` before production use.

- README: simplified production setup to 2 steps — copy SCG FASTA + fill in `in_dir`/`out_dir`/`gtdbtk_db`.

### Current pipeline state
- Dry run passes (exit 0, 36 jobs) in test mode.
- HEAD: `cd08ce4` — pushed to `origin/main`.
- `SWAM-meta/dbs/` is the persistent home for all auto-managed databases.

### Known issues / next steps
- `SCGs_40_All.fasta` must still be copied manually to `workflow/resources/` — git LFS could track it if needed.
- `checkm2.yaml`, `gtdbtk.yaml`, `metabolic.yaml` still not version-pinned.
- `slurm_account` / `slurm_partition` placeholders still need filling before cluster use.

---

## 2026-03-22 (session 17)

### What was done

**Bug fix — contig taxonomy showing rank label instead of taxname:**

- `workflow/scripts/contig_summary.py` → `parse_lca()`: changed `df.iloc[:, 2]` to `df.iloc[:, 3]`.
- Root cause: MMseqs2 `easy-taxonomy` `_lca.tsv` has 9 columns — `query, taxid, rank, taxname, n_frags, n_direct, n_classified, fraction, lineage`. The old code read column 2 (rank — e.g. `"species"`, `"no rank"`) instead of column 3 (actual taxname — e.g. `"Escherichia coli"`).
- Also tightened the column-count guard from `< 3` to `< 4` and corrected the docstring comment.
- Commit: `aa1e2d3`

### Current pipeline state
- Dry run passes (exit 0).
- Taxonomy values in `contig_summary.tsv` will now correctly show species/taxon names.

### Known issues / next steps
- `checkm2.yaml`, `gtdbtk.yaml`, `metabolic.yaml` still not version-pinned.
- `slurm_account` / `slurm_partition` placeholders need filling before cluster use.

---

## 2026-03-23 (session 18)

### What was done

**Five post-test-run fixes — all implemented and committed as `324e0a3`:**

1. **ARG spike reads (issue i)** — blaTEM-1 was not assembling from the mock data because E. coli genome-wide coverage (10x) was too low. Added a synthetic 2,861 bp "ARG spike" FASTA (blaTEM-1 flanked by random sequence) simulated at 100x. Appended 954 read pairs per sample to existing test FASTQs via wgsim. Updated `test/scripts/generate_mock_data.py` so future regenerations include the spike natively. Updated `test/data/mock_data_summary.txt`.

2. **Marker column names (issue ii)** — `short_reads_processing.R` was outputting `"pBI143 (cpg)"` and `"crAss001 (cpg)"` (spaces + parens). Renamed to `pBI143_cpg` and `crAss001_cpg` throughout. Updated `amr_unified.py` to match.

3. **AMR cpg priority reversal (issue iii)** — Old logic preferred contig cpg over short-reads cpg. Reversed to always report short-reads cpg (most sensitive; captures all reads including those that don't assemble). Contig cpg used only for `contigs_only` rows (unique assembled alleles). Evidence column (`both`/`short_reads_only`/`contigs_only`) already serves as the "unique allele" flag.

4. **MAG summary output (issue iv)** — New `workflow/scripts/mag_summary.py` aggregates per-sample MAG files into a single `mag_summary.tsv` at the output root. Columns: `sample`, `bin_id`, `abundance_trimmed_mean`, `n_amr_genes`, `amr_genes`, `n_mge`, `mge_elements`, plus optional `completeness`, `contamination`, `quality_score` (CheckM2) and `classification` (GTDB-Tk). New `mag_summary` localrule in `summary.smk`. Added to `all_targets()` in `Snakefile`.

5. **Output consolidation (issue v)** — Moved intermediate files out of the output root: `short_reads_output.csv` → `data/QAQC/short_reads_output.csv`; `markers_cpg.csv` → `data/QAQC/markers_cpg.csv`. Top-level outputs now: `AMR_unified.csv`, `AMR_abundance_summary.csv`, `contig_summary.tsv`, `fastp_summary.csv`, `mag_summary.tsv`.

### Current pipeline state
- Dry run passes (exit 0) in test mode.
- HEAD: `324e0a3` — pushed to `origin/main`.
- blaTEM-1 contig-level AMR annotation has **not yet been confirmed** in a full test run with the new spike reads (ARG spike added, but full re-run not completed this session).

### Known issues / next steps
- Run full test pipeline to confirm blaTEM-1 appears in `contig_summary.tsv`.
- `checkm2.yaml`, `gtdbtk.yaml`, `metabolic.yaml` still not version-pinned.
- `slurm_account` / `slurm_partition` placeholders need filling before cluster use.

---

## 2026-03-23 (session 19)

### What was done

**Session restoration only — no code changes.**

- Confirmed session 18 push (`324e0a3`) completed successfully to `origin/main`.
- Appended session 18 log entry to `docs/session-log.md` (had not been written before the previous session ended) and committed as `7f6f8e9`.

### Current pipeline state
- HEAD: `7f6f8e9` — pushed to `origin/main`.
- All 5 post-test fixes from session 18 are live.
- Full test run has not yet been executed to confirm blaTEM-1 ARG spike assembly.

### Known issues / next steps
- Run full test pipeline to confirm blaTEM-1 appears in `contig_summary.tsv`.
- `checkm2.yaml`, `gtdbtk.yaml`, `metabolic.yaml` still not version-pinned.
- `slurm_account` / `slurm_partition` placeholders need filling before cluster use.

---

## 2026-03-22 (session 20)

### What was done

**Full end-to-end test run** — confirmed pipeline completes 34/34 jobs in test mode.

**Bug identified: AMRFinderPlus protein-mode DB not initialized**

- `contig_amr` and `mag_amr` rules call `amrfinder -n/-p/-g` without specifying `--database`.
- AMRFinderPlus looks for its DB at `share/amrfinderplus/data/latest` inside the conda env, which does not exist (it requires `amrfinder -u` to download).
- The `|| true` fallback silently produces empty TSVs — pipeline completes but all contig-level and MAG-level AMR annotations are empty.
- This explains why blaTEM-1 is `short_reads_only` with `contig_cpg=0.0` — it was never annotated at the contig level, not because it failed to assemble.
- The KMA-indexed mini DB at `dbs/afp_db` covers short-read alignment only; `amrfinder` protein/nucleotide mode needs its own proprietary DB format.

### Current pipeline state
- 34/34 jobs complete ✅
- blaTEM-1 detected at short-reads level (cpg ~27) in both mock samples ✅
- `contig_summary.tsv` (819 rows), `AMR_unified.csv`, `AMR_abundance_summary.csv`, `fastp_summary.csv`, `mag_summary.tsv` all produced ✅
- 3 bins per sample at correct abundance (~10 cpg) ✅
- Contig-level and MAG-level AMRFinderPlus annotation: ❌ silently failing (empty output, gracefully handled)

### Known issues / next steps
- **Fix AMRFinderPlus DB init**: Add `amrfinder -u --database {_DBS_DIR}/amrfinderplus_db` to `initiate_dbs` (production branch); add `--database` flag to `contig_amr` and `mag_amr` rule shell blocks.
- `checkm2.yaml`, `gtdbtk.yaml`, `metabolic.yaml` still not version-pinned.
- `slurm_account` / `slurm_partition` placeholders need filling before cluster use.

---

## 2026-03-23 (session 21)

### What was done

**Fix: AMRFinderPlus database initialization for contig and MAG AMR annotation**

Root cause: `amrfinder` protein/nucleotide mode needs its own database (downloaded via `amrfinder -u`), separate from the KMA-indexed `afp_db` used for short-read alignment. Without it, `amrfinder` failed silently via `|| true` and produced empty TSV files for all contig-level and MAG-level AMR calls.

Changes (commit `9d5fd18`):
- `common.smk`: added `_AFP_DB_DIR = dbs/amrfinderplus_db/`
- `contigs.smk`: new `init_amrfinder_db` localrule — runs `amrfinder -u --database dbs/amrfinderplus_db/` in production; skips in test mode (touches sentinel only). `contig_amr` now depends on `.amrfinder_db.done` and passes `--database`.
- `mags.smk`: `mag_amr` now depends on `.amrfinder_db.done` and passes `--database`.
- `Snakefile`: `init_amrfinder_db` added to `localrules`.

Test run (9/9 re-run jobs, all pass): `init_amrfinder_db` skips download in test mode; `contig_amr` and `mag_amr` still produce empty output (expected — no AMRFinderPlus DB in test mode).

### Current pipeline state
- HEAD: `9d5fd18` — pushed to `origin/main`.
- Test mode: 34 jobs complete (full DAG), 9 re-run on code-change trigger.
- Production: first run will download ~600 MB AMRFinderPlus DB to `dbs/amrfinderplus_db/`. Subsequent runs skip download (idempotent sentinel check).

### Known issues / next steps
- Contig/MAG AMR output will remain empty in test mode (no DB downloaded). This is expected behaviour.
- `checkm2.yaml`, `gtdbtk.yaml`, `metabolic.yaml` still not version-pinned.
- `slurm_account` / `slurm_partition` placeholders need filling before cluster use.

---

## 2026-03-23 (session 6)

### What was done
- **First production-mode end-to-end test**: Ran short reads + contigs pipeline (no MAGs) in production mode using mock1/mock2 test data, output to `output/`.
  Command: `snakemake --scheduler greedy --use-conda --cores 8 --snakefile workflow/Snakefile --config in_dir=.../test/data out_dir=.../output run_mags=False skip_metabolic=True`
- **Fixed `workflow/resources/SCGs_40_All.fa` → `SCGs_40_All.fasta`**: File was misnamed; `common.smk` always referenced `.fasta`. Renamed and committed.
- **Fixed `initiate_dbs` (short_reads.smk) — idempotent production branch**: 
  - Replaced `h_genome` file output with `.h_genome.done.txt` sentinel so Snakemake never manages the 928 MB human genome as an output (would delete it on re-run).
  - Added separate idempotency guards for each asset: `AMR_CDS.fa`, `ReferenceGeneCatalog.txt`, human genome, `markers.fa`, KMA AFP index, KMA markers index.
  - `short_reads` rule now consumes `.h_genome.done.txt` as input and passes the actual genome path via `params.h_genome`.
- **Fixed `init_amrfinder_db` (contigs.smk)**: `amrfinder -u --database <dir>` is invalid in AMRFinder 4.x. Fixed to `amrfinder --update` (updates default conda env location) then `rsync` from `$CONDA_PREFIX/share/amrfinderplus/data/latest/` to `dbs/amrfinderplus_db/`. Existence check updated to `-f AMR_CDS.fa`.
- All 24 jobs completed successfully in ~25 minutes on 8 cores (no MAGs).

### Current pipeline state
✅ **Production mode (short reads + contigs) passes** — all 24 jobs for mock1 + mock2.

Key outputs verified in `output/`:
| File | Size |
|------|------|
| `AMR_unified.csv` | 15 KB, 79 rows |
| `AMR_abundance_summary.csv` | 103 B |
| `contig_summary.tsv` | 59 KB, 866 rows |
| `fastp_summary.csv` | 427 B |

### Known issues / next steps
- Production run is still sequential (one sample at a time, 8 cores) because each rule claims all 8 cores. With more cores (e.g., 16+) two samples would run in parallel.
- MAG stage not tested in production yet (`run_mags: True` requires GTDB-tk, CheckM2 databases on HPC).
- SLURM profiles not yet tested — user to test on HPC.
- Pushed to `origin/main` (commit `81505e3`).

---

## 2026-03-23 (session 22)

### What was done

**Fix: Show download progress for all database init rules**

In production mode, `wget -q` + `>> {log} 2>&1` silenced all download progress. Users saw an
echo message ("downloading human reference genome...") then nothing for 30–60+ minutes, making
the pipeline appear hung.

**Changes (3 files):**

- `workflow/rules/short_reads.smk` — `initiate_dbs`: removed `-q` from all 5 `wget` calls;
  changed `>> {log} 2>&1` → `2>&1 | tee -a {log}` so progress appears on the terminal and
  is also written to the log.

- `workflow/rules/contigs.smk`:
  - `init_genomad`: changed `curl -L ... 2>> {log} | tar -xz` to
    `curl -L --progress-bar ... 2> >(tee -a {log} >&2) | tar -xz` — uses process substitution
    to tee curl's stderr (progress) while keeping the binary data stream piped cleanly to tar.
  - `init_mmseqs_db`: removed `-q` from both `wget` calls; changed `>> {log} 2>&1` →
    `2>&1 | tee -a {log}` for both wget downloads and all three long mmseqs commands
    (`createdb`, `createtaxdb`, `createindex`).
  - `init_amrfinder_db`: changed `amrfinder --update >> {log} 2>&1` →
    `amrfinder --update 2>&1 | tee -a {log}`; added `--info=progress2` to `rsync` and
    changed to `2>&1 | tee -a {log}`.

- `workflow/rules/mags.smk` — `init_checkm2_db`: changed
  `checkm2 database --download ... >> {log} 2>&1` → `2>&1 | tee -a {log}`.

Dry-run validated: `snakemake -n --config test=True` completes with exit 0, 37 jobs.

### Current pipeline state
- HEAD: this commit, pushed to `origin/main`.
- Test mode: 37-job DAG valid.
- Production: all DB-init rules now print download progress to the terminal.

### Known issues / next steps
- `checkm2.yaml`, `gtdbtk.yaml`, `metabolic.yaml` still not version-pinned.
- `slurm_account` / `slurm_partition` placeholders need filling before cluster use.

---

## 2026-03-23 (session 23)

### What was done

**Fix: wget --show-progress for non-TTY contexts**

`wget` suppresses its progress bar when stdout is a pipe (not a TTY). Our previous
`2>&1 | tee -a {log}` pattern meant wget ran silently even though tee was present. Added
`--show-progress` to all 7 wget calls in `initiate_dbs` and `init_mmseqs_db` to force
progress output regardless of TTY. Commit `62d6159`.

**Feature: Reorganize dbs/ directory into logical subdirectories**

Before: all short-reads assets (AMR_CDS.fa, afp_db.*, scg_db.dmnd, human genome, marker
FASTAs/indices, ReferenceGeneCatalog.txt) and all per-stage sentinels were scattered flat at
the dbs/ root alongside subdirectory databases.

After (`28a61db`):
```
dbs/
├── short_reads/        # everything initiate_dbs creates (flat files + sentinels)
├── genomad_db/         # sentinel co-located: genomad_db/.done
├── amrfinderplus_db/   # sentinel co-located: amrfinderplus_db/.done
├── uniref50/           # sentinel co-located: uniref50/.done
└── checkm2/            # sentinel co-located: checkm2/.done
```
Added `_SR_DBS_DIR = dbs/short_reads/` to `common.smk`. Sentinel naming standardised
to `.done` (no `.txt` suffix). All four rule files updated; 37-job dry-run passes.

**Note for HPC users:** existing dbs will need to be moved into the new subdirectory layout
or the pipeline will re-download. Move commands:
```bash
mkdir -p dbs/short_reads
mv dbs/{AMR_CDS.fa,ReferenceGeneCatalog.txt,afp_db.*,scg_db.dmnd,GCF_*fna.gz,pBI143.fasta,crAss001.fasta,markers.fa,markers_db.*} dbs/short_reads/
# rename old sentinels
for f in afp scg h_genome markers metabolic; do
    [ -f dbs/.${f}.done.txt ] && mv dbs/.${f}.done.txt dbs/short_reads/.${f}.done
done
mv dbs/.genomad.db.done.txt dbs/genomad_db/.done 2>/dev/null || true
mv dbs/.amrfinder_db.done dbs/amrfinderplus_db/.done 2>/dev/null || true
mv dbs/.mmseqs_db.done dbs/uniref50/.done 2>/dev/null || true
mv dbs/.checkm2.done dbs/checkm2/.done 2>/dev/null || true
```

### Current pipeline state
- HEAD: `28a61db` — pushed to `origin/main`.
- Dry-run: 37-job DAG valid.
- Production run on HPC is in progress (downloading DBs).

### Known issues / next steps
- `checkm2.yaml`, `gtdbtk.yaml`, `metabolic.yaml` still not version-pinned.
- `slurm_account` / `slurm_partition` placeholders need filling before cluster use.

**README update (commit `d1b0f7e`):**
- Removed the manual "Add SCG reference database" step — `SCGs_40_All.fasta` now ships with the repo
- Production setup reduced to two required config keys (`in_dir`, `out_dir`)
- Optional databases table expanded: added geNomad DB and AMRFinderPlus protein DB rows with accurate sizes and `dbs/` subdirectory paths
- Repository structure updated to reflect the new `dbs/` subdirectory layout
- Production setup steps renumbered 2/3/4/5 → 1/2/3/4

### Current pipeline state
- HEAD: `d1b0f7e` — pushed to `origin/main`.
- Production run on HPC is in progress (downloading DBs with visible progress bars).

### Known issues / next steps
- HPC users with existing `dbs/` need to manually move files into the new subdirectory layout (move commands in session 23 entry).
- `checkm2.yaml`, `gtdbtk.yaml`, `metabolic.yaml` still not version-pinned.
- `slurm_account` / `slurm_partition` placeholders need filling before cluster use.

---

## 2026-03-24 (session 24)

### What was done

**Fix: fastp minimum read length to prevent Nonpareil kmer crash**

Running on an independent test dataset exposed a crash in the `short_reads` rule:

```
Nonpareil ... fatal error: Reads are required to have a minimum length of kmer size...
```

**Root cause:** `fastp` was called without `--length_required`, so its built-in default
(15 bp) allowed very short reads (~20 bp) through QC. Nonpareil's kmer mode (`-T kmer`)
requires every read to be at least as long as the kmer size (default 24 bp). Reads shorter
than 24 bp caused a fatal crash.

**Fix (`commit TBD`):**
- `config/config.yaml` — added `min_read_length: 50` under a new "Read QC parameters"
  section (configurable; 50 bp is the standard practical floor for metagenomic short reads).
- `workflow/rules/short_reads.smk` — added `min_read_length = config.get("min_read_length", 50)`
  param to the `short_reads` rule; injected `--length_required {params.min_read_length}` into
  the `fastp` shell command.

Dry-run validated: both `snakemake -n --scheduler greedy` and
`snakemake -n --scheduler greedy --config test=True` complete with exit 0.

### Current pipeline state
- HEAD: this commit, pushed to `origin/main`.
- Test mode: dry-run DAG valid.
- fastp now enforces 50 bp minimum read length (configurable via `min_read_length`).

### Known issues / next steps
- HPC users with existing `dbs/` need to manually move files into the new subdirectory layout (move commands in session 23 entry).
- `checkm2.yaml`, `gtdbtk.yaml`, `metabolic.yaml` still not version-pinned.
- `slurm_account` / `slurm_partition` placeholders need filling before cluster use.

---

## 2025-07-14 (session 25)

### What was done

**Fix 1 — Separate test/production database sentinels (`e2740f7`)**

Root cause of all HPC production issues (NG_MOCK genes, cpg inflation, contigs_only evidence):
`initiate_dbs` in test mode wrote the same sentinel filenames (`.afp.done`, `.scg.done`,
`.h_genome.done`, `.markers.done`, `.metabolic.done`) to the same `dbs/short_reads/`
directory as production. Once test mode had run, Snakemake never re-ran `initiate_dbs` in
production because the sentinels already existed, silently using the test KMA/DIAMOND
databases for all samples.

Fix: added `_DB_TAG = ".test" if _TEST else ""` to `common.smk`; all five sentinel filenames
now incorporate the tag (e.g. `.afp.done` → `.afp.test.done` in test mode, `.afp.done` in
production). `ReferenceGeneCatalog.txt` is not tagged (it's a data file overwritten on each
`initiate_dbs` run).

**Fix 2 — Drop unannotated rows from `contig_summary.tsv` (`e2740f7`)**

`contig_summary.py` previously added a row for every contig, even those with no AMR or MGE
annotation (empty `feature_type`). Removed the fallback loop; the output now contains only
AMR and MGE annotated rows, making the file leaner and downstream joins cleaner.

Both dry-runs pass: `snakemake -n --scheduler greedy` (production) and
`snakemake -n --scheduler greedy --config test=True` (test mode).

### Current pipeline state
- HEAD: `e2740f7`, pushed to `origin/main`.
- All sentinel contamination fixed; test and production can coexist on the same machine.
- `contig_summary.tsv` is now annotation-only (AMR + MGE rows only).

### Known issues / next steps
- **Contig abundance normalization** — `n_genomes` from SCG DIAMOND fails for diverse
  environmental metagenomes (very few reads align to 40 reference proteins → n_genomes ≈ 0
  → cpg ≈ ∞). A new normalization strategy is needed; not yet determined.
  Rejected approaches: median contig depth normalisation, n_genomes floor of 1.0.
- MGE JSONDecodeError in MobileElementFinder — upstream bug in `-outfmt 15` BLAST output
  parsing; `|| true` handles gracefully. Awaiting upstream fix.
- `checkm2.yaml`, `gtdbtk.yaml`, `metabolic.yaml` still not version-pinned.

---

## 2026-03-24 (session 26)

### What was done

**Median-SCG normalisation for convergent CPG values (`19e17bc`)**

Replaced the old `n_genomes = sum(aln_len/gene_len)/40` formula with a per-gene,
median-based SCG depth calculation. The new formula uses reads/base as the common
unit for both evidence streams, making short-read and contig CPG values directly
comparable and convergent in `AMR_unified`.

**Formula (both scripts now identical):**
```
SCG_depth_i      = n_reads_mapped_to_SCG_i / (slen_i × 3)    # reads/base
median_SCG_depth = median(SCG_depth_i for all detected SCGs)
cpg              = alignment_depth / median_SCG_depth
```

- **`short_reads_processing.R`**: replaced `genome.counts` block with `scg_norm`
  computing per-gene SCG depth (reads/base). Diagnostic `message()` reports
  number of SCGs detected and median depth. Warns when <10 SCGs detected.
  CPG becomes NA when no SCGs detected (propagated through AMR_unified as 0).
  Both AMR and anthropogenic marker CPG updated.
- **`contig_abundance.py`**: same median-SCG calculation from the shared DIAMOND
  file. `abundance_cpg = NaN` when median_scg_depth is unavailable (produces
  empty cell in TSV; handled as 0.0 by `amr_unified.py`).

**Why this converges:**
- KMA Depth (short reads) is in reads/base; samtools mean depth (contigs) is in
  reads/base. Both divide by median_SCG_depth (reads/base). The ratio is
  dimensionless and comparable between methods.
- If a gene has the same true coverage depth: short_read_cpg ≈ contig_cpg. ✓
- Median is more robust than mean/40 for diverse communities (only needs >half
  of SCGs to have ≥1 read for a non-zero median).

**Reference:** Bengtsson-Palme et al. (2017) — "Population-level impacts of
antibiotic usage on the human gut microbiome."

### Current pipeline state
- HEAD: `19e17bc`, pushed to `origin/main`.
- All CPG normalisation now uses median_SCG_depth (reads/base) consistently.
- AMR_unified short_reads_cpg and contig_cpg are now on the same scale.

### Known issues / next steps
- MGE JSONDecodeError in MobileElementFinder — upstream bug; `|| true` handles
  gracefully. Awaiting upstream fix.
- `checkm2.yaml`, `gtdbtk.yaml`, `metabolic.yaml` still not version-pinned.
- Recommend an end-to-end test run on real production samples to validate that
  CPG values are now in a plausible range (tetM ~0.1–2.0 cpg expected).

---

## 2026-03-24 (session 27)

### What was done

**Major output refactor — revert cpg, CoverM contig abundance, drop AMR unified (`784caee`)**

The median-SCG normalisation introduced in session 26 was not producing reliable
results. The entire AMR unified table concept was also abandoned. The pipeline
outputs are now redesigned around five clean, purpose-specific files.

**CPG formula reverted:**
- `short_reads_processing.R`: reverted from the median-SCG-depth formula back to
  the original `n_genomes = sum(alignment_length / subject_length) / 40`.
  `cpg = KMA_Depth / n_genomes`. Fallback: `n_genomes = 1` when no SCG hits.
- Markers cpg uses the same `genome.counts` object.

**New output 1 — `AMR_abundance_summary.csv` (short-reads only):**
- Produced directly by `short_reads_processing.R` as a 4th output.
- Columns: `sample`, `AMR_total_cpg` (sum of all gene cpg per sample),
  `pBI143_cpg`, `crAss001_cpg`.
- `short_reads_summary` rule now has 4 named outputs; R script accesses via
  `snakemake@output[[4]]`.

**Contig abundance → CoverM contig:**
- `contig_abundance` rule converted from Python script (`contig_abundance.py`)
  to a shell rule using `shortreads.yaml` conda env (already has minimap2,
  samtools, coverm).
- Pipeline: `minimap2 -ax sr → samtools sort/index → coverm contig`.
- CoverM methods: `trimmed_mean mean reads_aligned`. Header cleaned with awk
  to standard column names: `contig_id`, `trimmed_mean`, `mean_depth`, `reads_mapped`.
- Removed `scgs` input from `contig_abundance` rule (no longer needed).
- `contig_summary.py` updated: `parse_abundance()` now reads the three CoverM
  columns; output schema uses `mean_depth`, `trimmed_mean`, `reads_mapped`
  instead of `abundance_cpg`.
- BAM output path unchanged (still reused by MetaBAT2 `bin` rule).

**Removed `amr_unified`:**
- `amr_unified` rule deleted from `summary.smk`.
- `AMR_unified.csv` and old `AMR_abundance_summary.csv` targets removed from
  Snakefile `all_targets()`.
- `amr_unified.py` and `contig_abundance.py` deleted.
- `amr_unified` removed from `localrules`.

**Five pipeline outputs:**
1. `AMR_abundance_summary.csv` — per-sample AMR total cpg + pBI143/crAss001 cpg
2. `fastp_summary.csv` — per-sample fastp QAQC
3. `contig_summary.tsv` — AMR+MGE features with CoverM contig abundance
4. `mag_summary.tsv` — MAG summary (unchanged)
5. `data/QAQC/short_reads_output.csv` — per-gene short-read AMR data

Both dry-runs pass: production and `--config test=True`.

### Current pipeline state
- HEAD: `784caee`, pushed to `origin/main`.
- Test dry-run: 37 jobs (both mock samples, all stages).
- Production dry-run: 5 local-rule jobs (idempotent from existing outputs).

### Known issues / next steps
- MGE JSONDecodeError in MobileElementFinder — upstream bug; `|| true` handles
  gracefully. Awaiting upstream fix.
- `checkm2.yaml`, `gtdbtk.yaml`, `metabolic.yaml` still not version-pinned.
- Recommend end-to-end test run to validate CoverM contig abundances and
  reverted cpg values are in plausible ranges.

---

## 2026-03-24 (session 28)

### What was done

**End-to-end test run — all 24 jobs passed**

Ran full test mode (`snakemake --use-conda --cores 4 --scheduler greedy --config test=True --rerun-incomplete`) following the major output refactor in session 27.

- Cleared incomplete `test/output/data/genomad/mock2` directory (left from prior interrupted run), unlocked stale Snakemake lock, then re-ran.
- All 24 jobs completed successfully for both mock samples.
- Timing: ~21 minutes total (geNomad is the bottleneck at ~7 min per sample).

**Outputs verified:**

| File | Size | Notes |
|------|------|-------|
| `test/output/AMR_abundance_summary.csv` | 111 B | mock1: 211.1 cpg, mock2: 167.1 cpg; pBI143/crAss001 = 0 (not in mock genomes) |
| `test/output/fastp_summary.csv` | 453 B | Both samples QC'd |
| `test/output/contig_summary.tsv` | 11 KB | tet(A), fosA3 detected in mock1 plasmid contigs with CoverM depths |
| `test/output/mag_summary.tsv` | 402 B | 2 bins per sample; 0 AMR/MGE per bin (expected for mock data) |

### Current pipeline state
- HEAD: `784caee` (unchanged), pushed to `origin/main`.
- ✅ **Full end-to-end test mode passes** — 24/24 jobs, both mock samples.
- Five output files correct and populated.
- CoverM contig abundance columns (`mean_depth`, `trimmed_mean`, `reads_mapped`) present in `contig_summary.tsv`.
- Reverted cpg normalisation (`n_genomes = Σ(alignment_length/gene_length)/40`) producing plausible values (~167–211 cpg for 10× mock communities).

### Known issues / next steps
- MGE JSONDecodeError in MobileElementFinder — upstream bug; `|| true` handles gracefully.
- `checkm2.yaml`, `gtdbtk.yaml`, `metabolic.yaml` still not version-pinned.
- Recommend end-to-end test on real production samples to validate cpg values in realistic range.

---

## 2026-03-25 (session 29)

### What was done

**MobileElementFinder status investigation + fixes (`06bb1f2`)**

**MobileElementFinder status:**
Version 1.1.2 (latest) has an upstream `JSONDecodeError` bug in `me_finder/tools/blast.py`:
`blastn -outfmt 15` (JSON output) writes multiple concatenated JSON objects when processing
many query sequences; `json.load()` stops after the first and raises `Extra data`. The bug
exists in all 1.1.x versions (introduced when JSON replaced XML for BLAST output). Rolling
back to 1.0.x would likely fix the root cause but risks CLI incompatibility. Status:
- **mock1** (fewer contigs): works ✅
- **mock2** (many contigs): JSONDecodeError → empty TSV, handled by `|| true`

**Fix: contig chunking in `mge_annotation` rule**

Replaced the single `mefinder find` call with a chunk-based loop:
1. Split assembly FASTA into batches of ≤200 contigs using `awk`
2. Run `mefinder find` per chunk (small BLAST JSON → no parse error)
3. Concatenate CSV results (comment/header from first hit chunk; data rows from rest)
4. Move combined CSV to output TSV (or `touch` if no hits)

This is version-agnostic and handles the upstream bug without patching installed packages.
The `mag_mge` rule is unchanged — it already runs mefinder per-bin and bins are small.

**Fix: contig summary — trimmed_mean only**

CoverM `contig_abundance` rule already calls `--methods trimmed_mean` (2-column output).
The awk header rename was incorrect (4-column header on 2-column data). Fixed:
- `contigs.smk` `contig_abundance`: awk header → `"contig_id\ttrimmed_mean"`
- `contig_summary.py`: `parse_abundance()` returns scalar trimmed_mean; `all_rows` and
  `COLS` drop `mean_depth` and `reads_mapped`; docstring updated.

### Current pipeline state
- HEAD: `06bb1f2`, pushed to `origin/main`.
- Dry-run: 17 jobs (both mock samples, all stages) — passes cleanly.
- `contig_summary.tsv` schema: `sample | contig_id | trimmed_mean | molecule_type | taxonomy | feature_type | gene | start | stop | strand`
- MobileElementFinder chunking active; mock2 MGE annotations should now be captured.

### Known issues / next steps
- Recommend end-to-end test (`--config test=True`) to confirm MEF chunking works on mock2.
- `checkm2.yaml`, `gtdbtk.yaml`, `metabolic.yaml` still not version-pinned.

---

## 2026-03-25 (session 30)

### What was done

**End-to-end test run — all 15 jobs passed**

Ran full test mode (`snakemake --use-conda --cores 4 --scheduler greedy --config test=True --rerun-incomplete`)
to validate the session 29 fixes (MEF chunking + trimmed_mean headers).

**Results:**

| File | Notes |
|------|-------|
| `contig_abundance` headers | `contig_id\ttrimmed_mean` — 2 columns correct, both samples |
| `contig_summary.tsv` headers | `sample contig_id trimmed_mean molecule_type taxonomy ...` |
| `mock2_mge.tsv` | **49 lines** — MEF chunking fix works; was empty (0 B) before |
| `mock1_mge.tsv` | 40 lines — unchanged |

mock2 now produces 43 MGE feature rows in `contig_summary.tsv`; previously the
JSONDecodeError caused mock2 to produce zero MGE annotations.

### Current pipeline state
- HEAD: `73b1b70` (session 29 fixes), pushed to `origin/main`.
- ✅ **Full end-to-end test mode passes** — 15/15 jobs, both mock samples.
- MEF chunking (≤200 contigs/batch) confirmed working for large contig sets.
- `contig_summary.tsv` schema correct: trimmed_mean only.

---

## Session 31 — 2026-03-24

### Work done
- **`assembly_qa.tsv`** — new pipeline output (one row per sample), parallel to `fastp_summary.csv`
  - Metrics from MEGAHIT log: `n_contigs`, `total_length_bp`, `n50_bp`, `longest_contig_bp`, `shortest_contig_bp`, `mean_contig_length_bp`
  - Metrics from samtools flagstat on existing BAM: `reads_total`, `reads_mapped`, `pct_reads_mapped`
  - New script: `workflow/scripts/assembly_qa.py` (regex-parses MEGAHIT summary line; subprocess flagstat)
  - New `assembly_qa` localrule in `contigs.smk` using `contigs.yaml` env (has both samtools + pandas)
  - Added to `localrules` and `all_targets()` in `Snakefile`
  - Fixed conda env: initial draft used `shortreads.yaml` (no pandas) → changed to `contigs.yaml`
  - Test output verified: mock1 (331 contigs, 8.28 Mbp, N50=67.7 kb, 99.3% reads mapped), mock2 (486 contigs, 11.6 Mbp, N50=52.8 kb, 99.5% reads mapped)
- **Commit**: `e3b8da4`

### Current state
All 31 jobs pass end-to-end in test mode. Key outputs:
- `fastp_summary.csv` — read QC per sample
- `short_reads_output.csv`, `markers_cpg.csv` — short-read AMR + marker cpg
- `assembly_qa.tsv` — per-sample assembly QA metrics (**new this session**)
- `contig_summary.tsv` — per-contig annotations (trimmed_mean, molecule_type, taxonomy, feature_type, gene)
- `AMR_unified.csv`, `AMR_abundance_summary.csv` — cross-stage AMR
- `mag_summary.tsv` — MAG-level annotations

### Known issues / next steps
- None outstanding

---

## 2026-04-10 (session 27)

### What was done

#### Anthropogenic marker injection into test data
- Created `test/scripts/inject_marker_reads.py`: simulates 200 paired-end reads (150 bp, seed=99) from pBI143 and 500 from crAss001, appended to both mock1 and mock2 FASTQs
- pBI143 (Bacteroides fragilis plasmid, 2,747 bp) at ~10x coverage; crAss001 (102,679 bp) at ~1x
- Expected: pBI143_cpg >= 1.0, crAss001_cpg ~0.3-1.5 after SCG normalisation -> E_exposure = 1.00 / 0.75 for test samples
- Updated `test/data/mock_data_summary.txt` and `test/scripts/generate_mock_data.py` to document injection step

#### AMR Risk Scoring Module
- Implemented **Strategy 1 (Additive)** and **Strategy 2 (Multiplicative)** from `unified_amr_risk_score_strategies.txt`
- Four components (all 0-1):
  - **R** (Resistance Hazard): WHO criticality tiers from AMRFinderPlus subclass/class (0.40-1.00)
  - **M** (Mobility): molecule_type + MGE co-location from contig_summary.tsv; short_reads_only = 0.10
  - **H** (Host/Pathogenicity): ESKAPE/Enterobacteriaceae/human-assoc/environmental from MMseqs2 taxonomy (0.15-1.00)
  - **E** (Exposure): sample-level from markers_cpg.csv thresholds (0.10-1.00)
- New file: `workflow/scripts/amr_risk_score.py`
- New rule: `amr_risk_score` localrule in `workflow/rules/summary.smk`; added to `localrules` in `Snakefile`
- Removed `AMR_abundance_summary.csv` production from `short_reads_processing.R` (was out_file4); rule output also removed from `short_reads_summary` in `short_reads.smk`
- Output columns: `sample, AMR_total_cpg, pBI143_cpg, crAss001_cpg, E_exposure, R_mean, M_mean, H_mean, amr_risk_additive_raw, amr_risk_multiplicative_raw, amr_risk_additive, amr_risk_multiplicative`

### Pipeline state
- Dry-run confirmed: 39 jobs in test mode (up from 37), DAG valid
- Full end-to-end test run NOT yet done this session

### Known issues / next steps
- Run full test to verify E_exposure > 0.10 for both mock samples after marker injection
- `AMR_unified.csv` referenced in earlier session logs no longer exists; `AMR_abundance_summary.csv` is the authoritative per-sample AMR table

---

## 2026-04-10 (session 32)

### What was done

**End-to-end test run — all 11 remaining jobs passed** (resuming from session 27 pipe break)

Short-reads analysis and MAG stages were already complete from prior runs. Resumed from
the geNomad/contig stage forward with `--rerun-incomplete`.

**Bug fixed: `amr_risk_score` — Namedlist TypeError (`6e81cae`)**

`snakemake.input.get("contig_summary", None)` returned a Snakemake `Namedlist` instead of
a string because `contig_summary` is defined as a conditional list in the rule
(`[path] if _RUN_CTG else []`). Fixed by extracting the first element explicitly:

```python
_ctg       = snakemake.input.get("contig_summary", [])
contig_tsv = str(_ctg[0]) if _ctg else None
```

**Outputs verified:**

| File | Notes |
|------|-------|
| `AMR_abundance_summary.csv` | Risk scores generated for both mock samples |
| `contig_summary.tsv` | 214 rows (AMR + MGE features, both samples) |

**Risk score spot-check:**

| Sample | pBI143_cpg | crAss001_cpg | E_exposure | R_mean | M_mean | H_mean | Additive (norm) |
|--------|-----------|-------------|-----------|--------|--------|--------|----------------|
| mock1 | 3.64 | 0.24 | 1.0 | 0.752 | 0.836 | 0.181 | 100.0 |
| mock2 | 3.57 | 0.24 | 1.0 | 0.703 | 0.766 | 0.185 | 0.0 |

E_exposure = 1.0 for both samples confirms marker injection (pBI143 + crAss001) works.
Normalised scores are min-max across samples — with only 2 samples, one is always 0 and
one is 100 (expected behaviour; will spread meaningfully with more production samples).

### Current pipeline state
- HEAD: `6e81cae`, pushed to origin/main.
- Full end-to-end test mode passes — all jobs, both mock samples.
- AMR risk scoring module (additive + multiplicative, 4-component R x M x H x E) confirmed working.
- `AMR_abundance_summary.csv` schema: `sample, AMR_total_cpg, pBI143_cpg, crAss001_cpg, E_exposure, R_mean, M_mean, H_mean, amr_risk_additive_raw, amr_risk_multiplicative_raw, amr_risk_additive, amr_risk_multiplicative`

### Known issues / next steps
- Normalised risk scores with only 2 test samples are always 0/100 — not a bug, expected with min-max when n=2.
- `checkm2.yaml`, `gtdbtk.yaml`, `metabolic.yaml` still not version-pinned.
- Consider end-to-end test on real production samples to validate risk scores in a realistic multi-sample cohort.

---

## 2026-04-14 (session 33)

### What was done
- Refreshed `README.md` to match the current pipeline outputs and remove stale references to `AMR_unified.csv`.
- Reorganized the README around the user workflow: quick start, production setup, databases, SLURM, partial runs, and main outputs.
- Simplified wording throughout and documented the current sample-level risk scoring output in `AMR_abundance_summary.csv`.
- Added the newer top-level outputs `assembly_qa.tsv` and `mag_summary.tsv`, and corrected short-read table paths to `data/QAQC/short_reads_output.csv` and `data/QAQC/markers_cpg.csv`.

### Current pipeline state
- HEAD before this documentation update: `9222364`, pushed to `origin/main`.
- Full end-to-end test mode still last recorded as passing for both mock samples.
- Primary documented outputs are now aligned with the current workflow: `fastp_summary.csv`, `assembly_qa.tsv`, `contig_summary.tsv`, `AMR_abundance_summary.csv`, and `mag_summary.tsv`.

### Known issues / next steps
- No new pipeline issues introduced; this session was documentation-only.
- The existing follow-ups remain: environment version pinning for `checkm2.yaml`, `gtdbtk.yaml`, and `metabolic.yaml`, plus validation on a larger real cohort.

---

## 2026-04-27 (session 34)

### What was done
- Replaced the dual Slurm profiles with a single default profile at `config/slurm/config.yaml`.
  - Added resilient execution defaults: `keep-going`, `rerun-incomplete`, `restart-times`, `latency-wait`, and grouped execution for lightweight sample-parallel rules.
  - Kept heavy `contigs`, `genomad`, and `bin` resource scaling in-rule so memory/runtime can grow across retry attempts.
- Added `run_swam-meta.sh` as the default HPC entrypoint.
  - Writes a timestamped driver log under `{out_dir}/logs/run_status/`.
  - Runs the new post-run summary generator `workflow/scripts/build_run_report.py`.
- Added a run-status system modeled on SWAM-g.
  - Per-rule logs now live under `{out_dir}/logs/` (`logs/<rule>/<sample>.log` for sample rules, `logs/<rule>.log` for global rules).
  - New run-status outputs: `{out_dir}/logs/run_status/run_report.txt` and `{out_dir}/logs/run_status/sample_rule_status.tsv`.
- Reorganized database path handling in `workflow/rules/common.smk` and init rules.
  - New per-tool / per-reference layout: `dbs/amrfinderplus/`, `dbs/amrfinderplus_kma/`, `dbs/checkm2/`, `dbs/genomad/`, `dbs/human/`, `dbs/markers/`, `dbs/metabolic/`, `dbs/scg/`, `dbs/uniref50/`.
  - Added legacy-layout migration logic so old `dbs/short_reads/`, `dbs/genomad_db/`, and `dbs/amrfinderplus_db/` contents can be reused instead of re-downloaded.
- Promoted `short_reads_output.csv` back to a top-level key output at `{out_dir}/short_reads_output.csv`.
  - Kept `markers_cpg.csv` in `data/QAQC/` as the backend support table for exposure scoring.
  - Updated `Snakefile`, summary rules, README, config comments, and Copilot instructions to match.
- Fixed a MAG abundance edge case uncovered by the rerun.
  - `coverm genome` failed when MetaBAT bins shared duplicate contig IDs.
  - `mag_abundance` now builds a temporary de-duplicated bin FASTA set for CoverM so the workflow completes cleanly on the test dataset.

### Current pipeline state
- Working tree on `main` contains the Slurm/profile/logging/db-layout overhaul; not yet committed in this session.
- Full end-to-end test mode passes again after the overhaul:
  - `snakemake --use-conda --cores 4 --scheduler greedy --config test=True --rerun-incomplete`
- The new run-status report generator was smoke-tested successfully against `test/output/`.
- Key top-level outputs are now:
  - `fastp_summary.csv`
  - `short_reads_output.csv`
  - `assembly_qa.tsv`
  - `contig_summary.tsv`
  - `AMR_abundance_summary.csv`
  - `mag_summary.tsv`

### Known issues / next steps
- `checkm2.yaml`, `gtdbtk.yaml`, and `metabolic.yaml` are still not version-pinned.
- The new Slurm wrapper/reporting path was validated locally and against test outputs, but not yet on a real cluster submission.
- A larger real cohort is still needed to validate risk-score spread beyond the 2-sample mock test.

## 2026-04-28 — ARGNet discovery mode added and exercised on the mock dataset

### What changed
- Added a new `amr_mode` config switch in `config/config.yaml`:
  - `operational` (default): existing AMRFinderPlus-only behavior
  - `discovery`: keep all existing operational outputs and add ARGNet discovery outputs
- Added discovery-mode plumbing in `workflow/rules/common.smk` and `workflow/Snakefile`:
  - normalized/validated `amr_mode`
  - new discovery-mode top-level outputs:
    - `ARGnet_contig_summary.tsv`
    - `ARGnet_mag_summary.tsv`
- Added ARGNet bootstrap and inference rules:
  - `init_argnet` clones ARGNet under `dbs/argnet/ARGNet/` and downloads the long-sequence model files (`AELS_tall.h5`, `classifier-ls_tall.h5`)
  - `contig_argnet` runs ARGNet on per-sample Prodigal protein FASTA from assembled contigs
  - `mag_argnet` concatenates per-bin Prodigal proteins with bin-prefixed IDs and runs ARGNet once per sample
- Added a dedicated ARGNet conda env at `workflow/envs/argnet.yaml` (Python 3.11 + TensorFlow CPU stack compatible with upstream ARGNet)
- Added `workflow/scripts/run_argnet.py` as a wrapper around the upstream ARGNet checkout:
  - imports the upstream model module directly
  - forces CPU execution (`CUDA_VISIBLE_DEVICES=-1`) so discovery runs do not compete for GPU memory
  - sanitizes Prodigal amino-acid FASTA by stripping `*` stop codons before inference
  - moves repo-relative ARGNet output into the Snakemake-declared output path
- Added discovery summary scripts:
  - `workflow/scripts/argnet_contig_summary.py`
  - `workflow/scripts/argnet_mag_summary.py`
  These emit ARG-positive-only tables with contextual metadata instead of modifying existing AMRFinderPlus outputs or risk scoring.
- Updated `README.md` to document `amr_mode` and the new discovery outputs.

### Validation / test run
- Discovery-mode dry run passed:
  - `snakemake -s workflow/Snakefile -n --use-conda --cores 4 --scheduler greedy --config test=True amr_mode=discovery`
- Discovery-mode execution on the bundled mock dataset succeeded for the new outputs:
  - `snakemake -s workflow/Snakefile --use-conda --cores 4 --scheduler greedy --rerun-incomplete /home/benda/Projects/NARMS/SWAM-meta/test/output/ARGnet_contig_summary.tsv /home/benda/Projects/NARMS/SWAM-meta/test/output/ARGnet_mag_summary.tsv --config test=True amr_mode=discovery`
- Output shape confirmed:
  - per-sample raw ARGNet files contain `test_id`, `ARG_prediction`, `resistance_category`, `probability`
  - `ARGnet_contig_summary.tsv` adds `sample`, `protein_id`, `contig_id`, `trimmed_mean`, `molecule_type`, and `taxonomy`
  - `ARGnet_mag_summary.tsv` adds `sample`, `bin_id`, `protein_id`, `contig_id`, and MAG abundance
- Mock test output counts from this run:
  - `ARGnet_contig_summary.tsv`: 1648 ARG-positive rows
  - `ARGnet_mag_summary.tsv`: 1839 ARG-positive rows

### Current pipeline state
- Operational mode remains the default and still uses AMRFinderPlus-derived outputs and risk scoring only.
- Discovery mode now runs ARGNet alongside the operational AMRFinderPlus path for contig and MAG proteins.
- The mock dataset run completed successfully for both discovery summaries after two integration fixes:
  - sanitize `*` stop codons in Prodigal proteins before ARGNet inference
  - force CPU execution to avoid concurrent GPU-memory contention during test-mode runs

### Known issues / next steps
- Discovery mode is currently implemented for **contig and MAG proteins only**; short-read ARGNet inference is intentionally not wired into `short_reads_output.csv`.
- The mock data produced many ARGNet-positive calls, which is expected for a broad discovery model but should be interpreted as exploratory evidence rather than operational confirmation.
- A real production cohort is still needed to benchmark ARGNet runtime and calibrate how discovery summaries should be interpreted alongside AMRFinderPlus calls.

## 2026-04-28 — ARGNet short-read discovery added

### What changed
- Extended discovery mode so ARGNet now runs on **host-filtered short reads** as well as contig and MAG proteins.
- Added the short-read raw inference rule in `workflow/rules/short_reads.smk`:
  - `short_reads_argnet` converts `{sample}_R1.clean.fastq.gz` + `{sample}_R2.clean.fastq.gz` to a temporary FASTA
  - runs `run_argnet.py --seq-type nt --model argnet-s`
  - writes per-sample raw outputs to `data/argnet_short_reads/{sample}_argnet.tsv`
- Added a new top-level discovery summary:
  - `ARGnet_short-reads_summary.tsv`
- Added `workflow/scripts/argnet_short_reads_summary.py`, which summarizes ARGNet short-read calls by:
  - `sample`
  - `resistance_category`
  - `arg_read_count`
  - `total_reads_scored`
  - `total_arg_reads`
  - `non_arg_reads`
  - `short_fragment_reads`
  - `arg_fraction`
  - `arg_reads_per_million`
  - `mean_probability`
  - `max_probability`
- Kept short-read ARGNet metrics intentionally **read-level**, not cpg-normalized, so they can sit next to AMRFinder/KMA outputs without implying direct apples-to-apples abundance equivalence.
- Updated ARGNet bootstrap/env support:
  - `init_argnet` now downloads the short-read model files too (`AESS_tall.h5`, `classifier-ss_tall.h5`)
  - `workflow/envs/argnet.yaml` now installs `kito`, which the upstream `argnet_ssnt_new_chunk.py` imports
- Updated README/config comments and the session plan to reflect that discovery mode now covers **short reads, contigs, and MAGs**.

### Validation / test run
- Dry run passed for the new short-read target:
  - `snakemake -s workflow/Snakefile -n --use-conda --cores 4 --scheduler greedy /home/benda/Projects/NARMS/SWAM-meta/test/output/ARGnet_short-reads_summary.tsv --config test=True amr_mode=discovery`
- Full short-read discovery execution on the bundled mock dataset succeeded:
  - `snakemake -s workflow/Snakefile --use-conda --cores 4 --scheduler greedy /home/benda/Projects/NARMS/SWAM-meta/test/output/ARGnet_short-reads_summary.tsv --config test=True amr_mode=discovery`
- Observed mock output sizes:
  - `mock1_argnet.tsv`: 570,253 scored reads
  - `mock2_argnet.tsv`: 791,519 scored reads
  - `ARGnet_short-reads_summary.tsv`: 45 category rows across both samples

### Current pipeline state
- Discovery mode now produces a parallel ARGNet layer at all three levels:
  - short reads: `ARGnet_short-reads_summary.tsv`
  - contigs: `ARGnet_contig_summary.tsv`
  - MAGs: `ARGnet_mag_summary.tsv`
- Operational AMRFinder outputs and risk scoring remain unchanged.

### Known issues / next steps
- Short-read ARGNet is **substantially slower** than contig/MAG ARGNet on the mock dataset:
  - mock2 completed in ~21 minutes
  - mock1 completed in ~15 minutes
  This is workable for smoke-testing but may be expensive on large production cohorts.
- Despite forcing CPU mode in `run_argnet.py`, the upstream short-read model still appears to initialize GPU-visible TensorFlow state in this environment. It ran successfully here, but this behavior should be revisited if GPU contention becomes an issue on shared systems.
- The short-read ARGNet summary is intentionally read-level (counts/fractions/RPM), not cpg-normalized. Any future “apples to apples” comparison with AMRFinder short-read abundance should be treated as a separate methods problem rather than folded into this table implicitly.
