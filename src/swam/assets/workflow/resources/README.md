# workflow/resources

This directory contains small reference files that ship with SWAM.

## SCGs_40_All.fasta

A curated set of 40 universal single-copy genes (SCGs) used to estimate
the number of genomes in a sample (`n_genomes`), which normalises AMR
abundance to copies-per-genome (cpg).

**This file must be present before running SWAM in production mode.**
Copy it here from your working location:

```bash
cp /path/to/SCGs_40_All.fasta workflow/resources/SCGs_40_All.fasta
```

In test mode, the pipeline uses `test/dbs/scg/SCGs_40.fasta` (a mini subset)
and does not require this file.
