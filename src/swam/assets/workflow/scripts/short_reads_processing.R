#!/usr/bin/env Rscript

libraries <- c("tidyverse", "Nonpareil")

invisible(lapply(libraries, function(x) {
  suppressMessages(suppressWarnings(library(x, character.only = T)))
  }))
  
  
# Input files

afp_files       <- snakemake@input[["afp_files"]]
scgs            <- snakemake@input[["scgs"]]
npos            <- snakemake@input[["npos"]]
jsons           <- snakemake@input[["jsons"]]
afp_path        <- snakemake@input[["afp_metadata"]]
markers_files   <- snakemake@input[["markers_files"]]

# Output files

out_file1 <- snakemake@output[[1]]
out_file2 <- snakemake@output[[2]]
out_file3 <- snakemake@output[[3]]


###############################################
###      summarize metagenomic coverage data
###############################################

nps <- as.data.frame(summary(Nonpareil.set(npos)))
nps <- nps %>% rownames_to_column(var="sample")

#####################################################
###    generate summary table from fastp json files
#####################################################

# Robust pluck helper with default
pluck_default <- function(x, ..., .default = NA) {
  tryCatch(
    purrr::pluck(x, ..., .default = .default),
    error = function(e) .default
  )
}

# Summarize a single fastp JSON file into a one-row tibble
summarize_fastp <- function(json_path) {
  j <- jsonlite::fromJSON(json_path, simplifyVector = TRUE)
  
  # fastp's main summary is typically in j$summary$before_filtering and j$summary$after_filtering
  before <- pluck_default(j, "summary", "before_filtering", .default = list())
  after  <- pluck_default(j, "summary", "after_filtering",  .default = list())
  filtering <- pluck_default(j, "summary", "filtering_result", .default = list())
  
  # helper to safely get numeric fields (returns NA if missing)
  get_num <- function(lst, name) {
    val <- pluck_default(lst, name, .default = NA)
    if (is.null(val)) NA_real_ else as.numeric(val)
  }
  
  # Basic fields from before/after
  total_reads_before <- get_num(before, "total_reads")
  total_reads_after  <- get_num(after,  "total_reads")
  total_bases_before <- get_num(before, "total_bases")
  total_bases_after  <- get_num(after,  "total_bases")
  
  q30_before  <- get_num(before,  "q30_rate")
  q30_after   <- get_num(after,  "q30_rate")
  
  read1_len_after  <- get_num(after, "read1_mean_length")
  read2_len_after  <- get_num(after, "read2_mean_length")

  # Derived metrics: retention, percent changes, adapter trimmed percent
  pct_reads_retained <- if (!is.na(total_reads_before) && total_reads_before > 0) total_reads_after / total_reads_before * 100 else NA_real_
  pct_bases_retained <- if (!is.na(total_bases_before) && total_bases_before > 0) total_bases_after / total_bases_before * 100 else NA_real_
  avg_read_length    <-  (read1_len_after + read2_len_after) / 2
  
  tibble(
    sample = tools::file_path_sans_ext(basename(json_path)),
    
    total_reads_before = total_reads_before,
    total_reads_after  = total_reads_after,
    pct_reads_retained = pct_reads_retained,
    
    total_bases_before = total_bases_before,
    total_bases_after  = total_bases_after,
    pct_bases_retained = pct_bases_retained,
    
    q30_before = q30_before,
    q30_after  = q30_after,
    
    avg_read_length = avg_read_length
    
  )
}

fastp_summary <- purrr::map_dfr(jsons, summarize_fastp)

fastp_summary <- fastp_summary %>% left_join(nps, by="sample") %>%
                 rename(metagenomic_coverage=C, sequence_diversity=diversity) %>% 
                 select(sample, metagenomic_coverage, sequence_diversity, 
                        total_reads_before, total_reads_after, pct_reads_retained, 
                        total_bases_before, total_bases_after, pct_bases_retained,
                        q30_before, q30_after, 
                        avg_read_length)


####################################################
###      read in alignment files, normalize to cpg
####################################################

afp_metadata <- read_tsv(afp_path, show_col_types = FALSE)
message("AMRFinder catalog: ", nrow(afp_metadata), " genes loaded from ", afp_path)

catalog_long <- afp_metadata %>% 
                      pivot_longer(cols = c("refseq_nucleotide_accession", "genbank_nucleotide_accession"), values_to = "key") %>%
                      filter(!is.na(key) & key != "") %>%
                      distinct(key, .keep_all = TRUE)

afp <- do.call(bind_rows, lapply(afp_files, function(f) {
  x <- tryCatch(read_tsv(f, col_types = cols()),
                error=function(e) NULL)
  if (is.null(x) || nrow(x) == 0) return(NULL)
  x$sample <- sub(".afp.res", "", basename(f))
  x <- x %>% relocate(sample)
  x <- x %>% separate(`#Template`,
                       into = c("refseq_protein_accession", "refseq_nucleotide_accession",
                                "start", "stop", "allele_kma", "gene_family_kma", "desc_kma"),
                       sep = "\\|", extra = "merge", fill = "right") %>%
             # KMA includes the full FASTA description after the first space in the
             # sequence identifier; strip it so the accession join key is clean.
             mutate(refseq_nucleotide_accession = str_extract(refseq_nucleotide_accession, "^\\S+"))
  x
}))


scg <- do.call(bind_rows, lapply(scgs, function(f) {
  x <- tryCatch(read.table(f),
                error=function(e) NULL)
  if (is.null(x) || nrow(x) == 0) return(NULL)
  colnames(x) <- c("qseqid", "sseqid", "pident", "length", "evalue", "bitscore", "slen")
  x$sample <- sub(".scgs", "", basename(f))
  x <- x %>% relocate(sample)
  x
}))

# ---------------------------------------------------------------------------
# Estimate genome counts from 40 single-copy genes (SCGs).
#
# For each read-SCG alignment:
#   fraction covered = alignment_length / subject_length
# n_genomes = Σ(alignment_length / subject_length) / 40
#
# CPG = KMA_Depth / n_genomes
#
# Reference: Bengtsson-Palme et al. (2017).
# Falls back to n_genomes = 1 when no SCG alignments are detected.
# ---------------------------------------------------------------------------
if (!is.null(scg) && nrow(scg) > 0 && "sample" %in% colnames(scg)) {
  genome.counts <- scg %>%
    group_by(sample) %>%
    distinct(qseqid, .keep_all = TRUE) %>%   # one alignment per read (best hit)
    summarise(n_genomes = sum(length / slen) / 40, .groups = "drop")
} else {
  message("WARNING: No SCG alignments found. Falling back to n_genomes = 1.")
  genome.counts <- tibble(sample = sub(".scgs", "", basename(scgs)), n_genomes = 1)
}

# estimate copies per genome (cpg):  cpg = KMA_Depth / n_genomes
if (!is.null(afp) && nrow(afp) > 0) {
  afp_long <-  afp %>% 
                    left_join(genome.counts, by="sample") %>% 
                    left_join(fastp_summary %>% select(sample, avg_read_length), by="sample") %>%
                    mutate(read_count = round(Depth * Template_length / avg_read_length),
                           cpg = Depth / n_genomes) %>%
                    filter(read_count >= 1) %>% 
                    left_join(catalog_long, by=c("refseq_nucleotide_accession" = "key")) %>%
                    # Production catalog has allele/gene_family empty for many genes;
                    # fall back to values parsed directly from the KMA template header.
                    mutate(
                      allele      = coalesce(allele, allele_kma),
                      gene_family = coalesce(gene_family, gene_family_kma)
                    ) %>%
                    mutate(allele_pass = ifelse(Template_Identity >= 80 & Query_Identity >= 99 & !is.na(allele), "yes", "no")) %>%
                    select(sample, read_count, Depth, cpg, allele, allele_pass, gene_family, product_name, 
                           scope, type, subtype, class, subclass, refseq_nucleotide_accession)
} else {
  afp_long <- tibble(sample=character(), read_count=integer(), Depth=double(), cpg=double(),
                     allele=character(), allele_pass=character(), gene_family=character(),
                     product_name=character(), scope=character(), type=character(),
                     subtype=character(), class=character(), subclass=character(),
                     refseq_nucleotide_accession=character())
}


###############################################
###       write csv outputs
###############################################

write.csv(fastp_summary, out_file1, row.names=F)
write.csv(afp_long, out_file2, row.names=F)


###############################################
###   anthropogenic markers (pBI143, crAss001)
###############################################

markers <- do.call(bind_rows, lapply(markers_files, function(f) {
  x <- tryCatch(read_tsv(f, col_types = cols()),
                error = function(e) NULL)
  if (is.null(x) || nrow(x) == 0) return(NULL)
  x$sample <- sub("\\.markers\\.res$", "", basename(f))
  x <- x %>% relocate(sample) %>% rename(template = `#Template`)
  x
}))

if (!is.null(markers) && nrow(markers) > 0) {
  marker_cpg <- markers %>%
    # KMA includes the full FASTA description after the first space in the template
    # identifier; extract only the accession to allow reliable matching.
    mutate(seq_id = str_extract(template, "^\\S+")) %>%
    left_join(genome.counts, by = "sample") %>%
    left_join(fastp_summary %>% select(sample, avg_read_length), by = "sample") %>%
    mutate(cpg = Depth / n_genomes) %>%
    select(sample, seq_id, cpg)
} else {
  marker_cpg <- tibble(sample = character(), seq_id = character(), cpg = double())
}

###############################################
###   markers cpg output (pBI143 + crAss001)
###############################################

pbi143_cpg <- marker_cpg %>%
  filter(seq_id == "U30316.1") %>%
  select(sample, pBI143_cpg = cpg)

crass_cpg <- marker_cpg %>%
  filter(seq_id == "NC_049977.1") %>%
  select(sample, crAss001_cpg = cpg)

markers_out <- tibble(sample = unique(fastp_summary$sample)) %>%
  left_join(pbi143_cpg, by = "sample") %>%
  left_join(crass_cpg, by = "sample") %>%
  mutate(
    pBI143_cpg  = replace_na(pBI143_cpg, 0),
    crAss001_cpg = replace_na(crAss001_cpg, 0)
  )

write.csv(markers_out, out_file3, row.names = FALSE)
