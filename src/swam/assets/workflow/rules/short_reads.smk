# ------------------------------------------
#         initiate databases
# ------------------------------------------

rule initiate_dbs:
    output:
        afp_db         = os.path.join(_DB_AMR_KMA_DIR, f".done{_DB_TAG}"),
        scg_db         = os.path.join(_DB_SCG_DIR, f".done{_DB_TAG}"),
        afp_metadata   = os.path.join(_DB_AMR_KMA_DIR, "ReferenceGeneCatalog.txt"),
        h_genome_done  = os.path.join(_DB_HUMAN_DIR, f".done{_DB_TAG}"),
        markers_db     = os.path.join(_DB_MARKERS_DIR, f".done{_DB_TAG}"),
        metabolic_done = os.path.join(_DB_METABOLIC_DIR, f".done{_DB_TAG}")
    params:
        amr_kma_dir      = _DB_AMR_KMA_DIR,
        human_dir        = _DB_HUMAN_DIR,
        markers_dir      = _DB_MARKERS_DIR,
        scg_dir          = _DB_SCG_DIR,
        metabolic_root   = _DB_METABOLIC_DIR,
        h_genome_path    = os.path.join(_DB_HUMAN_DIR, "GCF_000001405.40_GRCh38.p14_genomic.fna.gz"),
        scg_db           = _SCG_DB,
        scg_source       = _SCG_SOURCE_FASTA,
        skip_metabolic   = config.get("skip_metabolic", False),
        metabolic_dir    = _METABOLIC_DIR,
        legacy_root      = _DBS_DIR,
        legacy_sr_dir    = _LEGACY_SR_DBS_DIR,
        test_mode        = _TEST,
        test_amr_fa      = os.path.join(_TEST_DIR, "dbs", "amrfinder", "AMR_CDS.fa"),
        test_meta        = os.path.join(_TEST_DIR, "dbs", "amrfinder", "ReferenceGeneCatalog.txt"),
        test_human       = os.path.join(_TEST_DIR, "dbs", "human", "human_mini.fna.gz"),
        test_markers_dir = os.path.join(_TEST_DIR, "dbs", "markers"),
    log:
        rule_log("initiate_dbs")
    conda: "../envs/shortreads.yaml"
    shell:
        """
        mkdir -p "$(dirname {log})" \
                 {params.amr_kma_dir} \
                 {params.human_dir} \
                 {params.markers_dir} \
                 {params.scg_dir} \
                 {params.metabolic_root}

        # Migrate legacy dbs/ layouts into per-tool subdirectories when needed.
        if [ -d "{params.legacy_sr_dir}" ]; then
            [ -f "{params.legacy_sr_dir}/AMR_CDS.fa" ] && mv -n "{params.legacy_sr_dir}/AMR_CDS.fa" "{params.amr_kma_dir}/AMR_CDS.fa" || true
            [ -f "{params.legacy_sr_dir}/ReferenceGeneCatalog.txt" ] && mv -n "{params.legacy_sr_dir}/ReferenceGeneCatalog.txt" "{params.amr_kma_dir}/ReferenceGeneCatalog.txt" || true
            [ -f "{params.legacy_sr_dir}/afp_db.name" ] && mv -n "{params.legacy_sr_dir}"/afp_db.* "{params.amr_kma_dir}/" || true
            [ -f "{params.legacy_sr_dir}/GCF_000001405.40_GRCh38.p14_genomic.fna.gz" ] && mv -n "{params.legacy_sr_dir}/GCF_000001405.40_GRCh38.p14_genomic.fna.gz" "{params.human_dir}/" || true
            [ -f "{params.legacy_sr_dir}/markers.fa" ] && mv -n "{params.legacy_sr_dir}/markers.fa" "{params.markers_dir}/markers.fa" || true
            [ -f "{params.legacy_sr_dir}/pBI143.fasta" ] && mv -n "{params.legacy_sr_dir}/pBI143.fasta" "{params.markers_dir}/" || true
            [ -f "{params.legacy_sr_dir}/crAss001.fasta" ] && mv -n "{params.legacy_sr_dir}/crAss001.fasta" "{params.markers_dir}/" || true
            [ -f "{params.legacy_sr_dir}/markers_db.name" ] && mv -n "{params.legacy_sr_dir}"/markers_db.* "{params.markers_dir}/" || true
            [ -f "{params.legacy_sr_dir}/scg_db.dmnd" ] && mv -n "{params.legacy_sr_dir}/scg_db.dmnd" "{params.scg_dir}/" || true
            [ -d "{params.legacy_sr_dir}/METABOLIC" ] && mv -n "{params.legacy_sr_dir}/METABOLIC" "{params.metabolic_root}/" || true
        fi

        [ -f "{params.legacy_root}/AMR_CDS.fa" ] && mv -n "{params.legacy_root}/AMR_CDS.fa" "{params.amr_kma_dir}/AMR_CDS.fa" || true
        [ -f "{params.legacy_root}/ReferenceGeneCatalog.txt" ] && mv -n "{params.legacy_root}/ReferenceGeneCatalog.txt" "{params.amr_kma_dir}/ReferenceGeneCatalog.txt" || true
        [ -f "{params.legacy_root}/afp_db.name" ] && mv -n "{params.legacy_root}"/afp_db.* "{params.amr_kma_dir}/" || true
        [ -f "{params.legacy_root}/GCF_000001405.40_GRCh38.p14_genomic.fna.gz" ] && mv -n "{params.legacy_root}/GCF_000001405.40_GRCh38.p14_genomic.fna.gz" "{params.human_dir}/" || true
        [ -f "{params.legacy_root}/markers.fa" ] && mv -n "{params.legacy_root}/markers.fa" "{params.markers_dir}/markers.fa" || true
        [ -f "{params.legacy_root}/pBI143.fasta" ] && mv -n "{params.legacy_root}/pBI143.fasta" "{params.markers_dir}/" || true
        [ -f "{params.legacy_root}/crAss001.fasta" ] && mv -n "{params.legacy_root}/crAss001.fasta" "{params.markers_dir}/" || true
        [ -f "{params.legacy_root}/markers_db.name" ] && mv -n "{params.legacy_root}"/markers_db.* "{params.markers_dir}/" || true
        [ -f "{params.legacy_root}/scg_db.dmnd" ] && mv -n "{params.legacy_root}/scg_db.dmnd" "{params.scg_dir}/" || true
        [ -d "{params.legacy_root}/METABOLIC" ] && mv -n "{params.legacy_root}/METABOLIC" "{params.metabolic_root}/" || true

        if [ "{params.test_mode}" = "True" ]; then

            echo "initiate_dbs: indexing AMRFinderPlus database (test)..."
            cp {params.test_amr_fa} {params.amr_kma_dir}/AMR_CDS.fa
            kma index -i {params.amr_kma_dir}/AMR_CDS.fa \
                      -o {params.amr_kma_dir}/afp_db >> {log} 2>&1
            touch {output.afp_db}
            cp {params.test_meta} {output.afp_metadata}

            echo "initiate_dbs: building SCG DIAMOND database (test)..."
            diamond makedb --in {params.scg_source} \
                           -d {params.scg_db} --quiet >> {log} 2>&1
            touch {output.scg_db}
            cp {params.test_human} {params.h_genome_path}
            touch {output.h_genome_done}

            echo "initiate_dbs: indexing anthropogenic markers (test)..."
            cat {params.test_markers_dir}/pBI143.fasta \
                {params.test_markers_dir}/crAss001.fasta \
                > {params.markers_dir}/markers.fa
            cp {params.test_markers_dir}/pBI143.fasta {params.markers_dir}/pBI143.fasta
            cp {params.test_markers_dir}/crAss001.fasta {params.markers_dir}/crAss001.fasta
            kma index -i {params.markers_dir}/markers.fa \
                      -o {params.markers_dir}/markers_db >> {log} 2>&1
            touch {output.markers_db}

            touch {output.metabolic_done}

        else

            if [ ! -f "{params.amr_kma_dir}/AMR_CDS.fa" ]; then
                echo "initiate_dbs: downloading AMRFinderPlus CDS..."
                wget --show-progress -P {params.amr_kma_dir} https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/latest/AMR_CDS.fa 2>&1 | tee -a {log}
            else
                echo "initiate_dbs: AMR_CDS.fa already present, skipping download"
            fi
            if [ ! -f "{output.afp_metadata}" ]; then
                echo "initiate_dbs: downloading ReferenceGeneCatalog.txt..."
                wget --show-progress -P {params.amr_kma_dir} https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/latest/ReferenceGeneCatalog.txt 2>&1 | tee -a {log}
            else
                echo "initiate_dbs: ReferenceGeneCatalog.txt already present, skipping"
            fi

            if [ ! -f "{params.amr_kma_dir}/afp_db.name" ]; then
                echo "initiate_dbs: indexing AMRFinderPlus database..."
                kma index -i {params.amr_kma_dir}/AMR_CDS.fa \
                          -o {params.amr_kma_dir}/afp_db >> {log} 2>&1
            else
                echo "initiate_dbs: KMA AFP index already present, skipping"
            fi
            touch {output.afp_db}

            if [ ! -f "{params.scg_dir}/scg_db.dmnd" ]; then
                echo "initiate_dbs: building SCG DIAMOND database..."
                diamond makedb --in {params.scg_source} \
                               -d {params.scg_db} --quiet >> {log} 2>&1
            else
                echo "initiate_dbs: SCG DIAMOND database already present, skipping"
            fi
            touch {output.scg_db}

            if [ ! -f "{params.h_genome_path}" ]; then
                echo "initiate_dbs: downloading human reference genome..."
                wget --show-progress -P {params.human_dir} https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.fna.gz 2>&1 | tee -a {log}
            else
                echo "initiate_dbs: human reference genome already present, skipping"
            fi
            touch {output.h_genome_done}

            if [ ! -f "{params.markers_dir}/markers.fa" ]; then
                echo "initiate_dbs: downloading anthropogenic markers (pBI143, crAss001)..."
                wget --show-progress -O {params.markers_dir}/pBI143.fasta \
                    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=U30316.1&rettype=fasta&retmode=text" 2>&1 | tee -a {log}
                wget --show-progress -O {params.markers_dir}/crAss001.fasta \
                    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_049977.1&rettype=fasta&retmode=text" 2>&1 | tee -a {log}
                cat {params.markers_dir}/pBI143.fasta {params.markers_dir}/crAss001.fasta \
                    > {params.markers_dir}/markers.fa
            else
                echo "initiate_dbs: markers.fa already present, skipping download"
            fi
            if [ ! -f "{params.markers_dir}/markers_db.name" ]; then
                echo "initiate_dbs: indexing anthropogenic markers..."
                kma index -i {params.markers_dir}/markers.fa \
                          -o {params.markers_dir}/markers_db >> {log} 2>&1
            else
                echo "initiate_dbs: KMA markers index already present, skipping"
            fi
            touch {output.markers_db}

            if [ "{params.skip_metabolic}" = "True" ]; then
                echo "initiate_dbs: METABOLIC skipped (skip_metabolic=True)"
            else
                if [ ! -d "{params.metabolic_dir}" ]; then
                    echo "initiate_dbs: cloning METABOLIC..."
                    git clone --depth 1 https://github.com/AnantharamanLab/METABOLIC.git \
                        {params.metabolic_dir} >> {log} 2>&1
                    echo "initiate_dbs: METABOLIC ready"
                else
                    echo "initiate_dbs: METABOLIC already present, skipping"
                fi
            fi
            touch {output.metabolic_done}

        fi
        echo "initiate_dbs: done"
        """

# ----------------------------------------------------
#   read preprocessing, host filtering, and alignment
# ----------------------------------------------------

rule short_reads:
    input:
        r1          = get_r1,
        r2          = get_r2,
        afp_db      = os.path.join(_DB_AMR_KMA_DIR, f".done{_DB_TAG}"),
        scg_db      = os.path.join(_DB_SCG_DIR, f".done{_DB_TAG}"),
        h_genome    = os.path.join(_DB_HUMAN_DIR, f".done{_DB_TAG}"),
        markers_db  = os.path.join(_DB_MARKERS_DIR, f".done{_DB_TAG}")
    params:
        h_genome         = os.path.join(_DB_HUMAN_DIR, "GCF_000001405.40_GRCh38.p14_genomic.fna.gz"),
        min_read_length  = config.get("min_read_length", 50),
    output:
        json        = os.path.join(output_dir, "data", "QAQC", "fastp_reports", "{sample}.json"),
        r1_clean    = os.path.join(output_dir, "data", "clean_reads", "{sample}_R1.clean.fastq.gz"),
        r2_clean    = os.path.join(output_dir, "data", "clean_reads", "{sample}_R2.clean.fastq.gz"),
        nonpareil   = os.path.join(output_dir, "data", "QAQC", "nonpareil", "{sample}.npo"),
        afp         = os.path.join(output_dir, "data", "alignments", "{sample}.afp.res"),
        scgs        = os.path.join(output_dir, "data", "alignments", "{sample}.scgs"),
        markers     = os.path.join(output_dir, "data", "alignments", "{sample}.markers.res")
    threads: lambda wc: res(16, 4)
    resources:
        mem_mb  = lambda wc: res(20000, 4000),
        threads = lambda wc: res(16, 4),
    benchmark:
        os.path.join(output_dir, "data", "QAQC", "benchmarks", "{sample}.short_reads.txt")
    log:
        sample_log("short_reads")
    conda: "../envs/shortreads.yaml"
    shell:
        """
        set -euo pipefail
        mkdir -p "$(dirname {log})" \
                 "$(dirname {output.json})" \
                 "$(dirname {output.r1_clean})" \
                 "$(dirname {output.nonpareil})" \
                 "$(dirname {output.afp})"
        export TMPDIR="$(dirname {output.r1_clean})"
        
        TMP_R1="$TMPDIR/{wildcards.sample}_R1.fastq.gz"
        TMP_R2="$TMPDIR/{wildcards.sample}_R2.fastq.gz"
        
        echo "[{wildcards.sample}] fastp: quality filtering reads..."
        fastp -i {input.r1} -I {input.r2} -o "$TMP_R1" -O "$TMP_R2" --html /dev/null --json {output.json} --length_required {params.min_read_length} >> {log} 2>&1
        echo "[{wildcards.sample}] fastp: done"

        echo "[{wildcards.sample}] minimap2: filtering host reads..."
        minimap2 -t {resources.threads} -ax sr {params.h_genome} "$TMP_R1" "$TMP_R2" 2>> {log} \
        | samtools view -u -f 12 -F 256 - 2>> {log} \
        | samtools fastq -n -1 >(pigz -p {resources.threads} > {output.r1_clean}) \
                 -2 >(pigz -p {resources.threads} > {output.r2_clean}) \
                 -0 /dev/null - 2>> {log}
        rm -f "$TMP_R1" "$TMP_R2"
        echo "[{wildcards.sample}] minimap2: done"
        
        echo "[{wildcards.sample}] KMA: aligning to AMRFinderPlus..."
        kma -ipe {output.r1_clean} {output.r2_clean} -o $(dirname {output.afp})/{wildcards.sample}.afp -t_db $(dirname {input.afp_db})/afp_db -1t1 -t {resources.threads} >> {log} 2>&1
        echo "[{wildcards.sample}] KMA: done (AMR)"
        
        echo "[{wildcards.sample}] KMA: aligning to anthropogenic markers..."
        kma -ipe {output.r1_clean} {output.r2_clean} -o $(dirname {output.markers})/{wildcards.sample}.markers -t_db $(dirname {input.markers_db})/markers_db -1t1 -t {resources.threads} >> {log} 2>&1
        echo "[{wildcards.sample}] KMA: done (markers)"
        
        echo "[{wildcards.sample}] DIAMOND: aligning to single-copy genes..."
        TMP_SE="$TMPDIR/{wildcards.sample}_SE.fastq.gz"
        cat {output.r1_clean} {output.r2_clean} > "$TMP_SE"
        diamond blastx --db $(dirname {input.scg_db})/scg_db \
            --query "$TMP_SE" \
            --out {output.scgs} \
            --outfmt 6 qseqid sseqid pident length evalue bitscore slen \
            --fast --max-target-seqs 1 --threads {resources.threads} --quiet >> {log} 2>&1
        rm -f "$TMP_SE"
        echo "[{wildcards.sample}] DIAMOND: done"
        
        echo "[{wildcards.sample}] nonpareil: estimating metagenomic coverage..."
        base=$(basename {output.r1_clean} .gz)
        gunzip -c {output.r1_clean} > $(dirname {output.nonpareil})/$base
        nonpareil -s $(dirname {output.nonpareil})/$base -T kmer -f fastq -b $(dirname {output.nonpareil})/{wildcards.sample} >> {log} 2>&1
        rm $(dirname {output.nonpareil})/$base
        echo "[{wildcards.sample}] nonpareil: done"
        
        """
         
# ---------------------------------------------------------------------------
#   Per-sample: ARGNet on host-filtered short reads in discovery mode
# ---------------------------------------------------------------------------

rule short_reads_argnet:
    input:
        r1          = get_clean_r1,
        r2          = get_clean_r2,
        argnet_done = os.path.join(_DB_ARGNET_DIR, ".done")
    output:
        tsv = os.path.join(output_dir, "data", "argnet_short_reads", "{sample}_argnet.tsv")
    params:
        repo_dir = _ARGNET_REPO_DIR
    threads: lambda wc: res(8, 4)
    resources:
        mem_mb  = lambda wc: res(32000, 8000),
        threads = lambda wc: res(8, 4),
    log:
        sample_log("short_reads_argnet")
    conda: "../envs/argnet.yaml"
    shell:
        """
        set -euo pipefail
        mkdir -p "$(dirname {log})" "$(dirname {output.tsv})"
        TMP_FASTA="$(mktemp "$(dirname {output.tsv})/{wildcards.sample}.short_reads_argnet.XXXXXX.fasta")"

        gunzip -c {input.r1} | awk '
            NR % 4 == 1 {{
                hdr = $0
                sub(/^@/, "", hdr)
                print ">R1|" hdr
                next
            }}
            NR % 4 == 2 {{ print }}
        ' > "$TMP_FASTA"

        gunzip -c {input.r2} | awk '
            NR % 4 == 1 {{
                hdr = $0
                sub(/^@/, "", hdr)
                print ">R2|" hdr
                next
            }}
            NR % 4 == 2 {{ print }}
        ' >> "$TMP_FASTA"

        echo "[{wildcards.sample}] ARGNet: annotating host-filtered short reads..."
        python workflow/scripts/run_argnet.py \
            --repo {params.repo_dir} \
            --input "$TMP_FASTA" \
            --output {output.tsv} \
            --seq-type nt \
            --model argnet-s >> {log} 2>&1
        rm -f "$TMP_FASTA"
        echo "[{wildcards.sample}] ARGNet: done"
        """

# ------------------------------------------
#           Summarize results 
# ------------------------------------------        

rule short_reads_summary:
    input:
        afp_files       = expand(os.path.join(output_dir, "data", "alignments", "{sample}.afp.res"), sample=samples),
        scgs            = expand(os.path.join(output_dir, "data", "alignments", "{sample}.scgs"), sample=samples),
        npos            = expand(os.path.join(output_dir, "data", "QAQC", "nonpareil", "{sample}.npo"), sample=samples),
        jsons           = expand(os.path.join(output_dir, "data", "QAQC", "fastp_reports", "{sample}.json"), sample=samples),
        afp_metadata    = os.path.join(_DB_AMR_KMA_DIR, "ReferenceGeneCatalog.txt"),
        markers_files   = expand(os.path.join(output_dir, "data", "alignments", "{sample}.markers.res"), sample=samples)
    output:
        fastp            = FASTP_SUMMARY,
        afp              = SHORT_READS_OUTPUT,
        markers_cpg      = MARKERS_CPG_OUTPUT
    conda: "../envs/Renv.yaml"
    log:
        rule_log("short_reads_summary")
    script:
        "../scripts/short_reads_processing.R"
