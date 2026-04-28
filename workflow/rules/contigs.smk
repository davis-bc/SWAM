# ------------------------------------------
#       prepare geNomad database 
# ------------------------------------------

CONTIG_MEM_MB_BY_ATTEMPT = (150000, 180000, 220000)
CONTIG_RUNTIME_BY_ATTEMPT = (4320, 5760, 7200)
GENOMAD_MEM_MB_BY_ATTEMPT = (150000, 180000, 220000)
GENOMAD_RUNTIME_BY_ATTEMPT = (4320, 5760, 7200)


def contigs_mem_mb(_, attempt):
    return res(attempt_value(CONTIG_MEM_MB_BY_ATTEMPT, attempt), 8000)


def contigs_runtime(_, attempt):
    return 240 if _TEST else attempt_value(CONTIG_RUNTIME_BY_ATTEMPT, attempt)


def genomad_mem_mb(_, attempt):
    return res(attempt_value(GENOMAD_MEM_MB_BY_ATTEMPT, attempt), 8000)


def genomad_runtime(_, attempt):
    return 240 if _TEST else attempt_value(GENOMAD_RUNTIME_BY_ATTEMPT, attempt)

rule init_genomad:
    output:
        genomad_db   = os.path.join(_DB_GENOMAD_DIR, ".done")
    params:
        db_dir = _DB_GENOMAD_DIR,
        legacy_dir = _LEGACY_GENOMAD_DB_DIR
    log:
        rule_log("init_genomad")
    conda: "../envs/genomad.yaml"
    shell:
        """
        mkdir -p "$(dirname {log})" "{params.db_dir}"
        if [ -d "{params.legacy_dir}" ] && [ ! -d "{params.db_dir}/genomad_db" ]; then
            mv -n "{params.legacy_dir}"/* "{params.db_dir}/" 2>/dev/null || true
        fi
        if [ ! -d "{params.db_dir}/genomad_db" ]; then
        
        echo "geNomad: downloading database..."
        cd {params.db_dir}
        curl -L --progress-bar https://zenodo.org/records/14886553/files/genomad_db_v1.9.tar.gz 2> >(tee -a {log} >&2) | tar -xz 2>> {log}
        touch {output.genomad_db}
        echo "geNomad: database ready"
        
        else 
        
        echo "geNomad: database already exists, skipping download"
        touch {output.genomad_db}
        
        fi
        """
# ------------------------------------------
#       de novo assembly + gene annotation 
# ------------------------------------------

rule contigs:
    input:
        r1       = get_clean_r1,
        r2       = get_clean_r2
    output:
        megahit  = os.path.join(output_dir, "data", "megahit", "{sample}.contigs.fa")
    threads: lambda wc: res(32, 4)
    resources:
        mem_mb  = contigs_mem_mb,
        runtime = contigs_runtime,
        threads = lambda wc: res(32, 4),
    benchmark:
        os.path.join(output_dir, "data", "QAQC", "benchmarks", "{sample}.contigs.txt")
    log:
        sample_log("contigs")
    conda: "../envs/contigs.yaml"
    shell:
        """
        set -euo pipefail
        mkdir -p "$(dirname {log})" "$(dirname {output.megahit})"
        echo "[{wildcards.sample}] MEGAHIT: assembling contigs..."

        # create temporary directory in same parent so moves are on same fs (atomic)
        tmpdir="$(mktemp -d "$(dirname {output.megahit})/megahit_{wildcards.sample}.tmp.XXXXXX")"

        # Run MEGAHIT into the temporary directory
        megahit -1 {input.r1} -2 {input.r2} \
                -t {resources.threads} \
                --presets meta-large \
                --min-contig-len 1000 \
                -o "$tmpdir" \
                --out-prefix {wildcards.sample} \
                --continue \
                -f >> {log} 2>&1

        # Move the contigs file into the final megahit directory.
        mv -f "$tmpdir/{wildcards.sample}.contigs.fa" {output.megahit}

        # Clean up the rest of temporary files
        rm -rf "$tmpdir"
        
        # rename fasta headers to include sample information
        sed -E -i 's/>/>{wildcards.sample}-/g' {output.megahit}
        echo "[{wildcards.sample}] MEGAHIT: done"
        
        """
    
# ------------------------------------------------------------------------
#       classify contigs into chromosome, plasmid, and phage with geNomad
# ------------------------------------------------------------------------

rule genomad:
    input:
        contigs      = os.path.join(output_dir, "data", "megahit", "{sample}.contigs.fa"),
        genomad_db   = os.path.join(_DB_GENOMAD_DIR, ".done")
    output:
        genomad      = directory(os.path.join(output_dir, "data", "genomad", "{sample}")),
        plas_contigs = os.path.join(output_dir, "data", "genomad", "{sample}", "{sample}.contigs_summary", "{sample}.contigs_plasmid.fna"),
        plas_summary = os.path.join(output_dir, "data", "genomad", "{sample}", "{sample}.contigs_summary", "{sample}.contigs_plasmid_summary.tsv"),
        vir_summary  = os.path.join(output_dir, "data", "genomad", "{sample}", "{sample}.contigs_summary", "{sample}.contigs_virus_summary.tsv")
    params:
        splits = config.get("genomad_splits", 1),
        db_dir = os.path.join(_DB_GENOMAD_DIR, "genomad_db")
    threads: lambda wc: res(32, 4)
    resources:
        mem_mb  = genomad_mem_mb,
        runtime = genomad_runtime,
        threads = lambda wc: res(32, 4),
    benchmark:
        os.path.join(output_dir, "data", "QAQC", "benchmarks", "{sample}.genomad.txt")
    log:
        sample_log("genomad")
    conda: "../envs/genomad.yaml"
    shell:
        """
        mkdir -p "$(dirname {log})"
        echo "[{wildcards.sample}] geNomad: classifying contigs (chromosome/plasmid/phage)..."
        genomad end-to-end {input.contigs} {output.genomad} {params.db_dir} \
            --relaxed --cleanup --splits {params.splits} >> {log} 2>&1
        # touch missing output files (geNomad omits them when no hits are found)
        touch {output.plas_contigs} {output.plas_summary} {output.vir_summary}
        echo "[{wildcards.sample}] geNomad: done"
        """

# ----------------------------------------------------------------------------------
#       determine cirularity of plasmid contigs, infer plasmid systems with MobMess 
# ----------------------------------------------------------------------------------

rule mobmess:
    input:
        r1           = get_clean_r1,
        r2           = get_clean_r2,
        plas_contigs = os.path.join(output_dir, "data", "genomad", "{sample}",  "{sample}.contigs_summary", "{sample}.contigs_plasmid.fna")
    output:
        circular_txt = os.path.join(output_dir, "data", "mobmess", "{sample}.contigs_circular.txt"),
        contig_bam   = os.path.join(output_dir, "data", "mobmess", "{sample}.contigs_plasmid.bam"),
        mobmess      = os.path.join(output_dir, "data", "mobmess", "{sample}-mobmess_contigs.txt")
    conda: "../envs/contigs.yaml"
    threads: lambda wc: res(32, 4)
    resources:
        mem_mb  = lambda wc: res(150000, 8000),
        threads = lambda wc: res(32, 4),
    log:
        sample_log("mobmess")
    shell:
        """
        mkdir -p "$(dirname {log})"
        mkdir -p $(dirname {output.contig_bam})
        # Skip processing if plasmid FASTA is empty (no plasmids found by geNomad)
        if [ ! -s {input.plas_contigs} ]; then
            echo "[{wildcards.sample}] MobMess: no plasmids found, skipping"
            touch {output.contig_bam} {output.circular_txt} {output.mobmess}
        else
            echo "[{wildcards.sample}] minimap2: mapping reads to plasmid contigs..."
            minimap2 -t {resources.threads} -ax sr {input.plas_contigs} {input.r1} {input.r2} 2>> {log} \
            | samtools sort -@4 -o {output.contig_bam} - 2>> {log} && samtools index {output.contig_bam} 2>> {log}
            echo "[{wildcards.sample}] minimap2: done"
            
            echo "[{wildcards.sample}] MobMess: detecting circular plasmid contigs..."
            workflow/scripts/detect_circular_contigs.py -b {output.contig_bam} -o {output.circular_txt} >> {log} 2>&1
            
            echo "[{wildcards.sample}] MobMess: inferring plasmid systems..."
            mobmess systems --sequences {input.plas_contigs} --complete {output.circular_txt} --output $(dirname {output.mobmess})/{wildcards.sample}-mobmess --threads {resources.threads} >> {log} 2>&1 || true
            # ensure output file always exists
            [ -f {output.mobmess} ] || touch {output.mobmess}
            echo "[{wildcards.sample}] MobMess: done"
        fi
        """

# ---------------------------------------------------------------------------
#   Build / validate MMseqs2 UniRef50 taxonomy database
#   Test mode: build from mini FASTA in test/dbs/uniref50/
#   Production: validate user-provided db path
# ---------------------------------------------------------------------------

rule init_mmseqs_db:
    output:
        done = os.path.join(_DB_UNIREF50_DIR, ".done")
    params:
        test_mode    = _TEST,
        db_prefix    = _UNIREF50_DB,
        dbs_dir      = _DB_UNIREF50_DIR,
        test_fasta   = os.path.join(_REPO, "test", "dbs", "uniref50", "uniref50_mini.fasta"),
        test_tax_dir = os.path.join(_REPO, "test", "dbs", "uniref50"),
        test_mapping = os.path.join(_REPO, "test", "dbs", "uniref50", "tax_mapping.tsv"),
    threads: lambda wc: res(16, 4)
    resources:
        mem_mb  = lambda wc: res(32000, 4000),
        threads = lambda wc: res(16, 4),
    conda: "../envs/contigs.yaml"
    log:
        rule_log("init_mmseqs_db")
    shell:
        """
        mkdir -p "$(dirname {log})"
        mkdir -p "$(dirname {params.db_prefix})"

        if [ "{params.test_mode}" = "True" ]; then
            echo "MMseqs2: building test taxonomy database..."
            mkdir -p {params.test_tax_dir}/tmp
            touch {params.test_tax_dir}/merged.dmp
            mmseqs createdb {params.test_fasta} {params.db_prefix} >> {log} 2>&1
            mmseqs createtaxdb {params.db_prefix} {params.test_tax_dir}/tmp \
                --ncbi-tax-dump {params.test_tax_dir} \
                --tax-mapping-file {params.test_mapping} >> {log} 2>&1
            mmseqs createindex {params.db_prefix} {params.test_tax_dir}/tmp --search-type 2 >> {log} 2>&1
            rm -rf {params.test_tax_dir}/tmp
            echo "MMseqs2: test database ready"
        else
            if [ ! -f "{params.db_prefix}.dbtype" ]; then
                UNIREF_FASTA="{params.dbs_dir}/uniref50.fasta"
                TAX_DIR="{params.dbs_dir}/taxdump"
                TMP_DIR="{params.dbs_dir}/tmp"
                mkdir -p "$TAX_DIR" "$TMP_DIR"

                echo "MMseqs2: downloading UniRef50 FASTA (~9 GB)..."
                wget --show-progress -O "$UNIREF_FASTA.gz" \
                    https://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref50/uniref50.fasta.gz 2>&1 | tee -a {log}
                echo "MMseqs2: decompressing UniRef50..."
                gunzip -f "$UNIREF_FASTA.gz" >> {log} 2>&1

                echo "MMseqs2: downloading NCBI taxonomy..."
                wget --show-progress -P "$TAX_DIR" \
                    https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz 2>&1 | tee -a {log}
                tar -xzf "$TAX_DIR/taxdump.tar.gz" -C "$TAX_DIR" >> {log} 2>&1

                echo "MMseqs2: building taxonomy database (this may take several hours)..."
                mmseqs createdb "$UNIREF_FASTA" {params.db_prefix} 2>&1 | tee -a {log}
                mmseqs createtaxdb {params.db_prefix} "$TMP_DIR" \
                    --ncbi-tax-dump "$TAX_DIR" \
                    --threads {resources.threads} 2>&1 | tee -a {log}
                mmseqs createindex {params.db_prefix} "$TMP_DIR" \
                    --search-type 2 \
                    --threads {resources.threads} 2>&1 | tee -a {log}

                rm -f "$UNIREF_FASTA"
                rm -rf "$TMP_DIR" "$TAX_DIR"
                echo "MMseqs2: production database ready"
            else
                echo "MMseqs2: database already exists at {params.db_prefix}, skipping"
            fi
        fi
        touch {output.done}
        """

# ---------------------------------------------------------------------------
#   Per-sample: gene calling with Prodigal (shared output for AMR + summary)
# ---------------------------------------------------------------------------

rule prodigal:
    input:
        contigs = os.path.join(output_dir, "data", "megahit", "{sample}.contigs.fa")
    output:
        faa = os.path.join(output_dir, "data", "prodigal", "{sample}.faa"),
        gff = os.path.join(output_dir, "data", "prodigal", "{sample}.gff")
    resources:
        mem_mb  = lambda wc: res(16000, 4000),
        threads = 1,
    log:
        sample_log("prodigal")
    conda: "../envs/contigs.yaml"
    shell:
        """
        mkdir -p "$(dirname {log})"
        mkdir -p $(dirname {output.faa})
        echo "[{wildcards.sample}] Prodigal: predicting genes..."
        prodigal -i {input.contigs} -a {output.faa} -f gff -o {output.gff} -p meta -q >> {log} 2>&1
        echo "[{wildcards.sample}] Prodigal: done"
        """

# ---------------------------------------------------------------------------
#   One-time: download AMRFinderPlus database (skipped in test mode)
# ---------------------------------------------------------------------------

rule init_amrfinder_db:
    """Download the AMRFinderPlus database to dbs/amrfinderplus/. Skipped in test mode."""
    output:
        done = os.path.join(_AFP_DB_DIR, ".done")
    params:
        db_dir      = _AFP_DB_DIR,
        legacy_dir  = _LEGACY_AFP_DB_DIR,
        test_mode   = _TEST
    log:
        rule_log("init_amrfinder_db")
    conda: "../envs/mags.yaml"
    shell:
        """
        mkdir -p "$(dirname {log})"
        mkdir -p {params.db_dir}
        if [ -d "{params.legacy_dir}" ] && [ ! -f "{params.db_dir}/AMR_CDS.fa" ]; then
            mv -n "{params.legacy_dir}"/* "{params.db_dir}/" 2>/dev/null || true
        fi
        if [ "{params.test_mode}" = "True" ]; then
            echo "init_amrfinder_db: test mode — skipping database download"
        elif [ -f "{params.db_dir}/AMR_CDS.fa" ]; then
            echo "init_amrfinder_db: database already exists, skipping download"
        else
            echo "init_amrfinder_db: downloading AMRFinderPlus database (~600 MB)..."
            amrfinder --update 2>&1 | tee -a {log}
            # amrfinder --update downloads to $CONDA_PREFIX/share/amrfinderplus/data/latest/
            AFP_SRC="$CONDA_PREFIX/share/amrfinderplus/data/latest"
            if [ -d "$AFP_SRC" ]; then
                rsync -a --info=progress2 "$AFP_SRC/" {params.db_dir}/ 2>&1 | tee -a {log}
                echo "init_amrfinder_db: database ready"
            else
                echo "ERROR: expected database at $AFP_SRC not found" >> {log}
                exit 1
            fi
        fi
        touch {output.done}
        """

# ---------------------------------------------------------------------------
#   One-time: clone ARGNet and fetch long-sequence models for discovery mode
# ---------------------------------------------------------------------------

rule init_argnet:
    output:
        done = os.path.join(_DB_ARGNET_DIR, ".done")
    params:
        db_dir   = _DB_ARGNET_DIR,
        repo_dir = _ARGNET_REPO_DIR
    log:
        rule_log("init_argnet")
    shell:
        """
        set -euo pipefail
        mkdir -p "$(dirname {log})" "{params.db_dir}"
        if [ ! -d "{params.repo_dir}/.git" ]; then
            echo "ARGNet: cloning repository..." | tee -a {log}
            if ! git clone --depth 1 https://github.com/id-bioinfo/ARGnet.git "{params.repo_dir}" >> {log} 2>&1; then
                git clone --depth 1 https://github.com/patience111/ARGNet.git "{params.repo_dir}" >> {log} 2>&1
            fi
        else
            echo "ARGNet: repository already present, skipping clone" | tee -a {log}
        fi
        mkdir -p "{params.repo_dir}/model" "{params.repo_dir}/results"
        if [ ! -f "{params.repo_dir}/model/AELS_tall.h5" ]; then
            echo "ARGNet: downloading AELS_tall.h5..." | tee -a {log}
            curl -L --progress-bar \
                -o "{params.repo_dir}/model/AELS_tall.h5" \
                https://github.com/patience111/ARGNet/releases/download/v1.0.1/AELS_tall.h5 2>> {log}
        fi
        if [ ! -f "{params.repo_dir}/model/classifier-ls_tall.h5" ]; then
            echo "ARGNet: downloading classifier-ls_tall.h5..." | tee -a {log}
            curl -L --progress-bar \
                -o "{params.repo_dir}/model/classifier-ls_tall.h5" \
                https://github.com/patience111/ARGNet/releases/download/v1.0.1/classifier-ls_tall.h5 2>> {log}
        fi
        if [ ! -f "{params.repo_dir}/model/AESS_tall.h5" ]; then
            echo "ARGNet: downloading AESS_tall.h5..." | tee -a {log}
            curl -L --progress-bar \
                -o "{params.repo_dir}/model/AESS_tall.h5" \
                https://github.com/patience111/ARGNet/releases/download/v1.0.1/AESS_tall.h5 2>> {log}
        fi
        if [ ! -f "{params.repo_dir}/model/classifier-ss_tall.h5" ]; then
            echo "ARGNet: downloading classifier-ss_tall.h5..." | tee -a {log}
            curl -L --progress-bar \
                -o "{params.repo_dir}/model/classifier-ss_tall.h5" \
                https://github.com/patience111/ARGNet/releases/download/v1.0.1/classifier-ss_tall.h5 2>> {log}
        fi
        touch {output.done}
        """

# ---------------------------------------------------------------------------
#   Per-sample: AMR annotation on contigs (protein + nucleotide + GFF mode)
# ---------------------------------------------------------------------------

rule contig_amr:
    input:
        contigs     = os.path.join(output_dir, "data", "megahit", "{sample}.contigs.fa"),
        faa         = os.path.join(output_dir, "data", "prodigal", "{sample}.faa"),
        gff         = os.path.join(output_dir, "data", "prodigal", "{sample}.gff"),
        afp_db_done = os.path.join(_AFP_DB_DIR, ".done")
    output:
        tsv = os.path.join(output_dir, "data", "amr_contigs", "{sample}_contig_amr.tsv")
    params:
        afp_db_dir = _AFP_DB_DIR
    threads: lambda wc: res(16, 4)
    resources:
        mem_mb  = lambda wc: res(16000, 4000),
        threads = lambda wc: res(16, 4),
    log:
        sample_log("contig_amr")
    conda: "../envs/mags.yaml"
    shell:
        """
        mkdir -p "$(dirname {log})"
        mkdir -p $(dirname {output.tsv})
        echo "[{wildcards.sample}] AMRFinderPlus: annotating contigs..."
        amrfinder \
            -n {input.contigs} \
            -p {input.faa} \
            -g {input.gff} \
            --database {params.afp_db_dir} \
            --threads {resources.threads} \
            --plus \
            -o {output.tsv} \
            --annotation_format prodigal >> {log} 2>&1 || true
        # Create empty file if no hits found or database not yet set up
        [ -f {output.tsv} ] || touch {output.tsv}
        echo "[{wildcards.sample}] AMRFinderPlus: done"
        """

# ---------------------------------------------------------------------------
#   Per-sample: ARGNet on Prodigal proteins for discovery mode
# ---------------------------------------------------------------------------

rule contig_argnet:
    input:
        faa         = os.path.join(output_dir, "data", "prodigal", "{sample}.faa"),
        argnet_done = os.path.join(_DB_ARGNET_DIR, ".done")
    output:
        tsv = os.path.join(output_dir, "data", "argnet_contigs", "{sample}_argnet.tsv")
    params:
        repo_dir = _ARGNET_REPO_DIR
    threads: lambda wc: res(8, 4)
    resources:
        mem_mb  = lambda wc: res(32000, 8000),
        threads = lambda wc: res(8, 4),
    log:
        sample_log("contig_argnet")
    conda: "../envs/argnet.yaml"
    shell:
        """
        set -euo pipefail
        mkdir -p "$(dirname {log})" "$(dirname {output.tsv})"
        echo "[{wildcards.sample}] ARGNet: annotating contig proteins..."
        python workflow/scripts/run_argnet.py \
            --repo {params.repo_dir} \
            --input {input.faa} \
            --output {output.tsv} \
            --seq-type aa \
            --model argnet-l >> {log} 2>&1
        echo "[{wildcards.sample}] ARGNet: done"
        """

# ---------------------------------------------------------------------------
#   Per-sample: MGE annotation on contigs with MobileElementFinder
# ---------------------------------------------------------------------------

rule init_mge_tool:
    """Install MobileElementFinder via pip --no-deps (conda provides all deps)."""
    output:
        done = os.path.join(output_dir, "data", "mge_contigs", ".mef_install.done")
    log:
        rule_log("init_mge_tool")
    conda: "../envs/contigs.yaml"
    shell:
        """
        mkdir -p "$(dirname {log})"
        mkdir -p $(dirname {output.done})
        echo "MobileElementFinder: installing..."
        pip install --no-deps --quiet mypy-extensions >> {log} 2>&1
        pip install --no-deps --quiet MGEdb==1.1.1 >> {log} 2>&1
        pip install --no-deps --quiet MobileElementFinder >> {log} 2>&1
        touch {output.done}
        echo "MobileElementFinder: installed"
        """

rule mge_annotation:
    input:
        contigs  = os.path.join(output_dir, "data", "megahit", "{sample}.contigs.fa"),
        mef_done = os.path.join(output_dir, "data", "mge_contigs", ".mef_install.done")
    output:
        tsv = os.path.join(output_dir, "data", "mge_contigs", "{sample}_mge.tsv"),
        tmp = temp(directory(os.path.join(output_dir, "data", "mge_contigs", "{sample}_tmp")))
    threads: lambda wc: res(16, 4)
    resources:
        mem_mb  = lambda wc: res(16000, 4000),
        threads = lambda wc: res(16, 4),
    log:
        sample_log("mge_annotation")
    conda: "../envs/contigs.yaml"
    shell:
        """
        mkdir -p "$(dirname {log})"
        mkdir -p {output.tmp}
        CHUNK_DIR={output.tmp}/chunks
        mkdir -p "$CHUNK_DIR"

        # Split assembly into chunks of ≤200 contigs.
        # Workaround for MEF ≥1.1 JSONDecodeError: blastn -outfmt 15 writes
        # multiple concatenated JSON objects for large inputs; json.load() fails.
        echo "[{wildcards.sample}] MobileElementFinder: splitting contigs into chunks..."
        awk -v CHUNK="$CHUNK_DIR" \
            'BEGIN{{n=0;c=0}} \
             /^>/{{c++;if(c==1||c%200==1){{n++;close(f);f=CHUNK"/"n".fa"}} print > f;next}} \
             {{print > f}}' \
            {input.contigs}

        TMP_CSV={output.tmp}/combined.csv
        HEADER_WRITTEN=0
        echo "[{wildcards.sample}] MobileElementFinder: annotating MGEs..."
        for chunk in "$CHUNK_DIR"/*.fa; do
            [ -f "$chunk" ] || continue
            base=$(basename "$chunk" .fa)
            CHUNK_PREFIX={output.tmp}/$base
            CHUNK_TMP={output.tmp}/tmp_$base
            mkdir -p "$CHUNK_TMP"
            mefinder find \
                --contig "$chunk" \
                --temp-dir "$CHUNK_TMP" \
                --threads {resources.threads} \
                "$CHUNK_PREFIX" >> {log} 2>&1 || true
            csv="${{CHUNK_PREFIX}}.csv"
            if [ -s "$csv" ]; then
                if [ $HEADER_WRITTEN -eq 0 ]; then
                    cat "$csv" > "$TMP_CSV"
                    HEADER_WRITTEN=1
                else
                    # Append only data rows (skip comment + header lines)
                    grep -v "^#" "$csv" | tail -n +2 >> "$TMP_CSV"
                fi
            fi
        done

        if [ -f "$TMP_CSV" ] && [ -s "$TMP_CSV" ]; then
            mv "$TMP_CSV" {output.tsv}
        else
            touch {output.tsv}
        fi
        echo "[{wildcards.sample}] MobileElementFinder: done"
        """

# ---------------------------------------------------------------------------
#   Per-sample: MMseqs2 taxonomy assignment on assembled contigs
# ---------------------------------------------------------------------------

rule mmseqs_taxonomy:
    input:
        contigs  = os.path.join(output_dir, "data", "megahit", "{sample}.contigs.fa"),
        db_done  = os.path.join(_DB_UNIREF50_DIR, ".done")
    output:
        lca = os.path.join(output_dir, "data", "mmseqs", "{sample}_lca.tsv")
    params:
        db_prefix = _UNIREF50_DB,
        tmp       = lambda wc: os.path.join(output_dir, "data", "mmseqs", f"{wc.sample}_tmp")
    threads: lambda wc: res(32, 4)
    resources:
        mem_mb  = lambda wc: res(64000, 4000),
        threads = lambda wc: res(32, 4),
    log:
        sample_log("mmseqs_taxonomy")
    conda: "../envs/contigs.yaml"
    shell:
        """
        mkdir -p "$(dirname {log})"
        mkdir -p $(dirname {output.lca}) {params.tmp}
        echo "[{wildcards.sample}] MMseqs2: taxonomy classification..."
        mmseqs easy-taxonomy \
            {input.contigs} \
            {params.db_prefix} \
            $(dirname {output.lca})/{wildcards.sample} \
            {params.tmp} \
            --threads {resources.threads} \
            --lca-ranks superkingdom,phylum,class,order,family,genus,species \
            --format-output "query,taxid,taxname,taxlineage" \
            --orf-filter 0 \
            -s 2 \
            >> {log} 2>&1 || true
        # easy-taxonomy writes {{prefix}}_lca.tsv; touch if absent (no hits or search failed)
        mv $(dirname {output.lca})/{wildcards.sample}_lca.tsv {output.lca} 2>/dev/null || true
        rm -rf {params.tmp}
        [ -f {output.lca} ] || touch {output.lca}
        echo "[{wildcards.sample}] MMseqs2: done"
        """

# ---------------------------------------------------------------------------
#   Per-sample: map clean reads to all contigs; compute abundance with CoverM
# ---------------------------------------------------------------------------

rule contig_abundance:
    input:
        r1      = get_clean_r1,
        r2      = get_clean_r2,
        contigs = get_contigs,
    output:
        bam = os.path.join(output_dir, "data", "contig_abundance", "{sample}_contigs.bam"),
        bai = os.path.join(output_dir, "data", "contig_abundance", "{sample}_contigs.bam.bai"),
        tsv = os.path.join(output_dir, "data", "contig_abundance", "{sample}_contig_abundance.tsv")
    threads: lambda wc: res(16, 4)
    resources:
        mem_mb  = lambda wc: res(20000, 4000),
        threads = lambda wc: res(16, 4),
    log:
        sample_log("contig_abundance")
    conda: "../envs/shortreads.yaml"
    shell:
        """
        set -euo pipefail
        mkdir -p "$(dirname {log})"
        mkdir -p $(dirname {output.bam})

        echo "[{wildcards.sample}] minimap2: mapping reads to contigs..."
        minimap2 -t {resources.threads} -ax sr {input.contigs} {input.r1} {input.r2} 2>> {log} \
        | samtools sort -@ 4 -o {output.bam} - >> {log} 2>&1
        samtools index {output.bam} >> {log} 2>&1
        echo "[{wildcards.sample}] minimap2: done"

        echo "[{wildcards.sample}] CoverM: computing per-contig abundance..."
        coverm contig \
            --bam-files {output.bam} \
            --methods trimmed_mean \
            --min-covered-fraction 0 \
            --threads {resources.threads} 2>> {log} \
        | awk 'NR==1{{print "contig_id\ttrimmed_mean"}} NR>1{{print}}' \
        > {output.tsv}
        echo "[{wildcards.sample}] CoverM: done"
        """

# ---------------------------------------------------------------------------
#   Aggregation: join all per-sample annotations into one summary TSV
# ---------------------------------------------------------------------------

rule contig_summary:
    input:
        plas_summaries = expand(os.path.join(output_dir, "data", "genomad", "{sample}", "{sample}.contigs_summary", "{sample}.contigs_plasmid_summary.tsv"), sample=samples),
        vir_summaries  = expand(os.path.join(output_dir, "data", "genomad", "{sample}", "{sample}.contigs_summary", "{sample}.contigs_virus_summary.tsv"),  sample=samples),
        lca_files      = expand(os.path.join(output_dir, "data", "mmseqs", "{sample}_lca.tsv"), sample=samples),
        abund_files    = expand(os.path.join(output_dir, "data", "contig_abundance", "{sample}_contig_abundance.tsv"), sample=samples),
        amr_files      = expand(os.path.join(output_dir, "data", "amr_contigs", "{sample}_contig_amr.tsv"), sample=samples),
        mge_files      = expand(os.path.join(output_dir, "data", "mge_contigs", "{sample}_mge.tsv"), sample=samples),
        prodigal_gffs  = expand(os.path.join(output_dir, "data", "prodigal", "{sample}.gff"), sample=samples)
    output:
        tsv = CONTIG_SUMMARY_OUTPUT
    params:
        samples    = samples,
        output_dir = output_dir
    conda: "../envs/contigs.yaml"
    log:
        rule_log("contig_summary")
    script:
        "../scripts/contig_summary.py"


# ---------------------------------------------------------------------------
#   Aggregation: per-sample assembly QA metrics
# ---------------------------------------------------------------------------

rule assembly_qa:
    """Per-sample assembly QA: N50, total length, read mapping stats."""
    input:
        bam_files    = expand(os.path.join(output_dir, "data", "contig_abundance", "{sample}_contigs.bam"), sample=samples),
        megahit_logs = expand(sample_log("contigs"), sample=samples)
    output:
        tsv = ASSEMBLY_QA_OUTPUT
    params:
        samples    = samples,
        output_dir = output_dir
    conda: "../envs/contigs.yaml"
    log:
        rule_log("assembly_qa")
    script:
        "../scripts/assembly_qa.py"
        
        
        
        
        
        
