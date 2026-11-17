
BASE_DIR="/home/metatranscriptomics"

# ==========================================
# PART 1: Sample Processing Loop (QC, Trimming, Assembly, Prediction, Annotation)
# ==========================================
for i in SRR37808534 SRR37808535 SRR37808538 SRR37808539
do
    echo "=== Processing sample: ${i} ==="

    # 1. Download the raw data
    mkdir -p "${BASE_DIR}/${i}"
    /home/sratoolkit.3.2.1-alma_linux64/bin/fastq-dump --split-files --gzip --outdir "${BASE_DIR}/${i}" "${i}"

    # 2. Quality control of the raw data (FastQC - optional/commented)
    # /home/FastQC/fastqc "${BASE_DIR}/${i}/${i}_1.fastq.gz" "${BASE_DIR}/${i}/${i}_2.fastq.gz" 
     
    # 3. Running Trimmomatic (using Nextera and polyG adapters)
    java -jar /home/BIB/Software/Trimmomatic-0.39/trimmomatic-0.39.jar \
    PE \
    -threads 4 \
    "${BASE_DIR}/${i}/${i}_1.fastq.gz" "${BASE_DIR}/${i}/${i}_2.fastq.gz" \
    "${BASE_DIR}/${i}/${i}_1p.fastq.gz" "${BASE_DIR}/${i}/${i}_1u.fastq.gz" \
    "${BASE_DIR}/${i}/${i}_2p.fastq.gz" "${BASE_DIR}/${i}/${i}_2u.fastq.gz" \
    ILLUMINACLIP:"${BASE_DIR}/nextera_polyG.fa":2:30:10 \
    LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:50

    # 4. SortMeRNA - Removal of ribosomal rRNA
    nohup sortmerna \
      --ref "${BASE_DIR}/sortmerna_db/smr_v4.3_default_db.fasta" \
      --reads "${BASE_DIR}/${i}/${i}_1p.fastq.gz" \
      --reads "${BASE_DIR}/${i}/${i}_2p.fastq.gz" \
      --workdir "${BASE_DIR}/${i}_sortmerna" \
      --fastx \
      --paired_in \
      --aligned "${BASE_DIR}/${i}_sortmerna/out/rRNA" \
      --other "${BASE_DIR}/${i}_sortmerna/out/non_rRNA" \
      > "${BASE_DIR}/${i}_sortmerna.log" 2>&1 &
      
    # 5. De novo assembly with MEGAHIT
    megahit \
      --12 "${BASE_DIR}/${i}_sortmerna/out/non_rRNA.fq.gz" \
      -o "${BASE_DIR}/${i}_megahit" \
      > "${BASE_DIR}/${i}_megahit.log" 2>&1 &
      
    # 6. Gene prediction by Prodigal
    mkdir -p "${BASE_DIR}/${i}_prodigal"
    nohup prodigal \
      -i "${BASE_DIR}/${i}_megahit/final.contigs.fa" \
      -a "${BASE_DIR}/${i}_prodigal/${i}_prodigal_proteins.faa" \
      -d "${BASE_DIR}/${i}_prodigal/${i}_prodigal_cds.fna" \
      -o "${BASE_DIR}/${i}_prodigal/${i}_prodigal_genes.gff" \
      -p meta 
     
    # 7. DIAMOND functional annotation
    mkdir -p "${BASE_DIR}/${i}_prodigal/diamond"
    diamond blastp \
      --query "${BASE_DIR}/${i}_prodigal/${i}_prodigal_proteins.faa" \
      --db /home/diamond/uniprot_sprot.dmnd \
      --out "${BASE_DIR}/${i}_prodigal/diamond/${i}_diamond.tsv" \
      --outfmt 6 \
      --evalue 1e-5 \
      --max-target-seqs 1 \
      --threads 8 \
      > "${BASE_DIR}/${i}_prodigal/diamond/${i}_diamond.log" 2>&1 &

    # 8. Renaming of the gene.fna from Prodigal (ensures unique gene IDs across samples)
    mkdir -p "${BASE_DIR}/renamed_genes"
    awk -v sample="${i}" '/^>/{sub(/^>/,">" sample "__");}1' \
      "${BASE_DIR}/${i}_prodigal/${i}_prodigal_cds.fna" \
    > "${BASE_DIR}/renamed_genes/${i}_renamed_genes.fna"

done 


# ==========================================
# PART 2: Global Catalog Creation & Indexing
# ==========================================

# 9. Combining all renamed genes.fna files for indexing
cat \
    "${BASE_DIR}/renamed_genes/SRR37808534_renamed_genes.fna" \
    "${BASE_DIR}/renamed_genes/SRR37808535_renamed_genes.fna" \
    "${BASE_DIR}/renamed_genes/SRR37808538_renamed_genes.fna" \
    "${BASE_DIR}/renamed_genes/SRR37808539_renamed_genes.fna" \
    > "${BASE_DIR}/combined_genes.fna"

# 10. Salmon indexing of the combined gene catalog
/home/salmon-1.9.0_linux_x86_64/bin/salmon index \
    -t "${BASE_DIR}/combined_genes.fna" \
    -i "${BASE_DIR}/combined_genes_salmon_index" \
    > "${BASE_DIR}/combined_salmon_index.log" 2>&1


# ==========================================
# PART 3: Quantification Loop 
# ==========================================
for i in SRR37808534 SRR37808535 SRR37808538 SRR37808539
do
    echo "=== Quantifying sample: ${i} ==="
    
    /home/salmon/salmon-1.9.0_linux_x86_64/bin/salmon quant \
        -i "${BASE_DIR}/combined_genes_salmon_index" \
        -l A \
        -r "${BASE_DIR}/${i}_sortmerna/out/non_rRNA.fq.gz" \
        -o "${BASE_DIR}/${i}_salmon_quant" \
        > "${BASE_DIR}/${i}_salmon_quant.log" 2>&1
done
