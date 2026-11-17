for i in SRR28387650 SRR28387651 SRR28387652 SRR28387657 SRR28387658 SRR28387660
do
    echo ${i}
    cd ${i}
    
    # 1. Quality Control
    /home/FastQC/fastqc ${i}_1p.fastq.gz ${i}_2p.fastq.gz

    # 2. HISAT2 Alignment
    /home/hisat2/hisat2 \
        -x /home/BIB/Databases/HumanGenome_HG37/hg19_NoChr_hisat2_index \
        -1 ${i}_1p.fastq.gz \
        -2 ${i}_2p.fastq.gz \
        -S ${i}.sam 

    # 3. SAM to BAM Conversion
    samtools view -b ${i}.sam -o ${i}.bam

    # 4. Sorting BAM
    samtools sort ${i}.bam -o ${i}_sorted.bam

    # 5. Indexing BAM
    samtools index ${i}_sorted.bam

    # 6. StringTie Quantification
    stringtie ${i}_sorted.bam \
        -G /home/Rnaseq/hg19.ncbiRefSeq.nochr.gtf \
        -o ${i}_transcripts.gtf \
        -A ${i}_gene_abund.tab \
        -e \
        -p 4

    cd ..
done


# 7. Create GTF location file for prepDE.py
echo "Creating GTF_location.txt..."
> GTF_location.txt
for i in SRR28387650 SRR28387651 SRR28387652 SRR28387657 SRR28387658 SRR28387660
do
    echo -e "${i}\t/home/Vaishnavi/Rnaseq/${i}/${i}_transcripts.gtf" >> GTF_location.txt
done


# 8. Generate Expression Count Matrices using prepDE.py
python3 /home/stringtie-2.2.1.Linux_x86_6/prepDE.py \
    -i GTF_location.txt \
    -g gene_count_matrix.csv \
    -t transcript_count_matrix.csv \
    -v
