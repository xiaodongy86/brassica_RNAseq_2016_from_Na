THREADS=10
GENOME_INDEX_DIR=/data5/xzy50/brassica/RNAseq_2016/STAR_INDEX
FASTA_FILES=/data5/xzy50/brassica/RNAseq_2016/Genome/Bju.genome.fasta
GTF_FILE=/data5/xzy50/brassica/RNAseq_2016/Genome/Bju.genome.gtf

BASE_DIR=/data5/xzy50/brassica/RNAseq_2016/trim_galore
OUT_DIR=//data5/xzy50/brassica/RNAseq_2016/STAR_align_out
ulimit -n 2048
ts STAR --genomeDir ${GENOME_INDEX_DIR} --readFilesCommand zcat --readFilesIn  ${BASE_DIR}/CMS_L_NP_1/0CL300001525_L01_12_1_trimmed.fq.gz --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 150000000000 --twopassMode Basic --outFilterMultimapNmax 1 --quantMode TranscriptomeSAM --runThreadN ${THREADS} --outFileNamePrefix ${OUT_DIR}/CMS_L_NP_1/
ts STAR --genomeDir ${GENOME_INDEX_DIR} --readFilesCommand zcat --readFilesIn  ${BASE_DIR}/CMS_L_NP_2/0CL300001525_L01_13_1_trimmed.fq.gz --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 150000000000 --twopassMode Basic --outFilterMultimapNmax 1 --quantMode TranscriptomeSAM --runThreadN ${THREADS} --outFileNamePrefix ${OUT_DIR}/CMS_L_NP_2/
ts STAR --genomeDir ${GENOME_INDEX_DIR} --readFilesCommand zcat --readFilesIn  ${BASE_DIR}/CMS_L_NP_3/0CL300001525_L01_14_1_trimmed.fq.gz --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 150000000000 --twopassMode Basic --outFilterMultimapNmax 1 --quantMode TranscriptomeSAM --runThreadN ${THREADS} --outFileNamePrefix ${OUT_DIR}/CMS_L_NP_3/
ts STAR --genomeDir ${GENOME_INDEX_DIR} --readFilesCommand zcat --readFilesIn  ${BASE_DIR}/CMS_L_P_1/0CL300001525_L01_15_1_trimmed.fq.gz --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 150000000000 --twopassMode Basic --outFilterMultimapNmax 1 --quantMode TranscriptomeSAM --runThreadN ${THREADS} --outFileNamePrefix ${OUT_DIR}/CMS_L_P_1/
ts STAR --genomeDir ${GENOME_INDEX_DIR} --readFilesCommand zcat --readFilesIn  ${BASE_DIR}/CMS_L_P_2/0CL300001525_L01_16_1_trimmed.fq.gz --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 150000000000 --twopassMode Basic --outFilterMultimapNmax 1 --quantMode TranscriptomeSAM --runThreadN ${THREADS} --outFileNamePrefix ${OUT_DIR}/CMS_L_P_2/
ts STAR --genomeDir ${GENOME_INDEX_DIR} --readFilesCommand zcat --readFilesIn  ${BASE_DIR}/CMS_L_P_3/0CL300001525_L01_17_1_trimmed.fq.gz --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 150000000000 --twopassMode Basic --outFilterMultimapNmax 1 --quantMode TranscriptomeSAM --runThreadN ${THREADS} --outFileNamePrefix ${OUT_DIR}/CMS_L_P_3/
