THREADS=6
GENOME_INDEX_DIR=/data5/xzy50/brassica/RNAseq_2016/Genome
FASTA_FILES=/data5/xzy50/brassica/RNAseq_2016/Genome/Bju.genome.fasta
GTF_FILE=/data5/xzy50/brassica/RNAseq_2016/Genome/Bju.genome.gff
ulimit -n 2048
#mkdir -p ${GENOME_INDEX_DIR}
STAR --runMode genomeGenerate --runThreadN ${THREADS} --genomeDir ${GENOME_INDEX_DIR} --genomeFastaFiles ${FASTA_FILES} --sjdbGTFfile ${GTF_FILE} 
