THREADS=6
GENOME_INDEX_DIR=/data5/xzy50/brassica/RNAseq_2016/Genome
FASTA_FILES=/data5/xzy50/brassica/RNAseq_2016/Genome/Bju.genome.fasta
GTF_FILE=/data5/xzy50/brassica/RNAseq_2016/Genome/Bju.genome.gff
ulimit -n 2048
#mkdir -p ${GENOME_INDEX_DIR}
#STAR --runMode genomeGenerate --runThreadN ${THREADS} --genomeDir ${GENOME_INDEX_DIR} --genomeFastaFiles ${FASTA_FILES} --sjdbGTFfile ${GTF_FILE} 
#Fatal INPUT FILE error, no exon lines in the GTF file: /data5/xzy50/brassica/RNAseq_2016/Genome/Bju.genome.gff
#Solution: check the formatting of the GTF file, it must contain some lines with exon in the 3rd column.
#          Make sure the GTF file is unzipped.
#          If exons are marked with a different word, use --sjdbGTFfeatureExon .

# the name of exon is mRNA
STAR --runMode genomeGenerate --runThreadN ${THREADS} --genomeDir ${GENOME_INDEX_DIR} --genomeFastaFiles ${FASTA_FILES} --sjdbGTFfile ${GTF_FILE} --sjdbGTFfeatureExon mRNA
