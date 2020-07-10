GTF_FILE=/data5/xzy50/brassica/RNAseq_2016/Genome/Bju.genome.gtf
OUT_DIR=/data5/xzy50/brassica/RNAseq_2016/QoRTs_out
BASE_DIR=/data5/xzy50/brassica/RNAseq_2016/STAR_align_out

ts QoRTs QC --stranded --minMAPQ 25 --SingleEnd --maxReadLength 50 --runFunctions GeneCalcs,JunctionCalcs,writeGeneCounts,writeDEXSeq,writeGeneBody,writeDESeq ${BASE_DIR}/CMS_L_NP_1/Aligned.sortedByCoord.out.bam ${GTF_FILE} ${OUT_DIR}/CMS_L_NP_1/

ts QoRTs QC --stranded --minMAPQ 25 --SingleEnd --maxReadLength 50 --runFunctions GeneCalcs,JunctionCalcs,writeGeneCounts,writeDEXSeq,writeGeneBody,writeDESeq ${BASE_DIR}/CMS_L_NP_2/Aligned.sortedByCoord.out.bam ${GTF_FILE} ${OUT_DIR}/CMS_L_NP_2/

ts QoRTs QC --stranded --minMAPQ 25 --SingleEnd --maxReadLength 50 --runFunctions GeneCalcs,JunctionCalcs,writeGeneCounts,writeDEXSeq,writeGeneBody,writeDESeq ${BASE_DIR}/CMS_L_NP_3/Aligned.sortedByCoord.out.bam ${GTF_FILE} ${OUT_DIR}/CMS_L_NP_3/

ts QoRTs QC --stranded --minMAPQ 25 --SingleEnd --maxReadLength 50 --runFunctions GeneCalcs,JunctionCalcs,writeGeneCounts,writeDEXSeq,writeGeneBody,writeDESeq ${BASE_DIR}/CMS_L_P_1/Aligned.sortedByCoord.out.bam ${GTF_FILE} ${OUT_DIR}/CMS_L_P_1/

ts QoRTs QC --stranded --minMAPQ 25 --SingleEnd --maxReadLength 50 --runFunctions GeneCalcs,JunctionCalcs,writeGeneCounts,writeDEXSeq,writeGeneBody,writeDESeq ${BASE_DIR}/CMS_L_P_2/Aligned.sortedByCoord.out.bam ${GTF_FILE} ${OUT_DIR}/CMS_L_P_2/

ts QoRTs QC --stranded --minMAPQ 25 --SingleEnd --maxReadLength 50 --runFunctions GeneCalcs,JunctionCalcs,writeGeneCounts,writeDEXSeq,writeGeneBody,writeDESeq ${BASE_DIR}/CMS_L_P_3/Aligned.sortedByCoord.out.bam ${GTF_FILE} ${OUT_DIR}/CMS_L_P_3/
