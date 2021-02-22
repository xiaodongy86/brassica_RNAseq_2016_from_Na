
.libPaths()
.libPaths(c("/usr/lib64/R/library","/usr/share/R/library","/home/xzy50/R/x86_64-redhat-linux-gnu-library/3.6"))
.libPaths()
library(MethylIT)
library(edgeR)
library(DESeq2)
library(genefilter)

setwd("/data5/xzy50/brassica/RNAseq_2016/QoRTs_out")
sampleFiles <- c("CMS_L_NP_1_QC.geneCounts.formatted.for.DESeq.txt.gz",
                 "CMS_L_NP_2_QC.geneCounts.formatted.for.DESeq.txt.gz",
                 "CMS_L_NP_3_QC.geneCounts.formatted.for.DESeq.txt.gz",
                 "CMS_L_P_1_QC.geneCounts.formatted.for.DESeq.txt.gz",
                 "CMS_L_P_2_QC.geneCounts.formatted.for.DESeq.txt.gz",   
                 "CMS_L_P_3_QC.geneCounts.formatted.for.DESeq.txt.gz")

sampleName <- c("CMS_L_NP_1","CMS_L_NP_2","CMS_L_NP_3","CMS_L_P_1","CMS_L_P_2","CMS_L_P_3")

sampleCondition <- factor( c( "CMS_L_NP", "CMS_L_NP", "CMS_L_NP",
                              "CMS_L_P", "CMS_L_P", "CMS_L_P"),
                           levels = c("CMS_L_NP", "CMS_L_P"))

sampleTable <- data.frame( sampleName = sampleName,
                           fileName = sampleFiles,
                           condition = sampleCondition )

deg.rnai <- DESeqDataSetFromHTSeqCount( sampleTable = sampleTable,
                                        directory = "./",
                                        design= ~ condition)

Bju_genome_gene_models <- import("/data5/xzy50/brassica/RNAseq_2016/Genome/Bju.genome.gtf")
head(Bju_genome_gene_models)
seqlevels(Bju_genome_gene_models)
gene = Bju_genome_gene_models[ Bju_genome_gene_models$type == "gene", c("ID","gene_id")]
#head(Tomato_ITAG3.0_gene_models)
#names(mcols(Tomato_ITAG3.0_gene_models))
#seqlevels( gene ) <- c( "0","1","2","3","4","5","6","7","8","9","10","11","12" )
#seqlevels(gene, pruning.mode = "coarse") <- c("0","1","2","3","4","5","6","7","8","9","10","11","12" )
#gene = sortBySeqnameAndStart(gene)
#gene$ID <- substr(gene$ID, start =6 ,stop =22)
#tail(rownames(deg.rnai))

idx <- match(rownames(deg.rnai), gene$gene_id)
sum(is.na(idx))
# 1900
GENES <- gene[na.omit(idx)]
#33925
idx <- match(gene$ID, rownames(deg.rnai))
sum(is.na(idx))
deg.rnai <- deg.rnai[na.omit(idx)]

genes <- data.frame(GENES)[, -4]

head(counts(deg.rnai))

# RRCtlP1 RRCtlP2 RRCtlP3 RDrP1 RDrP2 RDrP3
# Solyc01g005000.3    1651    1457    1057  3807  3771  4441
# Solyc01g005010.3     702     641     470   803   752   675
# Solyc01g005020.3    2604    1234    1224  2382  2403  3553
# Solyc01g005030.3     893     944     657  1108  1129  1509
# Solyc01g005060.3     115      86      28   173   131    68
# Solyc01g005070.3       0       0       0     0     0     0

exp_count <- DGEList(counts = counts(deg.rnai), group = sampleCondition,
                     genes = data.frame(genes = genes, length = width(GENES)))

# deg.rnai <- DESeqDataSetFromHTSeqCount( sampleTable = sampleTable,
#                                         directory = "./",
#                                         design= ~ condition)

exp_count$samples

keep <- filterByExpr(exp_count, min.count = 5)
sum(keep)
exp_count <- exp_count[keep, , keep.lib.sizes=FALSE]
dim(exp_count$counts)

#calculate normalization factor
exp_count <- calcNormFactors(exp_count)

#estimate dispersion
exp_count <- estimateDisp(exp_count)

#gene differential analysis
et_exp_count <- exactTest(exp_count, pair = c("RRCtl", "RDr"))


#  log2fc 0.5 cut off

idx <- which(abs(et_exp_count$table$logFC) >= 0.5 & abs(et_exp_count$table$PValue) <= 0.1)

length(idx) # 2483

et_exp_count$table <- et_exp_count$table[idx,]

dim(et_exp_count$table) # 2483

##  with p value adjustment ，BH

et_exp_count$table$adj.pval <- p.adjust(et_exp_count$table$PValue, method = "BH")
idx <- which(et_exp_count$table$adj.pval < 0.1)
DEG <- et_exp_count[idx,]
dim(DEG$table) # 2483

library(seqinr)
library(Biostrings)
library(BiocParallel)
db.fa <- "/Araport11_genes.201606.pep.fasta"
dir.fa <- "/data/users/xzy50/TomatoBlast"
dir.db <- "/data/users/xzy50/TomatoBlast"

### set up the query sequence 
file <- "/data/users/xzy50/TomatoBlast/ITAG3.0_proteins.fasta"
seqs <- readAAStringSet(filepath = file, format = "fasta")

DEGs_RR_RDR <- rownames(DEG$table)
#DEGs_RR_RDR<- rownames(DEGs_0.5_1043)
tmp <- dir.db
#grep(DMG_13_vs_12_CG_208, names(seqs))
#?AAStringSet
keepers<-c()  
i=1
for(item in DEGs_RR_RDR) {
  tempy <- grep(item,names(seqs))
  if(is.integer(tempy)){
    keepers[i]<-tempy
    i=i+1
  }else{
    print("problem, is.integer is not an integer")
  }
}

query.seq <- seqs[na.omit(keepers)]

# query.seq <- seqs[1:20]
seq.name  <- substr(x = names(query.seq), start = 1, stop = 16)


blastp <- function(query.seq, dtb = NULL, seq.name = NULL, db.fa = NULL,
                   
                   dir.fa = NULL, tmp = getwd(), maxTargetSeqs = 2, numcode = 1,
                   
                   num.cores = 1L, tasks = 0L) {
  
  if (missing(query.seq)) stop("Provide query sequence(s) in FASTa format")
  
  if (is.null(db.fa) && is.null(dtb))
    
    stop("Provide sequence database or a fasta file to create it")
  
  
  
  if (!is.null(dtb)) {
    
    db <- suppressWarnings(try(system(paste0("blastdbcmd -db ", dtb,
                                             
                                             " -info"),
                                      
                                      intern = TRUE, ignore.stderr = TRUE),
                               
                               silent = TRUE))
    
    if (inherits(db, "try-error") || length(db) == 0)
      
      stop("The database provided is not valid")
    
  }
  
  
  
  if (is.null(dir.fa) && is.null(dtb)) stop("Provide the fasta file directory")
  
  else {
    
    if (is.null(dtb) && !is.null(db.fa)) {
      
      db <- suppressWarnings(try(system(paste0("blastdbcmd -db ", dir.db,
                                               
                                               db.fa, ".prot -info"),
                                        
                                        intern = TRUE,
                                        
                                        ignore.stderr = TRUE),
                                 
                                 silent = TRUE))
      
      if (!inherits(db, "try-error") && length(db) > 0)
        
        dtb <-paste0(dir.db, db.fa, ".prot")
      
    }
    
  }
  
  
  
  if (is.null(dtb) && !is.null(db.fa) && !is.null(dir.fa)) {
    
    newdb <- paste0("makeblastdb -in ", dir.fa, db.fa, " -dbtype prot ",
                    
                    "-parse_seqids -out ", dir.db, db.fa, ".prot ",
                    
                    "-title ", db.fa)
    
    system(newdb)
    
    dtb <- try(system(paste0("blastdbcmd -db ", dir.db, db.fa, ".prot -info")
                      
                      , intern = TRUE))
    
    if (inherits(dtb, "try-error") || length(dtb) == 0)
      
      stop("Database creation failed")
    
    dtb <-paste0(dir.db, db.fa, ".prot")
    
  }
  
  
  
  # === Auxiliar function to perform blastp through OS command  ===
  
  blast <- function(k, query, dtb, tmp, seq.name, maxTargetSeqs) {
    
    sname <- seq.name[k]
    
    writeXStringSet(query[k], filepath = paste0(tmp, "tmp", sname, ".fasta"))
    
    str1 = paste0("blastp -db ", dtb, " -query ", tmp, "tmp", sname, ".fasta")
    
    str2 = paste0("-out ", tmp, "tmp", sname,
                  
                  ".txt -outfmt '6 qseqid sseqid pident ",
                  
                  "qcovs bitscore score evalue' -max_target_seqs ",
                  
                  maxTargetSeqs)
    
    system(paste( str1, str2, sep = " " ))
    
    tmp1 <- try( read.delim(paste0(tmp, "tmp", sname,".txt"),
                            
                            header = FALSE ), silent = TRUE)
    
    if (!inherits(tmp1, "try-error")) {
      
      res <- DataFrame(tmp1)
      
      colnames(res) <- c("qseqid", "sseqid", "pident", "qcovs",
                         
                         "bitscore", "score", "evalue")
      
      file.remove(c(paste0(tmp, "tmp", sname,".txt"),
                    
                    paste0(tmp, "tmp", sname, ".fasta" ) ) )
      
      return(res)
      
    } else {
      
      file.remove(paste0(tmp, "tmp", sname, ".fasta" ))
      
      res = DataFrame(qseqid = NA, sseqid = NA, pident = NA,
                      
                      qcovs = NA, bitscore = NA, score = NA, evalue = NA)
      
    }
    
    return(res)
    
  }
  
  # ------------------------------------------------------------------------- #
  
  
  
  if (length(query.seq) > 1) {
    
    # Set parallel computation
    
    if (Sys.info()['sysname'] == "Linux") {
      
      bpparam <- MulticoreParam(workers=num.cores, tasks = tasks)
      
    } else bpparam <- SnowParam(workers = num.cores, type = "SOCK")
    
    num.seq <- 1:length(query.seq)
    
    res <- bplapply(num.seq, blast, query=query.seq, dtb=dtb, tmp=tmp,
                    
                    seq.name=seq.name, maxTargetSeqs=maxTargetSeqs,
                    
                    BPPARAM = bpparam)
    
    res <- do.call(rbind, res)
    
  } else res <- blast(1, query=query.seq, dtb=dtb, tmp=tmp,
                      
                      seq.name=seq.name, maxTargetSeqs)
  
  
  
  return(res)
  
}


blastp_DEGs_RR_RDR <- blastp(query.seq=query.seq, dtb = NULL, seq.name = seq.name, db.fa = db.fa, maxTargetSeqs = 1,
                             dir.fa = dir.fa, tmp = tmp,  num.cores = 30L, tasks = 30L)


##################    log FC1  DEGs_logFC_1
sum(is.na(blastp_DEGs_RR_RDR))
# 210
blastp_DEGs_RR_RDR_na.omit <-na.omit(blastp_DEGs_RR_RDR)
# 2768
blastp_DEGs_RR_RDR_na.omit_score_100 = blastp_DEGs_RR_RDR_na.omit[blastp_DEGs_RR_RDR_na.omit$score >100, ] 
# 2423




blastpIDs <- unique(blastp_DEGs_RR_RDR_na.omit_score_100$qseqid)
blastp_DEGs_RR_RDR_na.omit_score_100_df <- as.data.frame(blastp_DEGs_RR_RDR_na.omit_score_100)

#  only take a match for each tomato gene 
unique_match_DEGs_RR_RDR <-as.data.frame(blastp_DEGs_RR_RDR_na.omit_score_100_df[1,])  
i=1
for(item in blastpIDs){
  unique_match_DEGs_RR_RDR[i,] = blastp_DEGs_RR_RDR_na.omit_score_100_df[blastp_DEGs_RR_RDR_na.omit_score_100_df$qseqid==item,][1,]; i=i+1}

dim(unique_match_DEGs_RR_RDR)
#2172 DEGs

unique_match_DEGs_RR_RDR$qseqid <- substr(x = unique_match_DEGs_RR_RDR$qseqid , start = 1, stop = 16)
unique_match_DEGs_RR_RDR$sseqid <- substr(x = unique_match_DEGs_RR_RDR$sseqid , start = 1, stop = 9)

unique_match_DEGs_RR_RDR_logFC <- data.frame(unique_match_DEGs_RR_RDR, DEG$table[match(unique_match_DEGs_RR_RDR$qseqid,rownames(DEG$table)),])

AG_gff3 = import("/data/TAIR10_gff3/Arabidopsis_thaliana.TAIR10.38.gff3.gz")
AG_anotation = AG_gff3[ AG_gff3$type == "gene",  c( "gene_id", "Name", "description" ) ]
seqlevels( AG_anotation ) <- c( "1","2","3","4","5","M","C" )
seqlevels(AG_anotation, pruning.mode = "coarse") <- c("1", "2", "3", "4", "5")
#AG_anotation_df <- as.data.frame(mcols(AG_anotation))
AG_anotation_df <- as.data.frame(AG_anotation)


unique_match_DEGs_RR_RDR_logFC_with_anotation <- data.frame(unique_match_DEGs_RR_RDR_logFC, AG_anotation_df[match(unique_match_DEGs_RR_RDR_logFC$sseqid,AG_anotation_df$gene_id),])
dim(unique_match_DEGs_RR_RDR_logFC_with_anotation)
###  2172
write.csv(unique_match_DEGs_RR_RDR_logFC_with_anotation,file = "/data5/tomato_RNAseq/unique_match_2172_DEGs_RR_RDR_logFC0.5_with_anotation_5_9_2020.csv")

library(EnrichmentBrowser)
library(GenomicRanges)
library(data.table )

# ============= KEGG datasets ============
kegg.gs <- getGenesets( "ath", db = "kegg" )
str(head(kegg.gs))

# GO terms of a selected ontology as gene sets
# ontology of interest (BP, MF or CC):
# BP: Biological Process
# MF: Molecular Function
# CC: Cellular Component
ath.bp.gs <- getGenesets(org = "ath", db = "go", go.onto = "BP", go.mode = "GO.db")
str(head(ath.bp.gs))

ath.mf.gs <- getGenesets(org = "ath", db = "go", go.onto = "MF", go.mode = "GO.db")
str(head(ath.mf.gs))

ath.cc.gs <- getGenesets(org = "ath", db = "go", go.onto = "CC", go.mode = "GO.db")
str(head(ath.cc.gs))

#### *** Compilation of a gene regulatory network from KEGG pathways (already
#### done) !!!
ath.grn <- compileGRN("ath")
head(ath.grn)
#      FROM        TO          TYPE
# [1,] "AT1G48970" "AT1G04170" "+" 
# [2,] "AT1G48970" "AT2G18720" "+" 
# [3,] "AT1G48970" "AT2G40290" "+" 
# [4,] "AT1G48970" "AT3G07920" "+" 
# [5,] "AT1G48970" "AT5G01940" "+" 
# [6,] "AT1G48970" "AT5G05470" "+"



#########################################################################################
#
# --------------- Network Enrichment Analysis Test . Based on GO -----------------------
#
#########################################################################################
# Compute NEAT (Signorelli et al., 2016) with small correction, a test for network enrichment
# analysis between/from a first list of sets ('A sets') and/to a second list of sets ('B sets').
# neat is the R package that implements NEAT, the Network Enrichment Analysis Test which is 
# presented in Signorelli, M., Vinciotti, V., Wit, E. C. (2016). NEAT: an efficient network 
# enrichment analysis test. BMC Bioinformatics, 17:352.

library( neat ) # https://cran.r-project.org/web/packages/neat/vignettes/neat.html
source("/data/R_functions/.neat.R" )
# ATH_GOs <- read.delim("/media/robersy/OS/Users/Robersy/Documents/Work/TAIR10_gff3/GOs/ATH_GO_GOSLIM.txt", header = FALSE)
# # ATH_GOs <- read.delim("C:/Users/Robersy/Documents/Work/TAIR10_gff3/GOs/ATH_GO_GOSLIM.txt", header = FALSE)
# ATH_GOs = data.table( ATH_GOs )
# colnames( ATH_GOs ) <- c( "locus", "TAIR.accession", "object.name", "relationship.type", "GO.term", "GO.ID",
#                           "TAIR.Keyword.ID", "Aspect", "GOslim term", "Evidence.code", "Evidence.description", 
#                           "Evidence with", "Reference", "Annotator", "Date.annotated" )
# 
# # *** To build custom annotations ****
# locus = ATH_GOs$locus
# GOs = paste0( ATH_GOs$GO.ID, "_", ATH_GOs$GO.term )
# GO2Genes = by( locus, GOs, function(s) unique( as.character(s) ), simplify = T )
# GO2Genes = GO2Genes[1:length(GO2Genes)]
#ATn = rownames( DMG )

unique_match_DEGs_RR_RDR$sseqid <- substr(x = unique_match_DEGs_RR_RDR$sseqid , start = 1, stop = 9)
ATn = unique_match_DEGs_RR_RDR$sseqid
# ath.grn <- compile.grn.from.kegg( "/JOB/Work/KEGG/ath.zip" )
# head( ath.grn )
genesNet = unique( as.vector( ath.grn[ , 1:2 ] ) ) 

test = .neat( alist = list('AT.DEG' = ATn ), blist = ath.bp.gs, network = ath.grn[ , 1:2 ],
              nettype = 'directed', nodes = genesNet, alpha = 0.05 )
sum( test$pvalue <= 0.05 )
# ntest = test[ , ]
ntest = test[( test$nab > 0 & test$pvalue <= 0.05 ), ]
ntest = ntest[ (ntest$conclusion == "Overenrichment" ), ]

write.csv(ntest, file = "/data5/tomato_RNAseq/NEAT_GO_blastp_DEGs_RR_RDR_na.omit_score_100_total_2172_5_9_2020.csv")

GO.NetAB = ath.bp.gs[ match( as.character( ntest$B ), names( ath.bp.gs ) ) ]
GeneGO.Net = sapply( as.list( GO.NetAB ), function(s) s[ na.omit( match( ATn, s ) ) ] ) 


Netw = names( GeneGO.Net )
GenesGO.Net = c()
for( k in 1:length( GeneGO.Net ) ) {
  if( length( GeneGO.Net[[ k ]]) > 0 ) 
    GenesGO.Net = rbind( GenesGO.Net, data.frame( GeneID = GeneGO.Net[[ k ]], Network = Netw[ k ] ) )
}

GenesGO.Net = data.table( GenesGO.Net )
GenesGO.Nets = GenesGO.Net[ , list( Network = list( unique( as.character( Network ) ) ) ), by = GeneID ]
#DEGlist = rownames( DMG )
DEGlist = unique_match_DEGs_RR_RDR$sseqid








library(topGO)
library(EnrichmentBrowser)  
library(dplyr)
# sig.Genes <- function( topGOdata, goID ) {
#    go.genes <- genesInTerm( topGOdata, goID)[[ 1 ]]
#    sig.genes <- sigGenes( topGOdata )
#    sig.genes[ na.omit( match( go.genes, sig.genes ) ) ]
#}

# preparing allGenes object 
AG = import("/data/TAIR10_gff3/Arabidopsis_thaliana.TAIR10.38.gtf.gz")
gene = AG[ AG$type == "gene", c( "gene_id", "gene_biotype" ) ]
gene = gene[ gene$gene_biotype == "protein_coding", "gene_id" ]
seqlevels(gene, pruning.mode = "coarse") <- c("1", "2", "3", "4", "5")
seqlevels(gene)<- paste0("Chr", 1:5)
gene = sortBySeqnameAndStart(gene)

allgenes <- mcols(gene)
allgenes$seletion <- rep(0, nrow(allgenes))


#DMGs <- read.csv( "/data/users/xzy50/msh1_transgen_10_10_2018_final_output/updatedGLM_11_6_2018/overlap/DMG_gen1_all_gen_both_gb_promoter_with_DIMP_info.csv", header = TRUE )
rownames(unique_match_DEGs_RR_RDR_logFC) <- unique_match_DEGs_RR_RDR_logFC$qseqid
DEG = unique_match_DEGs_RR_RDR_logFC[ , c( "logFC", "adj.pval","sseqid" ) ]
dim(DEG)
#colnames( DMG ) <- c( "FC", "ADJ.PVAL" )
#rownames( DMG ) <- DMGs$X 
#geneList = DMG$FC
#names( geneList ) <- rownames(DMG)

idx <- na.omit(match(DEG$sseqid, allgenes$gene_id))
sum(is.na(idx))


allgenes$seletion[idx] <- 1 

# check the number of row has 1 , make sure it fits DMG number
allgenes[allgenes$seletion == 1,]
#DataFrame with 97 rows and 2 columns

geneList <- allgenes$seletion
names( geneList ) <- allgenes$gene_id

#' gene Selection Fun: function to specify which genes are interesting based on
#' the gene scores/ The function must have one argument, named 'allScore' and
#' must not depend on any attributes of this object.



topDiffGenes <- function(allScore) {
  # allScore: log2FC or p-value
  return( allScore > 0 )
}


sum(topDiffGenes(geneList)) # Selected top genes

# To shortening ath.bp.gs GO-names
# GO terms of a selected ontology as gene sets
# ontology of interest (BP, MF or CC):
# BP: Biological Process
# MF: Molecular Function
# CC: Cellular Component
ath.bp.gs <- getGenesets(org = "ath", db = "go", go.onto = "BP", go.mode = "GO.db")
str(head(ath.bp.gs))

ath.mf.gs <- getGenesets(org = "ath", db = "go", go.onto = "MF", go.mode = "GO.db")
str(head(ath.mf.gs))

ath.cc.gs <- getGenesets(org = "ath", db = "go", go.onto = "CC", go.mode = "GO.db")
str(head(ath.cc.gs))

# List of 6
# $ GO:0000002_mitochondrial_genome_maintenance: chr [1:3] "AT1G47720" "AT3G10140" "AT3G24320"
# $ GO:0000003_reproduction                    : chr [1:18] "AT1G18450" "AT1G54490" "AT2G06210" "AT2G38440" ...
# $ GO:0000012_single_strand_break_repair      : chr "AT1G08130"
# $ GO:0000023_maltose_metabolic_process       : chr [1:150] "AT1G02560" "AT1G03310" "AT1G03630" "AT1G03680" ...
# $ GO:0000024_maltose_biosynthetic_process    : chr "AT4G17090"
# $ GO:0000025_maltose_catabolic_process       : chr [1:2] "AT2G40840" "AT5G64860"

#kegg.gs <- getGenesets( "ath", db = "kegg" )
#str(head(kegg.gs))

x = sapply(names(ath.bp.gs),
           function(s) strsplit(s, split = "_")[[1]][1])
names(ath.bp.gs) <- unname(x)

str(head(ath.bp.gs))

# List of 6
# $ GO:0000002: chr [1:3] "AT1G47720" "AT3G10140" "AT3G24320"
# $ GO:0000003: chr [1:18] "AT1G18450" "AT1G54490" "AT2G06210" "AT2G38440" ...
# $ GO:0000012: chr "AT1G08130"
# $ GO:0000023: chr [1:150] "AT1G02560" "AT1G03310" "AT1G03630" "AT1G03680" ...
# $ GO:0000024: chr "AT4G17090"
# $ GO:0000025: chr [1:2] "AT2G40840" "AT5G64860"

# genes to GOs object , use in gene2GO argument


geneID2GO = inverseList(ath.bp.gs)
str(head(geneID2GO))


dmgs <- new("topGOdata", description = "DMGs all gen relax cutoff",
            ontology = "BP", allGenes = geneList, geneSel = topDiffGenes,
            annot = annFUN.gene2GO, gene2GO = geneID2GO, nodeSize = 10 )


# Enrichment tests. whichAlgorithms()  whichTests()

dmg.Fisher <- runTest( dmgs, algorithm = "classic", statistic = "fisher" )

dmg.KS <- runTest( dmgs, algorithm = "classic", statistic = "ks" )

dmg.KS.w01 <- runTest( dmgs, algorithm = "weight01", statistic = "ks" )

dmg.F.elim <- runTest( dmgs, algorithm = "elim", statistic = "fisher" )



# ======================== Analysis of results ============================= #

dmg.allRes <- GenTable( dmgs, classicFisher = dmg.Fisher, classicKS = dmg.KS,
                        
                        w01KS = dmg.KS.w01, F.elim = dmg.F.elim,
                        
                        orderBy = "w01KS", ranksOf = "classicFisher",
                        
                        topNodes = 50 )

class(dmg.allRes)

dmg.allRes[dmg.allRes$Significant >0,]


#GO_0010050 <- ath.bp.gs[match("GO:0010050", names(ath.bp.gs))][[1]]

dmg.allRes <- GenTable( dmgs, classicFisher = dmg.Fisher, classicKS = dmg.KS,
                        
                        w01KS = dmg.KS.w01, F.elim = dmg.F.elim,
                        
                        orderBy = "classicKS", ranksOf = "classicKS",
                        
                        topNodes = 50 )



dmg.allRes <- GenTable( dmgs, classicFisher = dmg.Fisher, classicKS = dmg.KS,
                        
                        w01KS = dmg.KS.w01, F.elim = dmg.F.elim,
                        
                        orderBy = "classicFisher", ranksOf = "classicFisher",
                        
                        topNodes = 50 )


dmg.allRes <- GenTable( dmgs, classicFisher = dmg.Fisher, classicKS = dmg.KS,
                        
                        w01KS = dmg.KS.w01, F.elim = dmg.F.elim,
                        
                        orderBy = "w01KS", ranksOf = "w01KS",
                        
                        topNodes = 50 )

dmg.allRes[dmg.allRes$Significant >2,]


dmg.allRes <- GenTable( dmgs, classicFisher = dmg.Fisher, classicKS = dmg.KS,
                        
                        w01KS = dmg.KS.w01, F.elim = dmg.F.elim,
                        
                        orderBy = "F.elim", ranksOf = "F.elim",
                        
                        topNodes = 25 )

dmg.allRes[dmg.allRes$Significant >2,]

dmg.allRes[dmg.allRes$Significant >3,]



GO_0009725<- ath.bp.gs[match("GO:0009725", names(ath.bp.gs))][[1]]
DMG[na.omit(match(GO_0009785,rownames(DMG))),]












