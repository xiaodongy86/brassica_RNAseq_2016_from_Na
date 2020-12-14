install.packages("seqinr")
install.packages("Biostrings")

library(seqinr)
library(Biostrings)
library(BiocParallel)

dir.create("/Users/mackenzie/RProject/Brassica_juncea/")
setwd("/Users/mackenzie/RProject/Brassica_juncea")
getwd()
# download arabidopsis pep file 
download.file(url ="ftp://ftp.ensemblgenomes.org/pub/plants/release-44/fasta/arabidopsis_thaliana/pep/Arabidopsis_thaliana.TAIR10.cds.all.fa.gz",
              destfile = "/Users/mackenzie/RProject/Brassica_juncea/Arabidopsis_thaliana.TAIR10.pep.fasta", quiet=FALSE)

db.fa <- "/Araport11_genes.201606.pep.fasta"
dir.fa <- "/Users/mackenzie/RProject/Brassica_juncea"
dir.db <- "/Users/mackenzie/RProject/Brassica_juncea"

tmp <- dir.db


#download Brassica_juncea pep file 
download.file(url ="http://brassicadb.org/brad/datasets/pub/Genomes/Brassica_juncea/V1.5/Bju.chr.modified_id.pep.fa.gz",
              destfile = "/Users/mackenzie/RProject/Brassica_juncea/Bju.chr.modified_id.pep.fa", quiet=FALSE)
file <- "/Users/mackenzie/RProject/Brassica_juncea/Bju.chr.modified_id.pep.fa"
seqs <- readAAStringSet(filepath = file, format = "fasta")
names(seqs)

# [1] "BjuA000007" "BjuA000013" "BjuA000130" "BjuA000221" "BjuA000307" "BjuA000308" "BjuA000309" "BjuA000310" "BjuA000311"
# [10] "BjuA000312" "BjuA000313" "BjuA000328" "BjuA000337" "BjuA000338" "BjuA000339" "BjuA000340" "BjuA000341" "BjuA000342"
# [19] "BjuA000343" "BjuA000344" "BjuA000345" "BjuA000350" "BjuA000351" "BjuA000352" "BjuA000353" "BjuA000354" "BjuA000355"


seqAt <- readAAStringSet(filepath = "/Users/mackenzie/RProject/Brassica_juncea/Araport11_genes.201606.pep.fasta", format = "fasta")
names(seqAt)

#Bj_DEGs_np_vs_p <- read.csv("/Users/mackenzie/RProject/Brassica_juncea/Bj_DEGs_np_vs_p.csv")

head(DEGs_FC1_Pvalue0.05_CMS_NP_vs_P)

# flowering_time_genes <- read.csv("/Users/mackenzie/RProject/Brassica_juncea/flowering_time_genes.csv")
# head(flowering_time_genes)
# 
# flowering_time_genes$Gene.number <- substr(x = flowering_time_genes$Gene.number, start = 1, stop = 9)
# CR_1_vs_WT.DEG_with_AT_ID$sseqid <- substr(x = CR_1_vs_WT.DEG_with_AT_ID$sseqid, start = 1, stop = 9)
# 
# 
# IDX = na.omit(match(flowering_time_genes$Gene.number,CR_1_vs_WT.DEG_with_AT_ID$sseqid))

# CR_1_vs_WT.DEG_with_AT_ID_Flower_Genes <- CR_1_vs_WT.DEG_with_AT_ID[IDX,]

# CR_1_vs_WT.DEG_with_AT_ID_Flower_72_Genes <- data.frame(CR_1_vs_WT.DEG_with_AT_ID_Flower_Genes,flowering_time_genes[match(CR_1_vs_WT.DEG_with_AT_ID_Flower_Genes$sseqid,flowering_time_genes$Gene.number),])
# write.csv(CR_1_vs_WT.DEG_with_AT_ID_Flower_72_Genes, file = "/Users/mackenzie/RProject/Brassica_juncea/CR_1_vs_WT.DEG_with_AT_ID_Flower_72_Genes.csv")


dim(DEGs_FC1_Pvalue0.05_CMS_NP_vs_P)

DEGs_FC1_Pvalue0.05_CMS_NP_vs_P_name  <-rownames(DEGs_FC1_Pvalue0.05_CMS_NP_vs_P)

keepers<-c()  
i=1
for(item in DEGs_FC1_Pvalue0.05_CMS_NP_vs_P_name) {
  tempy <- grep(item,names(seqs))
  if(is.integer(tempy)){
    keepers[i]<-tempy
    i=i+1
  }else{
    print("problem, is.integer is not an integer")
  }
}


query.seq <- seqs[keepers]
seq.name <- names(seqs)


#query.seq <- seqs[1:2]
# seq.name  <- substr(x = names(query.seq), start = 27, stop = 49)

# 
# seq.name  <- substr(x = names(query.seq), start = 1, stop = 8)


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


blastp(query.seq=query.seq, dtb = NULL, seq.name = seq.name, db.fa = db.fa, dir.fa = dir.fa, tmp = tmp,
       maxTargetSeqs = 2)




