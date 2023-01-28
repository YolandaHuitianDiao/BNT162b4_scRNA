####################   Reused codes     ####################
# Last update: Oct24, 2022
# By: Huitian Diao
############################################################

library(stringr)

source("/mnt/efs/hdiao/Projects/Tools/function_GSEA.R")
base_dir <- "/mnt/efs/hdiao/Projects/BNT162b4_scRNA"

####################     Parameters     ####################
genome <- "mm"
key_name <- "1_T"
############################################################
sum.dir <- file.path(base_dir, "1_scanpy_scirpy", key_name, "1_DE_GSEA")
dir.create(sum.dir, showWarnings = FALSE, recursive=TRUE)

####################      Summarize     ####################
use.signatures.df <- read_csv(file.path(sum.dir, "use_signatures_and_simplified_names.csv"))
use.signatures <- use.signatures.df$use.signatures
use.signatures.simp <- use.signatures.df$use.signatures.simp

##########---------- Group v.s. Group
key <- 'group_vs_group'
gsea.file.base.dir <- file.path(sum.dir, key)
setwd(gsea.file.base.dir)

###----- Summarize signatures
gsea.cats <- c("Exp391_CD8", "CD8", "IL2")

cp.label.list <- list(c('T_multi_G9', 'T_multi_G7'), c('T_multi_G12', 'T_multi_G7'), 
                   c('T_multi_G12', 'T_multi_G9'))

# List of summary files - for plotting
cp.gsea.sum.out.files <- c()
cp.gsea.sum.out.files.simp.names <- c()
for (cp.label.vec in cp.label.list) {
  cp.file.base <- paste(cp.label.vec[1], "_vs_", cp.label.vec[2], 
                        "_differential_clean_sep", sep="")
  
  cp.file.df <- data.frame()
  
  for (gsea.cat in gsea.cats){
    cp.file.cat <- paste(cp.label.vec[1], "/", cp.file.base, "/", gsea.cat, "/", 
                         gsub("_differential_clean_sep", "_clean---", cp.file.base), 
                         gsea.cat, "_", genome,"_sigs.csv", sep="")
    cp.file.cat <- file.path(gsea.file.base.dir, cp.file.cat)
    cp.file.cat.df <- read_csv(cp.file.cat)
    cp.file.df <- rbind(cp.file.df, cp.file.cat.df)
  }
  
  cp.file.df <- cp.file.df %>% distinct(ID, .keep_all=TRUE) %>% arrange(p.adjust)
  
  cp.csv.name <- paste(gsub("T_multi_", "", cp.label.vec[1]), "_vs_", 
                       gsub("T_multi_", "", cp.label.vec[2]), "_Slt-sigs.csv", 
                       sep="")
  cp.gsea.sum.out.files <- c(cp.gsea.sum.out.files, cp.csv.name)
  cp.gsea.sum.out.files.simp.names <- c(cp.gsea.sum.out.files.simp.names, 
                                        gsub("_Slt-sigs.csv","",cp.csv.name))
  
  write.csv(cp.file.df, cp.csv.name, row.names=FALSE)
}


###----- Summary & Plots
GSEA_sum(cp.gsea.sum.out.files, key, 22, 12, 
         TRUE, cp.gsea.sum.out.files.simp.names,
         TRUE, cp.gsea.sum.out.files.simp.names,
         TRUE, use.signatures, 
         TRUE, use.signatures.simp,
         FALSE)


##########---------- Cluster v.s. Cluster
key <- 'cluster_vs_rest'
gsea.file.base.dir <- file.path(sum.dir, key)
setwd(gsea.file.base.dir)

###----- Summarize signatures
gsea.cats <- c("Exp391_CD8", "CD8", "IL2")

cp.label.list <- list('0','1','2','3','4','5','6')

# List of summary files - for plotting
cp.gsea.sum.out.files <- c()
cp.gsea.sum.out.files.simp.names <- c()
for (cp.label.vec in cp.label.list) {
  cp.file.base <- paste(cp.label.vec,"_differential_clean_sep", sep="")
  
  cp.file.df <- data.frame()
  
  for (gsea.cat in gsea.cats){
    cp.file.cat <- paste(cp.file.base, "/", gsea.cat, "/", 
                         gsub("_differential_clean_sep", "_clean---", cp.file.base), 
                         gsea.cat, "_", genome,"_sigs.csv", sep="")
    cp.file.cat <- file.path(gsea.file.base.dir, cp.file.cat)
    cp.file.cat.df <- read_csv(cp.file.cat)
    cp.file.df <- rbind(cp.file.df, cp.file.cat.df)
  }
  
  cp.file.df <- cp.file.df %>% distinct(ID, .keep_all=TRUE) %>% arrange(p.adjust)
  
  cp.csv.name <- paste(gsub("T_multi_", "", cp.label.vec[1]), "_vs_", 
                       gsub("T_multi_", "", cp.label.vec[2]), "_Slt-sigs.csv", 
                       sep="")
  cp.gsea.sum.out.files <- c(cp.gsea.sum.out.files, cp.csv.name)
  cp.gsea.sum.out.files.simp.names <- c(cp.gsea.sum.out.files.simp.names, 
                                        gsub("_Slt-sigs.csv","",cp.csv.name))
  
  write.csv(cp.file.df, cp.csv.name, row.names=FALSE)
}


###----- Summary & Plots
GSEA_sum(cp.gsea.sum.out.files, key, 22, 15, 
         TRUE, cp.gsea.sum.out.files.simp.names,
         TRUE, cp.gsea.sum.out.files.simp.names,
         TRUE, use.signatures, 
         TRUE, use.signatures.simp,
         FALSE)













