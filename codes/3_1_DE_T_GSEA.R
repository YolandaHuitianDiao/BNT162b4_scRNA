####################   Reused codes     ####################
# Last update: Sept 13, 2020
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
setwd(sum.dir)

###----- Find all de files
de.out.dir <- file.path(base_dir, "1_scanpy_scirpy", key_name, "1_DE")
de.files <- list.files(de.out.dir, pattern="*differential_clean.csv", full.names=TRUE,recursive=TRUE)
de.files <- rev(de.files)

use.cols <- c("gene_names", "t-test_overestim_var_logfc", "t-test_overestim_var_padj")

###----- Select gene signature
gs.files <- c("https://raw.githubusercontent.com/TCellResearchTeam/T_cell_signature_Reference/master/X_GeneSignatures_mm/2020_mouse_CD8_mm_sigs.csv",
              "https://raw.githubusercontent.com/TCellResearchTeam/T_cell_signature_Reference/master/X_GeneSignatures_mm/CD4-ifn_mm_sigs.csv",
              "https://raw.githubusercontent.com/TCellResearchTeam/T_cell_signature_Reference/master/X_GeneSignatures_mm/CD4_mm_sigs.csv",
              "https://raw.githubusercontent.com/TCellResearchTeam/T_cell_signature_Reference/master/X_GeneSignatures_mm/CD8_mm_sigs.csv",
              "https://raw.githubusercontent.com/TCellResearchTeam/T_cell_signature_Reference/master/X_GeneSignatures_mm/Exp391_CD8_mm_sigs.csv",
              "https://raw.githubusercontent.com/TCellResearchTeam/T_cell_signature_Reference/master/X_GeneSignatures_mm/IL2_mm_sigs.csv",
              "https://raw.githubusercontent.com/TCellResearchTeam/T_cell_signature_Reference/master/X_GeneSignatures_mm/Wolski_curated_CD8_mm_sigs.csv")
gs.pattern <- "_mm_sigs.csv"
fc.cutoff <- 1


completed.name <- file.path(sum.dir, "completed.csv")
if (file.exists(completed.name)){
  completed.df <- read_csv(completed.name) %>% column_to_rownames('...1')
} else {
  completed.df <- data.frame()
  for (gs.file.i in c("file", gs.files)) {
    gs.file.base <- gsub(gs.pattern, "",basename(gs.file.i))
    completed.df[gs.file.base] <- character(0)
  }
}


for (file.i in de.files) {
  if (! file.i %in% completed.df$file) {
    file.i.qc.row <- c(file.i)
    
    ###----- Create output names & output directory
    file.i.outpath <- gsub(de.out.dir, "", file.i)
    file.i.outpath <- gsub(".csv", "", file.i.outpath)
    file.i.outpath <- gsub("^/", "",file.i.outpath)
    
    file.i.simp.name <- tail(unlist(strsplit( file.i.outpath, "/")), 1)
    file.i.simp.name <- gsub("_differential", "", file.i.simp.name)
    file.i.simp.vec <- unlist(strsplit(file.i.simp.name, "_vs_"))
    
    # Create output directory
    file.i.dir <- file.path(sum.dir, paste(file.i.outpath, "_sep", sep=""))
    dir.create(file.i.dir, recursive=TRUE, showWarnings=FALSE)
    
    for (gs.file.i in gs.files) {
      gs.file.base <- gsub(gs.pattern, "",basename(gs.file.i))
      file.i.gs.dir <- file.path(file.i.dir, gs.file.base)
      dir.create(file.i.gs.dir, recursive=TRUE, showWarnings=FALSE)
      setwd(file.i.gs.dir)
      
      # Run GSEA analysis
      useGroup <- "t-test_overestim_var_logfc"
      outName <- file.i.simp.name
      
      success <- "yes"
      gsea <- tryCatch({
        GSEA_analysis(file.i, useGroup, outName, gs.file.i, fc.cutoff)
      }, error=function(cond) {
        message(paste(file.i, gs.file.i,sep='; '))
        return("no")
      })
      
      if (!is.null(gsea)){
        success <- ""
      }
      
      file.i.qc.row <- c(file.i.qc.row, success)
    }
    
    completed.df[nrow(completed.df) +1, ] <- file.i.qc.row
    write.csv(completed.df, completed.name, row.names=TRUE)
  }
}


