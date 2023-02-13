#!/usr/bin/env Rscript
# Theo Tricou

require("data.table")
require("phangorn")
require("parallel")


# read tree speceies tree with ale ID
tree <- read.tree('aletree')

# create list of node and leaves nem
nodes_list <- c(tree$tip.label, tree$node.label)

# list reconciliaiton files
rec_list <- Sys.glob('*uml_rec')

ancestral_duplication <- function(leaf, rec, fam){
  dup_score = sum(rec[c(grep(leaf,rec[,1]),
  as.numeric(nodes_list[Ancestors(tree, which(nodes_list == leaf))])+1), 2])
  return(c(ID = leaf, score_dup = dup_score))
}



info_by_gene <- function(x){
  rec <- as.data.frame(fread(x, skip= "S_terminal_branch"))[,-1]
  # sum duplication frequency across all ancestor nodes of a leaf
  res <- as.data.frame(do.call('rbind',lapply(tree$tip.label, function(leaf) ancestral_duplication(leaf, rec, fam))))
  # retrive gene ID from rec file name
  res$FAM <- paste(strsplit(strsplit(x, "[.]")[[1]][1], "_")[[1]][c(3,4)], collapse="_")
  return(res)
}

RES <- do.call('rbind', mclapply(rec_list, info_by_gene, mc.cores = 4))

write.table(RES_anc, file = "Outlier_Duplication.txt", quote = F)
