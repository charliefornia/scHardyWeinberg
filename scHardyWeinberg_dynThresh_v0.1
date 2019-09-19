library(Seurat)
library(matrixStats)

# raw expression data (UMIs) is binarized as: lower third => not expressed (ignoring top 1% of values)

# STEP1
# load seurat and choose replicate
load("~/Seurat.Robj")
# scan for names using: unique(seurat_obj@meta.data$orig.ident)

# STEP2
# subset Seurat object
seurat_obj_rep <- subset(INPUT_seurat_obj, subset = orig.ident == "rep1")

# combine _LTR with the main TE (note that additional LTRs like e.g. ...A_LTR will NOT be combined)
seurat_obj_rep@assays$RNA@counts[grep("TE-", rownames(seurat_obj_rep@assays$RNA@counts), value = T),] -> TE_mat
seurat_obj_rep@assays$RNA@counts[grep("TE-", rownames(seurat_obj_rep@assays$RNA@counts), value = T, invert = T),] -> GENE_mat
rownames(TE_mat) <- gsub("-LTR", "", rownames(TE_mat))
t(sapply(by(TE_mat,rownames(TE_mat),colSums),identity)) -> new_TE_mat
rbind(GENE_mat, new_TE_mat) -> new_mat
rm(TE_mat, GENE_mat, new_TE_mat)


# STEP3
# define functions

# this function gets the threshold that separates lower third of values from the rest.
# [inputvalue] is the number of top-hits that is ignored when calculating the lower third

get_nth <- function(inputvector, inputvalue) {
  return(max(round((tail(head(sort(inputvector, decreasing = T), inputvalue),1))/3,0),1))
}

run_scHW <- function(seurat_mat) {
  genes_all <- rownames(seurat_mat)
  # number of cells that are ignored when caculating the lower third
  num_cells <- round((ncol(seurat_mat)/100),digits = 0)
  # apply function to row-by-row, with number of cells to ignore as variable input
  t(apply(seurat_mat,1,get_nth, num_cells)) -> tmp.cutoff
  for(i in 1:length(tmp.cutoff)) {
    # set all values BELOW cutoff value to 0
    seurat_mat[i,seurat_mat[i,] <= tmp.cutoff[i]] <- 0
    # set all values ABOVE cutoff value to 1
    seurat_mat[i,seurat_mat[i,] > tmp.cutoff[i]] <- 1
  }
  # now, the rowMeans is the same as the proportion of cells expressing each gene at a level above the cutoff value
  tmp.P1.genes <- rowMeans(seurat_mat)
  names(tmp.P1.genes) <- genes_all
  # prepare object that will be proportion of cells expressing each pair of genes
  tmp.X11 <- matrix(NA, ncol = length(genes_all), nrow = length(genes_all))
  colnames(tmp.X11) <- genes_all
  rownames(tmp.X11) <- genes_all
  # loop through matric row-by-row
  for(k in 1:length(genes_all)) {
    # create matrix where the k-th line is added to every line of gene matrix
    tmp.mat01 <- sweep(seurat_mat, 2, seurat_mat[k,], "+")
    # now, everytime both genes are expressed, the value will be 2
    tmp.mat01[tmp.mat01 == 1] <- 0
    tmp.mat01[tmp.mat01 == 2] <- 1
    # again, the rowmeans will be the proportion of cells expressing both genes at a level above the cutoff value
    tmp.X11[genes_all[k],] <- rowMeans(tmp.mat01)
    rm(tmp.mat01)
  }
  # for easier calculations, create a n_gene x n_gene matrix - this provides the P of the first gene (P1), when replicated by row:
  tmp.P1.genes.matrix <- matrix(tmp.P1.genes, nrow = length(genes_all), ncol = length(genes_all), byrow = T)
  # and of the second gene (P2) when replicated by column:
  tmp.P2.genes.matrix <- matrix(tmp.P1.genes, nrow = length(genes_all), ncol = length(genes_all), byrow = F)
  # this line is the caculation of the CD coefficient
  tmp.HARDY_coef <- (tmp.X11-(tmp.P1.genes.matrix*tmp.P2.genes.matrix))^2/(tmp.P1.genes.matrix*(1 - tmp.P1.genes.matrix)*tmp.P2.genes.matrix*(1 - tmp.P2.genes.matrix))
  rownames(tmp.HARDY_coef) <- genes_all
  colnames(tmp.HARDY_coef) <- genes_all
  diag(tmp.HARDY_coef) <- NaN
  return(tmp.HARDY_coef)
}

# STEP4
# run function and save output to file

run_scHW(seurat_obj = new_mat) -> SAMPLENAME_HWCoef_rep1
save(SAMPLENAME_HWCoef_rep1, file = "~/SAMPLENAME_HWCoef_rep1.Robj")

# STEP5
# combine all HWCoef objects and run this
load("~/SAMPLENAME_HWCoef_rep1.Robj")
load("~/SAMPLENAME_HWCoef_rep2.Robj")
load("~/SAMPLENAME_HWCoef_rep3.Robj")
load("~/SAMPLENAME_HWCoef_rep4.Robj")
load("~/SAMPLENAME_HWCoef_rep5.Robj")
load("~/SAMPLENAME_HWCoef_rep6.Robj")
load("~/SAMPLENAME_HWCoef_rep7.Robj")
load("~/SAMPLENAME_HWCoef_rep8.Robj")
X <- list(rowRanks(SAMPLENAME_HWCoef_rep1), rowRanks(SAMPLENAME_HWCoef_rep2), rowRanks(SAMPLENAME_HWCoef_rep3), rowRanks(SAMPLENAME_HWCoef_rep4), rowRanks(SAMPLENAME_HWCoef_rep5), rowRanks(SAMPLENAME_HWCoef_rep6), rowRanks(SAMPLENAME_HWCoef_rep7), rowRanks(SAMPLENAME_HWCoef_rep8)) 
Y <- do.call(cbind, X)
SAMPLENAME_HWCoef_combined <- array(Y, dim=c(dim(X[[1]]), length(X)))
rownames(SAMPLENAME_HWCoef_combined) <- rownames(SAMPLENAME_HWCoef_rep1.Robj)
colnames(SAMPLENAME_HWCoef_combined) <- colnames(SAMPLENAME_HWCoef_rep1.Robj)
save(SAMPLENAME_HWCoef_combined, file = "~/SAMPLENAME_HWCoef_combined.Robj")

# STEP6
# test a gene pair

load("~/SAMPLENAME_HWCoef_combined.Robj")
mu_vec <- apply(SAMPLENAME_HWCoef_combined, c(1), mean, na.rm = TRUE)

read.csv("~/Dropbox/CloudDesktop/input.tsv", header = F, sep = "\t") -> input_list

for(i in 1:nrow(input_list)) {
  if(length(unique(NOV2017_DS_combined_HWCoef_combined[input_list[i,1], input_list[i,2],]))>2 & !is.na(as.numeric(mu_vec[as.character(input_list[i,1])]))){  
    t.test(NOV2017_DS_combined_HWCoef_combined[input_list[i,1], input_list[i,2],], mu = as.numeric(mu_vec[as.character(input_list[i,1])]), alternative = "greater")$p.value -> input_list[i,3]
  }
}

p.adjust(input_list[,3], method = "BH") -> input_list[,4]

# for random gene-te pairs
for(r in seq(5,23,2)){
  for(i in 1:nrow(input_list)) {
    sample(1:length(NOV2017_DS_combined_HWCoef_combined[,1,1]), size = 1) -> random_number
    rownames(NOV2017_DS_combined_HWCoef_combined)[random_number] -> input_list[i,r]
    if(length(unique(NOV2017_DS_combined_HWCoef_combined[input_list[i,1], rownames(NOV2017_DS_combined_HWCoef_combined)[random_number],]))>2 & !is.na(as.numeric(mu_vec[as.character(input_list[i,1])]))){  
      t.test(NOV2017_DS_combined_HWCoef_combined[input_list[i,1], rownames(NOV2017_DS_combined_HWCoef_combined)[random_number],], mu = as.numeric(mu_vec[as.character(input_list[i,1])]), alternative = "greater")$p.value -> input_list[i,r+1]
    }
  }
}

write.table(input_list, file = "~/Dropbox/CloudDesktop/output.tsv", quote = F, row.names = F, col.names = FALSE, sep = "\t")
