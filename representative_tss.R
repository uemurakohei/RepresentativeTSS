# this is the representative tss selection using FANTOM CAT robust transcript data (see materials methods)
# the coordinates of TSS were downloaded from the FANTOM website (https://fantom.gsc.riken.jp/5/suppl/Hon_et_al_2016/data/assembly/lv3_robust/)
# the information of protein codoig genes was obtainded from https://fantom.gsc.riken.jp/cat/v1/#/genes

# loading source data -------------------

# loading FANTOM CAT robust transcript data
FANTOM_robust <- read.table("source_data/FANTOM_CAT.lv3_robust.bed")

# loading FANTOM CAT robust transcript data
coding_genes_fantom_robust <- read.csv("source_data/coding_genes_fantom_robust.csv")

# obtain gene ids
genes_list <- sapply(FANTOM_robust[,4], function(x) unlist(strsplit(x, split="\\|"))[1])

# choosing representative tss data using TIE score cutoff -------------------------------

all_ids_representative_tss <- c()
for(o_gene in unique(genes_list)){
  o_FANTOM_robust <- FANTOM_robust[grep(o_gene, genes_list), ]
  o_maxTIE <- o_FANTOM_robust[which.max(o_FANTOM_robust[, "V5"]), ]
  all_ids_representative_tss <- c(all_ids_representative_tss, which(o_maxTIE[,"V4"] == FANTOM_robust[,"V4"]))
}
representative_tss_from_TIE_score <- FANTOM_robust[all_ids_representative_tss, ]

# obtain protein coding genes data from the representative TSSs data
gene_id <- sapply(representative_tss_from_TIE_score[,4], function(x){unlist(strsplit(x, "\\|"))[1]})
pos_mRNA <- which(gene_id %in% coding_genes_fantom_robust[,"GeneID"])

# save the data
write.table(representative_tss_from_TIE_score[pos_mRNA, ],
  file="representative_tss_only_coding_genes.txt",
  col.names=F,
  row.names=F)
