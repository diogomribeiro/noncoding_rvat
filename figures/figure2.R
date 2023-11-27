
library(data.table)
source("~/git/noncoding_rvat/figures/manhattan.R")

################
# ABC
################
png("~/git/noncoding_rvat/figures/manhattanABC.png",1800,350)
mergedData = fread("~/git/noncoding_rvat/figures/data/fig2_ABC_skato_manhattan.tsv.gz")
manhattan(mergedData, chr="chro",bp="tss",snp="GENENAME",p="PVAL", col = c("#bae4b3", "#006d2c"), annotatePval = 1e-9, ylim = c(0,36), suggestiveline = 5.062) # ABC: green
dev.off()

################
# CRD
################
png("~/git/noncoding_rvat/figures/manhattanCRD.png",1800,350)
mergedData = fread("~/git/noncoding_rvat/figures/data/fig2_CRD_skato_manhattan.tsv.gz")
manhattan(mergedData, chr="chro",bp="tss",snp="GENENAME",p="PVAL", col = c("#bcbddc", "#756bb1"), annotatePval = 1e-9, ylim = c(0,51), suggestiveline = 4.434) # CRD: purple
dev.off()

################
# HIC
################
png("~/git/noncoding_rvat/figures/manhattanHIC.png",1800,350)
mergedData = fread("~/git/noncoding_rvat/figures/data/fig2_HIC_skato_manhattan.tsv.gz")
manhattan(mergedData, chr="chro",bp="tss",snp="GENENAME",p="PVAL", col = c("#fcae91", "#a50f15"), annotatePval = 1e-9, ylim = c(0,51), suggestiveline = 4.217) # HIC: red
dev.off()

################
# NC
################
png("~/git/noncoding_rvat/figures/manhattanNC.png",1800,350)
mergedData = fread("~/git/noncoding_rvat/figures/data/fig2_NC_skato_manhattan.tsv.gz")
manhattan(mergedData, chr="chro",bp="tss",snp="GENENAME",p="PVAL", col = c("#9ecae1", "#2171b5"), annotatePval = 1e-9, ylim = c(0,51), suggestiveline = 4.360) # NC: blue
dev.off()

################
# NC - conditional
################
png("~/git/noncoding_rvat/figures/manhattanNC_conditional.png",1800,350)
mergedData = fread("~/git/noncoding_rvat/figures/data/fig2_NC_conditional_skato_manhattan.tsv.gz")
manhattan(mergedData, chr="chro",bp="tss",snp="GENENAME",p="PVAL", col = c("#9ecae1", "#2171b5"), annotatePval = 1e-9, ylim = c(0,51), suggestiveline = 4.360) # NC: blue
dev.off()

################
# CDS
################
png("~/git/noncoding_rvat/figures/manhattanCDS.png",1800,350)
mergedData = fread("~/git/noncoding_rvat/figures/data/fig2_CDS_skato_manhattan.tsv.gz")
manhattan(mergedData, chr="chro",bp="tss",snp="GENENAME",p="PVAL", col = c("#d9d9d9", "#636363"), annotatePval = 1e-9, ylim = c(0,51), suggestiveline = 4.146) # CDS: grey
dev.off()

################
# CDS - conditional
################
png("~/git/noncoding_rvat/figures/manhattanCDS_conditional.png",1800,350)
mergedData = fread("~/git/noncoding_rvat/figures/data/fig2_CDS_conditional_skato_manhattan.tsv.gz")
manhattan(mergedData, chr="chro",bp="tss",snp="GENENAME",p="PVAL", col = c("#d9d9d9", "#636363"), annotatePval = 1e-9, ylim = c(0,51), suggestiveline = 4.146) # CDS: grey
dev.off()


