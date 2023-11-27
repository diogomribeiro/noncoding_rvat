
library(data.table)
library(ggplot2)

wantedTraits = fread("~/git/noncoding_rvat/supporting_data/42_traits.txt", header = F)

normalData = fread("~/burdenNC/data/results/200k/NC_SKATO_process_all.tsv.gz")
# normalData = fread("~/burdenNC/data/results/200k/CDS_SKATO_process_all.tsv.gz")

condData = fread("~/burdenNC/data/results/200k/conditional/NC_conditional_SKATO_process_all.tsv.gz")
# condData = fread("~/burdenNC/data/results/200k/conditional/CDS_conditional_SKATO_process_all.tsv.gz")

mergedData = merge(normalData[,.(CHROM,GENE,LOG10P,PVAL,FDR,TRAIT)],condData[,.(CHROM,GENE,LOG10P,PVAL,FDR,TRAIT)], by = c("CHROM","GENE","TRAIT"))

mergedData = mergedData[TRAIT %in% wantedTraits$V1]

# unique(mergedData[,.(GENE,TRAIT)])
# unique(mergedData[,.(TRAIT)])
# unique(mergedData[,.(GENE)])

# Gene-trait associations
# Before conditional
nrow(mergedData[FDR.x < 0.05])
nrow(mergedData[PVAL.x < 1e-9])
# After conditional
nrow(mergedData[FDR.y < 0.05])
nrow(mergedData[PVAL.y < 1e-9])

# % remaining
remainFDR5 = nrow(mergedData[FDR.x < 0.05][FDR.y < 0.05])
remainP1e9 = nrow(mergedData[PVAL.x < 1e-9][PVAL.y < 1e-9])
round(remainFDR5 * 100 / nrow(mergedData[FDR.x < 0.05]),2)
round(remainP1e9) * 100 / nrow(mergedData[PVAL.x < 1e-9])

png("fig4c_conditional_nc.png",750,700)
# png("fig4d_conditional_cds.png",750,700)

# Plot
test = cor.test(mergedData$LOG10P.x,mergedData$LOG10P.y, method = "spearman")
text = paste("Spearman R:",round(test$estimate,2), "P-value",format.pval(test$p.value,2))
ggplot(data = mergedData[LOG10P.y > 1 | LOG10P.x > 1], aes(LOG10P.x, LOG10P.y)) + 
  geom_point() +
  # geom_smooth(method = "lm") + 
  geom_abline(slope = 1, intercept = 0, color = "grey", linetype = "dashed", size = 2) +
  annotate("text", x = Inf, y = Inf, label = text, hjust = 1.2, vjust = 1.5, size = 8, fontface = "bold"  ) +
  ylab("-Conditional P (-log10)") +
  xlab("No conditional P (-log10)") +
  theme_minimal() + 
  theme(text = element_text(size=30), legend.position = c(0.86,0.85),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=3), aspect.ratio = 1)  

dev.off()

############
## WRITE CONDITIONAL SUMMARY
############
# # add gene names mapping and trait names
# traitMapping = fread("/home/dribeiro/git/burdenNC/src/regenie/analysis/wanted_traits.txt", header = F)
# traitMapping$TRAIT = data.table(unlist(lapply(traitMapping$V1, function(x) unlist(strsplit(x,"[_]"))[1])))$V1
# geneMapping = fread("/home/dribeiro/git/burdenNC/src/saige/annotations/gencode_id_mapping.tsv")
# mergedData = merge(mergedData, traitMapping[,.(TRAIT,V2)], by = "TRAIT")
# mergedData = merge(mergedData, geneMapping[,.(V4,V6)], by.x = "GENE", by.y = "V4")
# colnames(mergedData) = c("GENE","TRAIT","CHROM","NOCOND_LOG10P","NOCOND_PVAL","NOCOND_FDR","COND_LOG10P","COND_PVAL","COND_FDR","TRAIT_NAME","GENE_NAME")
#write.table(mergedData,"~/git/burdenNC/src/regenie/analysis/conditional/NC.conditional.summary.tsv",quote=F,sep="\t",row.names=F)
#write.table(mergedData[NOCOND_FDR < 0.05 | COND_FDR < 0.05],"~/git/burdenNC/src/regenie/analysis/conditional/NC.conditional.FDR5.tsv",quote=F,sep="\t",row.names=F)
# write.table(mergedData,"~/git/burdenNC/src/regenie/analysis/conditional/CDS.conditional.summary.tsv",quote=F,sep="\t",row.names=F)
# write.table(mergedData[NOCOND_FDR < 0.05 | COND_FDR < 0.05],"~/git/burdenNC/src/regenie/analysis/conditional/CDS.conditional.FDR5.tsv",quote=F,sep="\t",row.names=F)


