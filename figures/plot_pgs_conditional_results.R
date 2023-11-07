
library(data.table)
library(ggplot2)

noCond = fread("zcat ~/burdenNC/data/results/200k/NC_SKATO_process_all.tsv.gz")
# noCond = fread("zcat ~/burdenNC/data/results/200k/CDS_SKATO_process_all.tsv.gz")

cond = fread("zcat ~/burdenNC/data/results/200k/conditional/PGS/NC.pgs_3traits.regenie.gz")
# cond = fread("zcat ~/burdenNC/data/results/200k/conditional/PGS/CDS.pgs_3traits.regenie.gz")

genoCond = fread("zcat ~/burdenNC/data/results/200k/conditional/NC_conditional_SKATO_process_all.tsv.gz")

cond$GENE = data.table(unlist(lapply(cond$ID, function(x) unlist(strsplit(x,"[.]"))[1])))$V1
cond$TRAIT = data.table(unlist(lapply(cond$TRAIT, function(x) unlist(strsplit(x,"[.]"))[2])))$V1
cond$TRAIT = data.table(unlist(lapply(cond$TRAIT, function(x) unlist(strsplit(x,"[_]"))[2])))$V1
cond = cond[TEST == "ADD-SKATO"][,.(CHROM,GENPOS,GENE,LOG10P,TRAIT)]

mergedData = merge(noCond[,.(CHROM,GENE,LOG10P,TRAIT)],cond[,.(CHROM,GENE,LOG10P,TRAIT)], by = c("CHROM","GENE","TRAIT"))

# Gene-trait associations
# Before conditional
nrow(mergedData[LOG10P.x > 5])
nrow(mergedData[LOG10P.x > 9])
# After conditional
nrow(mergedData[LOG10P.y > 5])
nrow(mergedData[LOG10P.y > 9])

100 * nrow(mergedData[LOG10P.y > LOG10P.x]) / nrow(mergedData)
100 * nrow(mergedData[LOG10P.x > 9][LOG10P.y > LOG10P.x]) / nrow(mergedData[LOG10P.x > 9])

##################
## PGS vs normal
##################
png("sup10a_pgs_conditional.png",750,700)
test = cor.test(mergedData$LOG10P.x,mergedData$LOG10P.y, method = "spearman")
text = paste("Spearman R:",round(test$estimate,2), "P-value",format.pval(test$p.value,2))
ggplot(data = mergedData[LOG10P.x > 1], aes(LOG10P.x, LOG10P.y)) + 
  geom_point() +
  # geom_smooth(method = "lm") + 
  geom_abline(slope = 1, intercept = 0, color = "grey", linetype = "dashed", size = 2) +
  annotate("text", x = Inf, y = Inf, label = text, hjust = 1.5, vjust = 1.5, size = 7, fontface = "bold"  ) +
  ylab("PGS conditional (-log10)") +
  xlab("No conditional (-log10)") +
  theme_minimal() + 
  theme(text = element_text(size=32), legend.position = c(0.86,0.85),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=3), aspect.ratio = 1)  

dev.off()

##################
## GENO vs normal
##################
png("sup10b_geno_conditional.png",750,700)

mergedData2 = merge(noCond[,.(CHROM,GENE,LOG10P,TRAIT)],genoCond[,.(CHROM,GENE,LOG10P,TRAIT)], by = c("CHROM","GENE","TRAIT"))
mergedData2 = mergedData2[TRAIT %in% unique(mergedData$TRAIT)]
100 * nrow(mergedData2[LOG10P.y > LOG10P.x]) / nrow(mergedData2)
100 * nrow(mergedData2[LOG10P.x > 9][LOG10P.y > LOG10P.x]) / nrow(mergedData2[LOG10P.x > 9])
100 * nrow(mergedData2[LOG10P.x > 9][LOG10P.y < 9]) / nrow(mergedData2[LOG10P.x > 9])

test = cor.test(mergedData2$LOG10P.x,mergedData2$LOG10P.y, method = "spearman")
text = paste("Spearman R:",round(test$estimate,2), "P-value",format.pval(test$p.value,2))
ggplot(data = mergedData2[LOG10P.x > 1], aes(LOG10P.x, LOG10P.y)) + 
  geom_point() +
  # geom_smooth(method = "lm") + 
  geom_abline(slope = 1, intercept = 0, color = "grey", linetype = "dashed", size = 2) +
  annotate("text", x = Inf, y = Inf, label = text, hjust = 1.5, vjust = 1.5, size = 7, fontface = "bold"  ) +
  ylab("Conditional (-log10)") +
  xlab("No conditional (-log10)") +
  theme_minimal() + 
  theme(text = element_text(size=32), legend.position = c(0.86,0.85),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=3), aspect.ratio = 1)  

dev.off()
