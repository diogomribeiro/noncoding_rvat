
library(data.table)
library(ggplot2)
library(patchwork)

png("~/fig4.png",1800,1700)

##############
# PANEL A
##############

# mergedData = mergedData[LOG10P.y > 1 | LOG10P.x > 1]

mergedData = fread("~/git/noncoding_rvat/figures/data/fig4_conditional_NC.tsv.gz")

# test = cor.test(mergedData$LOG10P.x,mergedData$LOG10P.y, method = "spearman")
# text = paste("Spearman R:",round(test$estimate,2), "P-value",format.pval(test$p.value,2))
g1 = ggplot(data = mergedData, aes(LOG10P.x, LOG10P.y)) + 
  geom_point() +
  geom_point(x = mergedData[GENE == "ENSG00000101439"][TRAIT == "p30720"]$LOG10P.x, y = mergedData[GENE == "ENSG00000101439"][TRAIT == "p30720"]$LOG10P.y, size=14, shape=1, color="#006d2c") + 
  geom_text(x = mergedData[GENE == "ENSG00000101439"][TRAIT == "p30720"]$LOG10P.x + 12, y = mergedData[GENE == "ENSG00000101439"][TRAIT == "p30720"]$LOG10P.y - 3, size=11, label = "CST3 / Cystatin C", color="#006d2c") + 
  # geom_smooth(method = "lm") +
  geom_abline(slope = 1, intercept = 0, color = "grey", linetype = "dashed", size = 2) +
  # annotate("text", x = Inf, y = Inf, label = text, hjust = 1.2, vjust = 1.5, size = 8, fontface = "bold"  ) +
  ylab("Conditional p-value (-log10)") +
  xlab("No conditional p-value (-log10)") +
  labs(tag = "a") + 
  theme_minimal() + 
  theme(text = element_text(size=36), legend.position = c(0.86,0.85),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=3), aspect.ratio = 1)  

##############
# PANEL B
##############

# mergedData = mergedData[LOG10P.y > 1 | LOG10P.x > 1]

mergedData = fread("~/git/noncoding_rvat/figures/data/fig4_conditional_CDS.tsv.gz")

# test = cor.test(mergedData$LOG10P.x,mergedData$LOG10P.y, method = "spearman")
# text = paste("Spearman R:",round(test$estimate,2), "P-value",format.pval(test$p.value,2))
g2 = ggplot(data = mergedData, aes(LOG10P.x, LOG10P.y)) + 
  geom_point() +
  # geom_smooth(method = "lm") +
  geom_abline(slope = 1, intercept = 0, color = "grey", linetype = "dashed", size = 2) +
  # annotate("text", x = Inf, y = Inf, label = text, hjust = 1.2, vjust = 1.5, size = 8, fontface = "bold"  ) +
  ylab("Conditional p-value (-log10)") +
  xlab("No conditional p-value (-log10)") +
  labs(tag = "b") +
  theme_minimal() + 
  theme(text = element_text(size=36), legend.position = c(0.86,0.85),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=3), aspect.ratio = 1)  


##############
# PANEL C
##############


fdr5Data = fread("~/git/noncoding_rvat/figures/data/fig4_single_burden_fdr5.tsv")
fdr5NCSum = sum(fdr5Data[1:3,]$values)
fdr5CDSSum = sum(fdr5Data[4:6,]$values)

g3 = ggplot(data = fdr5Data, aes(x = annotation, y = values, fill = association)) + 
  geom_bar(stat = "identity", color = "black", size = 1.5) +
  annotate(geom = "text", x = "CDS", y = fdr5Data[annotation == "CDS"][association == "Single"]$values, label = paste0(fdr5Data[annotation == "CDS"][association == "Single"]$values," (", round(100 * fdr5Data[annotation == "CDS"][association == "Single"]$values/fdr5CDSSum),"%)"), vjust = 2, size = 9) +
  annotate(geom = "text", x = "CDS", y = fdr5Data[annotation == "CDS"][association == "Single"]$values + 180, label = paste0(fdr5Data[annotation == "CDS"][association == "Group-Single"]$values," (", round(100 * fdr5Data[annotation == "CDS"][association == "Group-Single"]$values/fdr5CDSSum),"%)"), vjust = 2, size = 9) +
  annotate(geom = "text", x = "CDS", y = fdr5Data[annotation == "CDS"][association == "Single"]$values + 500, label = paste0(fdr5Data[annotation == "CDS"][association == "Group"]$values," (", round(100 * fdr5Data[annotation == "CDS"][association == "Group"]$values/fdr5CDSSum),"%)"), vjust = 2, size = 9) +
  annotate(geom = "text", x = "NC", y = fdr5Data[annotation == "NC"][association == "Single"]$values, label = paste0(fdr5Data[annotation == "NC"][association == "Single"]$values," (", round(100 * fdr5Data[annotation == "NC"][association == "Single"]$values/fdr5NCSum),"%)"), vjust = 2, size = 9) +
  annotate(geom = "text", x = "NC", y = fdr5Data[annotation == "NC"][association == "Single"]$values + 100, label = paste0(fdr5Data[annotation == "NC"][association == "Group-Single"]$values," (", round(100 * fdr5Data[annotation == "NC"][association == "Group-Single"]$values/fdr5NCSum),"%)"), vjust = 2, size = 9) +
  annotate(geom = "text", x = "NC", y = fdr5Data[annotation == "NC"][association == "Single"]$values + 350, label = paste0(fdr5Data[annotation == "NC"][association == "Group"]$values," (", round(100 * fdr5Data[annotation == "NC"][association == "Group"]$values/fdr5NCSum),"%)"), vjust = 2, size = 9) +
  # geom_text(aes(label = values)) +
  scale_fill_manual(values = c("#238b45","#74c476","#fc9272")) +
  ylab("gene-trait associations") +
  xlab("Annotation") +
  labs(tag = "c") + 
  theme_minimal() + 
  theme(text = element_text(size=36), legend.position = c(0.3,0.85),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=3), aspect.ratio = 1)  



##############
# PANEL D
##############

p1e9Data = fread("~/git/noncoding_rvat/figures/data/fig4_single_burden_p1e9.tsv")
p1e9NCSum = sum(p1e9Data[1:3,]$values)
p1e9CDSSum = sum(p1e9Data[4:6,]$values)

g4 = ggplot(data = p1e9Data, aes(x = annotation, y = values, fill = association)) + 
  geom_bar(stat = "identity", color = "black", size = 1.5) +
  annotate(geom = "text", x = "CDS", y = p1e9Data[annotation == "CDS"][association == "Single"]$values, label = paste0(p1e9Data[annotation == "CDS"][association == "Single"]$values," (", round(100 * p1e9Data[annotation == "CDS"][association == "Single"]$values/p1e9CDSSum),"%)"), vjust = 3, size = 9) +
  annotate(geom = "text", x = "CDS", y = p1e9Data[annotation == "CDS"][association == "Single"]$values + 150, label = paste0(p1e9Data[annotation == "CDS"][association == "Group-Single"]$values," (", round(100 * p1e9Data[annotation == "CDS"][association == "Group-Single"]$values/p1e9CDSSum),"%)"), vjust = 2, size = 9) +
  annotate(geom = "text", x = "CDS", y = p1e9Data[annotation == "CDS"][association == "Single"]$values + 200, label = paste0(p1e9Data[annotation == "CDS"][association == "Group"]$values," (", round(100 * p1e9Data[annotation == "CDS"][association == "Group"]$values/p1e9CDSSum),"%)"), vjust = 2, size = 9) +
  annotate(geom = "text", x = "NC", y = p1e9Data[annotation == "NC"][association == "Single"]$values, label = paste0(p1e9Data[annotation == "NC"][association == "Single"]$values," (", round(100 * p1e9Data[annotation == "NC"][association == "Single"]$values/p1e9NCSum),"%)"), vjust = 3, size = 9) +
  annotate(geom = "text", x = "NC", y = p1e9Data[annotation == "NC"][association == "Single"]$values + 80, label = paste0(p1e9Data[annotation == "NC"][association == "Group-Single"]$values," (", round(100 * p1e9Data[annotation == "NC"][association == "Group-Single"]$values/p1e9NCSum),"%)"), vjust = 2, size = 9) +
  annotate(geom = "text", x = "NC", y = p1e9Data[annotation == "NC"][association == "Single"]$values + 120, label = paste0(p1e9Data[annotation == "NC"][association == "Group"]$values," (", round(100 * p1e9Data[annotation == "NC"][association == "Group"]$values/p1e9NCSum),"%)"), vjust = 2, size = 9) +
  # geom_text(aes(label = values)) +
  scale_fill_manual(values = c("#238b45","#74c476","#fc9272")) +
  ylab("gene-trait associations") +
  xlab("Annotation") +
  labs(tag = "d") + 
  theme_minimal() + 
  theme(text = element_text(size=36), legend.position = c(0.3,0.85),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=3), aspect.ratio = 1)  



g1 + g2 + g3 + g4
dev.off()

