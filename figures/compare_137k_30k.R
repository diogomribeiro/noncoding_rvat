
library(data.table)
library(ggplot2)

cdsData137k = fread("zcat ~/burdenNC/data/results/137k_30k/CDS_137k_SKATO_process_all.tsv.gz")
cdsData30k = fread("zcat  ~/burdenNC/data/results/137k_30k/CDS_30k_SKATO_process_all.tsv.gz")
ncData137k = fread("zcat ~/burdenNC/data/results/137k_30k/NC_137k_SKATO_process_all.tsv.gz")
ncData30k = fread("zcat  ~/burdenNC/data/results/137k_30k/NC_30k_SKATO_process_all.tsv.gz")
shuffle137k = fread("zcat ~/burdenNC/data/results/137k_30k/NC_137k_shuffle_process_all.tsv.gz")

###########
## CDS
###########
png("~/sup7b_CDS_replication.png",800,700)
mergedCDS = merge(cdsData30k,cdsData137k,by=c("GENE","TRAIT"), all.x = T, all.y = T)
paste(nrow(unique(mergedCDS[is.na(PVAL.x)])), "CDS gene-traits not evaluated in 30k")
paste(nrow(unique(mergedCDS[is.na(PVAL.y)])), "CDS gene-traits not evaluated in 137k")
paste(nrow(unique(mergedCDS[!is.na(PVAL.x)][!is.na(PVAL.y)])),"NC gene-traits compared")

mergedCDS = mergedCDS[!is.na(PVAL.x)][!is.na(PVAL.y)]
test = cor.test(-log10(mergedCDS[-log10(PVAL.y) > 3]$PVAL.x),-log10(mergedCDS[-log10(PVAL.y) > 3]$PVAL.y), method = "spearman")
text = paste("Spearman R:",round(test$estimate,2), "P-value",format.pval(test$p.value,2))
ggplot(data = mergedCDS[-log10(PVAL.y) > 3], aes(-log10(PVAL.x), -log10(PVAL.y))) + 
  geom_point() +
  geom_smooth(method = "lm") + 
  geom_abline(slope = 1, intercept = 0, color = "grey", linetype = "dashed", size = 2) +
  annotate("text", x = Inf, y = Inf, label = text, hjust = 1.5, vjust = 1.5, size = 7, fontface = "bold"  ) +
  ylab("-log10(137k)") +
  xlab("-log10(30k)") +
  theme_minimal() + 
  theme(text = element_text(size=32), legend.position = c(0.86,0.85),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=3), aspect.ratio = 1)  

dev.off()

###########
# NC
###########
png("~/sup7a_NC_replication.png",800,700)
mergedNC = merge(ncData30k,ncData137k,by=c("GENE","TRAIT"), all.x = T, all.y = T)
paste(nrow(unique(mergedNC[is.na(PVAL.x)])), "NC gene-traits not evaluated in 30k")
paste(nrow(unique(mergedNC[is.na(PVAL.y)])), "NC gene-traits not evaluated in 137k")
paste(nrow(unique(mergedNC[!is.na(PVAL.x)][!is.na(PVAL.y)])),"NC gene-traits compared")

mergedNC = mergedNC[!is.na(PVAL.x)][!is.na(PVAL.y)]
test = cor.test(-log10(mergedNC[-log10(PVAL.y) > 3]$PVAL.x),-log10(mergedNC[-log10(PVAL.y) > 3]$PVAL.y), method = "spearman")
text = paste("Spearman R:",round(test$estimate,2), "P-value",format.pval(test$p.value,2))
ggplot(data = mergedNC[-log10(PVAL.y) > 3], aes(-log10(PVAL.x), -log10(PVAL.y))) + 
  geom_point() +
  geom_smooth(method = "lm") + 
  geom_abline(slope = 1, intercept = 0, color = "grey", linetype = "dashed", size = 2) +
  annotate("text", x = Inf, y = Inf, label = text, hjust = 1.5, vjust = 1.5, size = 7, fontface = "bold"  ) +
  ylab("-log10(137k)") +
  xlab("-log10(30k)") +
  theme_minimal() + 
  theme(text = element_text(size=32), legend.position = c(0.86,0.85),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=3), aspect.ratio = 1)  

dev.off()


###########
# Shuffled NC
###########
mergedSH = merge(ncData30k,shuffle137k,by=c("GENE","TRAIT"), all.x = T, all.y = T)
paste(nrow(unique(mergedSH[is.na(PVAL.x)])), "NC gene-traits not evaluated in 30k")
paste(nrow(unique(mergedSH[is.na(PVAL.y)])), "NC gene-traits not evaluated in 137k")
paste(nrow(unique(mergedSH[!is.na(PVAL.x)][!is.na(PVAL.y)])),"NC gene-traits compared")
mergedSH = mergedSH[!is.na(PVAL.x)][!is.na(PVAL.y)]

###########
# P<0.05 replication
###########
ncFDR0.05 = nrow(mergedNC[FDR.y < 0.05][PVAL.x < 0.05]) * 100 / nrow(mergedNC[FDR.y < 0.05])
ncP1e9 = nrow(mergedNC[PVAL.y < 1e-9][PVAL.x < 0.05]) * 100 / nrow(mergedNC[PVAL.y < 1e-9])
cdsFDR0.05 = nrow(mergedCDS[FDR.y < 0.05][PVAL.x < 0.05]) * 100 / nrow(mergedCDS[FDR.y < 0.05])
cdsP1e9 = nrow(mergedCDS[PVAL.y < 1e-9][PVAL.x < 0.05]) * 100 / nrow(mergedCDS[PVAL.y < 1e-9])
shFDR0.05 = nrow(mergedSH[FDR.y < 0.05][PVAL.x < 0.05]) * 100 / nrow(mergedSH[FDR.y < 0.05])
shP1e9 = nrow(mergedSH[PVAL.y < 1e-9][PVAL.x < 0.05]) * 100 / nrow(mergedSH[PVAL.y < 1e-9])

data = data.table(annotation = c("NC","NC","CDS", "CDS","EXPECTED","EXPECTED"), 
                  percentage = c(ncFDR0.05,ncP1e9,cdsFDR0.05,cdsP1e9,shFDR0.05,shP1e9), 
                  cutoff = c("FDR 5%","P<1e-9","FDR 5%","P<1e-9","FDR 5%","P<1e-9"))

png("fig3c_replication.png",800,700)
ggplot(data, aes(x=annotation, y=percentage, alpha = cutoff, fill = annotation )) +
  geom_bar( stat = "identity", position = "dodge", size = 1, color = "black") +
  geom_text(aes(label = paste0(round(percentage,1),"%"), group = cutoff), position = position_dodge(width = .9), vjust = -0.5, alpha = 1, size = 9, fontface = "bold" ) +
  annotate(geom = "text", x = "NC", y = -3, label = paste0("N=",nrow(mergedNC[FDR.y < 0.05]), " N=",nrow(mergedNC[PVAL.y < 1e-9])), size = 8) +
  annotate(geom = "text", x = "CDS", y = -3, label =  paste0("N=",nrow(mergedCDS[FDR.y < 0.05]), " N=",nrow(mergedCDS[PVAL.y < 1e-9])), size = 8) +
  annotate(geom = "text", x = "EXPECTED", y = -3, label = paste0("N=",nrow(mergedSH[FDR.y < 0.05]), " N=",nrow(mergedSH[PVAL.y < 1e-9])), size = 8) +
  scale_x_discrete(limits = c("NC","CDS","EXPECTED")) +
  labs(x="Annotation", y = "% replicated gene-trait") +
  guides(fill = "none") +
  ylim(c(-3, 100)) +
  scale_fill_manual(values = c("#2171b5","#525252","grey"), limits = c("NC","CDS","EXPEC.")) +
  theme_linedraw() +
  theme(text = element_text(size=32), legend.position = c(0.86,0.85),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=3), aspect.ratio = 1)


dev.off()

