
library(data.table)
library(ggplot2)
library(patchwork)


png("~/fig3.png",2400,1700)

##############
# PANEL A
##############

data = fread("~/git/noncoding_rvat/figures/data/fig3_subsample_replication.tsv")

g1 = ggplot(data, aes(x=annotation, y=percentage, alpha = cutoff, fill = annotation )) +
  geom_bar( stat = "identity", position = "dodge", size = 1, color = "black") +
  geom_text(aes(label = paste0(round(percentage,1),"%"), group = cutoff), position = position_dodge(width = .9), vjust = -0.5, alpha = 1, size = 8, fontface = "bold" ) +
  annotate(geom = "text", x = "NC", y = -3, label = paste0("N=",nrow(mergedNC[FDR.y < 0.05]), " N=",nrow(mergedNC[PVAL.y < 1e-9])), size = 8) +
  annotate(geom = "text", x = "CDS", y = -3, label =  paste0("N=",nrow(mergedCDS[FDR.y < 0.05]), " N=",nrow(mergedCDS[PVAL.y < 1e-9])), size = 8) +
  annotate(geom = "text", x = "EXPECTED", y = -3, label = paste0("N=",nrow(mergedSH[FDR.y < 0.05]), " N=",nrow(mergedSH[PVAL.y < 1e-9])), size = 8) +
  scale_x_discrete(limits = c("NC","CDS","EXPECTED")) +
  labs(x="Annotation", y = "% replicated gene-trait", tag = "a") +
  guides(fill = "none") +
  ylim(c(-3, 100)) +
  scale_fill_manual(values = c("#2171b5","#525252","grey"), limits = c("NC","CDS","EXPEC.")) +
  theme_linedraw() +
  theme(text = element_text(size=40), legend.position = c(0.86,0.85),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=3), aspect.ratio = 1)

##############
# PANEL B
##############

data = fread("~/git/noncoding_rvat/figures/data/fig3_sample_size_estimate_fdr5.tsv")

cor.test(data[annotation == "NC"]$results, data[annotation == "NC"]$sample_size)
cor.test(data[annotation == "CDS"]$results, data[annotation == "CDS"]$sample_size)

# NC
modelNC = lm(results ~ sample_size, data = data[annotation == "NC"])
predicted_value <- predict(modelNC, data.table(sample_size = 400000))
data = rbind(data, data.table(sample_size = 400000, annotation = "NC", results = predicted_value))
# CDS
modelCDS = lm(results ~ sample_size, data = data[annotation == "CDS"])
predicted_value <- predict(modelCDS, data.table(sample_size = 400000))
data = rbind(data, data.table(sample_size = 400000, annotation = "CDS", results = predicted_value))

options(scipen = 1)
g2 = ggplot(data = data, aes(x = sample_size, y = results, color = annotation)) + 
  geom_abline(intercept = coef(modelNC)[1], slope = coef(modelNC)[2], color = "#D95F02", size = 1.5) +
  geom_abline(intercept = coef(modelCDS)[1], slope = coef(modelCDS)[2], color = "#1B9E77", size = 1.5) +
  geom_point(size = 3) +
  geom_point(data = data[sample_size == 400000], aes(x = sample_size, y = results, color = annotation), color = "black", fill = "black", shape = 21, size = 3) +
  scale_color_manual(values = c("#1B9E77","#D95F02")) +
  labs(tag = "b") +
  xlim(c(0,410000)) +
  ylab("gene-trait associations") +
  xlab("Sample size") +
  theme_minimal() + 
  theme(text = element_text(size=40), legend.position = c(0.36,0.85),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=3), aspect.ratio = 1)  


##############
# PANEL C
##############

data = fread("~/git/noncoding_rvat/figures/data/fig3_known_gene_trait_direct_fdr5.tsv")

g3 = ggplot(data[replication == "direct_replication"][pval_cutoff == "F0"], aes(x=annotation, y=percentage, fill = annotation )) + 
  geom_bar( stat = "identity", position = "dodge", size = 1, color = "black") +
  geom_text(aes(label = paste0(round(percentage,1),"%")), size = 8, fontface= "bold", vjust = -0.5) +
  annotate(geom = "text", x = "CRD", y = -1.8, label = "N=165|33", size = 8) +
  annotate(geom = "text", x = "ABC", y = -1.8, label = "N=450|131", size = 8) +
  annotate(geom = "text", x = "HIC", y = -1.8, label = "N=1074|335", size = 8) +
  annotate(geom = "text", x = "NC", y = -1.8, label = "N=1513|399", size = 8) +
  annotate(geom = "text", x = "CDS", y = -1.8, label = "N=1159|311", size = 8) +
  annotate(geom = "text", x = "CONTROL", y = -1.8, label = "N=352|104", size = 8) +
  annotate(geom = "text", x = "EXPEC.", y = -1.8, label = "N=1513|399", size = 8) +
  scale_x_discrete(limits = c("ABC","CRD","HIC","NC","CDS","CONTROL","EXPEC.")) +
  labs(x="Annotation", y = "% known gene-trait", tag = "c") +
  guides(fill = "none") +
  ylim(c(0,100)) +
  scale_fill_manual(values = c("#238b45","#6a51a3","#cb181d","#2171b5","#525252","#d9d9d9","#fed976"), limits = c("ABC","CRD","HIC","NC","CDS","CONTROL","EXPEC.")) +
  theme_linedraw() +
  theme(text = element_text(size=40), legend.position = c(0.40,0.85),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=3))


##############
# PANEL D
##############

data = fread("~/git/noncoding_rvat/figures/data/fig3_known_gene_trait_window_fdr5.tsv")

g4 = ggplot(data[replication == "window_replication"][pval_cutoff == "F0"], aes(x=annotation, y=percentage, fill = annotation )) + 
  geom_bar( stat = "identity", position = "dodge", size = 1, color = "black") +
  geom_text(aes(label = paste0(round(percentage,1),"%")), size = 8, fontface= "bold", vjust = -0.5) +
  annotate(geom = "text", x = "CRD", y = -1.8, label = "N=165|33", size = 8) +
  annotate(geom = "text", x = "ABC", y = -1.8, label = "N=450|131", size = 8) +
  annotate(geom = "text", x = "HIC", y = -1.8, label = "N=1074|335", size = 8) +
  annotate(geom = "text", x = "NC", y = -1.8, label = "N=1513|399", size = 8) +
  annotate(geom = "text", x = "CDS", y = -1.8, label = "N=1159|311", size = 8) +
  annotate(geom = "text", x = "CONTROL", y = -1.8, label = "N=352|104", size = 8) +
  annotate(geom = "text", x = "EXPEC.", y = -1.8, label = "N=1513|399", size = 8) +
  scale_x_discrete(limits = c("ABC","CRD","HIC","NC","CDS","CONTROL","EXPEC.")) +
  labs(x="Annotation", y = "% known gene-trait", tag = "d") +
  guides(fill = "none") +
  ylim(c(0,100)) +
  scale_fill_manual(values = c("#238b45","#6a51a3","#cb181d","#2171b5","#525252","#d9d9d9","#fed976"), limits = c("ABC","CRD","HIC","NC","CDS","CONTROL","EXPEC.")) +
  theme_linedraw() +
  theme(text = element_text(size=40), legend.position = c(0.40,0.85),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=3))


g1 + g3 + g2 + g4

dev.off()

