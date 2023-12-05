# Plot number of gene-trait associations depending on sample size

library(data.table)
library(ggplot2)

png("~/git/noncoding_rvat/figures/png",720,650)

## P1e9
nc_fdr5_30k = fread("~/burdenNC/data/results/137k_30k/NC_30k_SKATO_process_PVAL1e9.tsv")
nc_fdr5_83k = fread("~/burdenNC/data/results/83k/NC_83k_SKATO_process_PVAL1e9.tsv")
nc_fdr5_137k = fread("~/burdenNC/data/results/137k_30k/NC_137k_SKATO_process_PVAL1e9.tsv")
nc_fdr5_200k = fread("~/burdenNC/data/results/200k/NC_SKATO_process_PVAL1e9.tsv")
 
cds_fdr5_30k = fread("~/burdenNC/data/results/137k_30k/CDS_30k_SKATO_process_PVAL1e9.tsv")
cds_fdr5_83k = fread("~/burdenNC/data/results/83k/CDS_83k_SKATO_process_PVAL1e9.tsv")
cds_fdr5_137k = fread("~/burdenNC/data/results/137k_30k/CDS_137k_SKATO_process_PVAL1e9.tsv")
cds_fdr5_200k = fread("~/burdenNC/data/results/200k/CDS_SKATO_process_PVAL1e9.tsv")


results = c(nrow(nc_fdr5_30k),nrow(nc_fdr5_83k),nrow(nc_fdr5_137k),nrow(nc_fdr5_200k),nrow(cds_fdr5_30k),nrow(cds_fdr5_83k),nrow(cds_fdr5_137k),nrow(cds_fdr5_200k))

data = data.table(sample_size = c(26212,83370,136740,166740,26212,83370,136740,166740), annotation = c("NC","NC","NC","NC","CDS","CDS","CDS","CDS"), results)

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
ggplot(data = data, aes(x = sample_size, y = results, color = annotation)) + 
  geom_abline(intercept = coef(modelNC)[1], slope = coef(modelNC)[2], color = "#D95F02", size = 1.5) +
  geom_abline(intercept = coef(modelCDS)[1], slope = coef(modelCDS)[2], color = "#1B9E77", size = 1.5) +
  geom_point(size = 3) +
  geom_point(data = data[sample_size == 400000], aes(x = sample_size, y = results, color = annotation), color = "black", fill = "black", shape = 21, size = 3) +
  scale_color_manual(values = c("#1B9E77","#D95F02")) +
  ylab("gene-trait associations") +
  xlab("Sample size") +
  theme_minimal() + 
  theme(text = element_text(size=32), legend.position = c(0.36,0.85),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=3), aspect.ratio = 1)  

dev.off()
