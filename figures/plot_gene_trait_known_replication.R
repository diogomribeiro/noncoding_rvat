

library(data.table)
library(ggplot2)

data = fread("~/burdenNC/data/known_replication/genebass_gwas_full.summary", header = F)
# data = fread("~/burdenNC/data/known_replication/gwas_full.summary", header = F)
# data = fread("~/burdenNC/data/known_replication/genebass_full.summary", header = F)
colnames(data) = c("replication","percentage","V3")
data$annotation = data.table(unlist(lapply(data$V3, function(x) unlist(strsplit(x,"[.]"))[1])))$V1
data$pval_cutoff = data.table(unlist(lapply(data$V3, function(x) unlist(strsplit(x,"[.]"))[2])))$V1

data[annotation == "NC_shuffle"]$annotation = "EXPEC."
data[annotation == "JAVIERRE"]$annotation = "HIC"

png("~/known_direct_fdr0.05.png",1200,800)
ggplot(data[replication == "direct_replication"][pval_cutoff == "F0"], aes(x=annotation, y=percentage, fill = annotation )) + 
  geom_bar( stat = "identity", position = "dodge", size = 1, color = "black") +
  geom_text(aes(label = paste0(round(percentage,1),"%")), size = 7, fontface= "bold", vjust = -0.5) +
  annotate(geom = "text", x = "CRD", y = -1.8, label = "N=165|33", size = 8) +
  annotate(geom = "text", x = "ABC", y = -1.8, label = "N=450|131", size = 8) +
  annotate(geom = "text", x = "HIC", y = -1.8, label = "N=1074|335", size = 8) +
  annotate(geom = "text", x = "NC", y = -1.8, label = "N=1513|399", size = 8) +
  annotate(geom = "text", x = "CDS", y = -1.8, label = "N=1159|311", size = 8) +
  annotate(geom = "text", x = "CONTROL", y = -1.8, label = "N=352|104", size = 8) +
  annotate(geom = "text", x = "EXPEC.", y = -1.8, label = "N=1513|399", size = 8) +
  scale_x_discrete(limits = c("ABC","CRD","HIC","NC","CDS","CONTROL","EXPEC.")) +
  labs(x="Annotation", y = "% known gene-trait") +
  guides(fill = "none") +
  ylim(c(0,100)) +
  scale_fill_manual(values = c("#238b45","#6a51a3","#cb181d","#2171b5","#525252","#d9d9d9","#fed976"), limits = c("ABC","CRD","HIC","NC","CDS","CONTROL","EXPEC.")) +
  theme_linedraw() +
  theme(text = element_text(size=40), legend.position = c(0.40,0.85),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=3))
dev.off()

png("~/known_direct_p1e9.png",1200,800)
ggplot(data[replication == "direct_replication"][pval_cutoff == "P1e9"], aes(x=annotation, y=percentage, fill = annotation )) + 
  geom_bar( stat = "identity", position = "dodge", size = 1, color = "black") +
  geom_text(aes(label = paste0(round(percentage,1),"%")), size = 7, fontface= "bold", vjust = -0.5) +
  annotate(geom = "text", x = "CRD", y = -1.8, label = "N=165|33", size = 8) +
  annotate(geom = "text", x = "ABC", y = -1.8, label = "N=450|131", size = 8) +
  annotate(geom = "text", x = "HIC", y = -1.8, label = "N=1074|335", size = 8) +
  annotate(geom = "text", x = "NC", y = -1.8, label = "N=1513|399", size = 8) +
  annotate(geom = "text", x = "CDS", y = -1.8, label = "N=1159|311", size = 8) +
  annotate(geom = "text", x = "CONTROL", y = -1.8, label = "N=352|104", size = 8) +
  annotate(geom = "text", x = "EXPEC.", y = -1.8, label = "N=1513|399", size = 8) +
  scale_x_discrete(limits = c("ABC","CRD","HIC","NC","CDS","CONTROL","EXPEC.")) +
  labs(x="Annotation", y = "% known gene-trait") +
  guides(fill = "none") +
  ylim(c(0,100)) +
  scale_fill_manual(values = c("#238b45","#6a51a3","#cb181d","#2171b5","#525252","#d9d9d9","#fed976"), limits = c("ABC","CRD","HIC","NC","CDS","CONTROL","EXPEC.")) +
  theme_linedraw() +
  theme(text = element_text(size=40), legend.position = c(0.40,0.85),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=3))
dev.off()


png("~/known_window500kb_fdr0.05.png",1200,800)
ggplot(data[replication == "window_replication"][pval_cutoff == "F0"], aes(x=annotation, y=percentage, fill = annotation )) + 
  geom_bar( stat = "identity", position = "dodge", size = 1, color = "black") +
  geom_text(aes(label = paste0(round(percentage,1),"%")), size = 7, fontface= "bold", vjust = -0.5) +
  annotate(geom = "text", x = "CRD", y = -1.8, label = "N=165|33", size = 8) +
  annotate(geom = "text", x = "ABC", y = -1.8, label = "N=450|131", size = 8) +
  annotate(geom = "text", x = "HIC", y = -1.8, label = "N=1074|335", size = 8) +
  annotate(geom = "text", x = "NC", y = -1.8, label = "N=1513|399", size = 8) +
  annotate(geom = "text", x = "CDS", y = -1.8, label = "N=1159|311", size = 8) +
  annotate(geom = "text", x = "CONTROL", y = -1.8, label = "N=352|104", size = 8) +
  annotate(geom = "text", x = "EXPEC.", y = -1.8, label = "N=1513|399", size = 8) +
  scale_x_discrete(limits = c("ABC","CRD","HIC","NC","CDS","CONTROL","EXPEC.")) +
  labs(x="Annotation", y = "% known gene-trait") +
  guides(fill = "none") +
  ylim(c(0,100)) +
  scale_fill_manual(values = c("#238b45","#6a51a3","#cb181d","#2171b5","#525252","#d9d9d9","#fed976"), limits = c("ABC","CRD","HIC","NC","CDS","CONTROL","EXPEC.")) +
  theme_linedraw() +
  theme(text = element_text(size=40), legend.position = c(0.40,0.85),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=3))
dev.off()


png("~/known_window500kb_p1e9.png",1200,800)
ggplot(data[replication == "window_replication"][pval_cutoff == "P1e9"], aes(x=annotation, y=percentage, fill = annotation )) + 
  geom_bar( stat = "identity", position = "dodge", size = 1, color = "black") +
  geom_text(aes(label = paste0(round(percentage,1),"%")), size = 7, fontface= "bold", vjust = -0.5) +
  annotate(geom = "text", x = "CRD", y = -1.8, label = "N=165|33", size = 8) +
  annotate(geom = "text", x = "ABC", y = -1.8, label = "N=450|131", size = 8) +
  annotate(geom = "text", x = "HIC", y = -1.8, label = "N=1074|335", size = 8) +
  annotate(geom = "text", x = "NC", y = -1.8, label = "N=1513|399", size = 8) +
  annotate(geom = "text", x = "CDS", y = -1.8, label = "N=1159|311", size = 8) +
  annotate(geom = "text", x = "CONTROL", y = -1.8, label = "N=352|104", size = 8) +
  annotate(geom = "text", x = "EXPEC.", y = -1.8, label = "N=1513|399", size = 8) +
  scale_x_discrete(limits = c("ABC","CRD","HIC","NC","CDS","CONTROL","EXPEC.")) +
  labs(x="Annotation", y = "% known gene-trait") +
  guides(fill = "none") +
  ylim(c(0,100)) +
  scale_fill_manual(values = c("#238b45","#6a51a3","#cb181d","#2171b5","#525252","#d9d9d9","#fed976"), limits = c("ABC","CRD","HIC","NC","CDS","CONTROL","EXPEC.")) +
  theme_linedraw() +
  theme(text = element_text(size=40), legend.position = c(0.40,0.85),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=3))
dev.off()
