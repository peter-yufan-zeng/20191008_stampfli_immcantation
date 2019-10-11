# Import required packages
library(alakazam)
library(shazam)
library(dplyr)
library(ggplot2)
library(Cairo)
library(cowplot)
library(RColorBrewer)
library(colorRamps)

setwd("C:/Users/zyfni/Desktop/20191006_Stampfli_immcantation_analysis/20191008_stampfli_immcantation/")

ra1 <- read.table("../data/ra_ctrl1_db-pass_clone-pass_germ-pass.tab", header = T, sep = "\t")
ra1$SAMPLE <- "RA_CTRL1"
ra2 <- read.table("../data/ra_ctrl2_db-pass_clone-pass_germ-pass.tab", header = T, sep = "\t")
ra2$SAMPLE <- "RA_CTRL2"
ra3 <- read.table("../data/ra_ctrl3_db-pass_clone-pass_germ-pass.tab", header = T, sep = "\t")
ra3$SAMPLE <- "RA_CTRL3"
ra4 <- read.table("../data/ra_imm4_db-pass_clone-pass_germ-pass.tab", header = T, sep = "\t")
ra4$SAMPLE <- "RA_IMM4"
ra5 <- read.table("../data/ra_imm5_db-pass_clone-pass_germ-pass.tab", header = T, sep = "\t")
ra5$SAMPLE <- "RA_IMM5"
ra6 <- read.table("../data/ra_imm6_db-pass_clone-pass_germ-pass.tab", header = T, sep = "\t")
ra6$SAMPLE <- "RA_IMM6"
ra7 <- read.table("../data/ra_imm7_db-pass_clone-pass_germ-pass.tab", header = T, sep = "\t")
ra7$SAMPLE <- "RA_IMM7"
ra8 <- read.table("../data/ra_imm8_db-pass_clone-pass_germ-pass.tab", header = T, sep = "\t")
ra8$SAMPLE <- "RA_IMM8"
ra9 <- read.table("../data/ra_imm9_db-pass_clone-pass_germ-pass.tab", header = T, sep = "\t")
ra9$SAMPLE <- "RA_IMM9"
ra11 <- read.table("../data/ra_imm11_db-pass_clone-pass_germ-pass.tab", header = T, sep = "\t")
ra11$SAMPLE <- "RA_IMM11"
ra12 <- read.table("../data/ra_imm12_db-pass_clone-pass_germ-pass.tab", header = T, sep = "\t")
ra12$SAMPLE <- "RA_IMM12"
cs14 <- read.table("../data/cs_ctrl14_db-pass_clone-pass_germ-pass.tab", header = T, sep = "\t")
cs14$SAMPLE <- "CS_CTRL14"
cs15 <- read.table("../data/cs_ctrl15_db-pass_clone-pass_germ-pass.tab", header = T, sep = "\t")
cs15$SAMPLE <- "CS_CTRL15"
cs16 <- read.table("../data/cs_ctrl16_db-pass_clone-pass_germ-pass.tab", header = T, sep = "\t")
cs16$SAMPLE <- "CS_CTRL16"
cs17 <- read.table("../data/cs_imm17_db-pass_clone-pass_germ-pass.tab", header = T, sep = "\t")
cs17$SAMPLE <- "CS_IMM17"
cs18 <- read.table("../data/cs_imm18_db-pass_clone-pass_germ-pass.tab", header = T, sep = "\t")
cs18$SAMPLE <- "CS_IMM18"
cs19 <- read.table("../data/cs_imm19_db-pass_clone-pass_germ-pass.tab", header = T, sep = "\t")
cs19$SAMPLE <- "CS_IMM19"
cs21 <- read.table("../data/cs_imm21_db-pass_clone-pass_germ-pass.tab", header = T, sep = "\t")
cs21$SAMPLE <- "CS_IMM21"
cs22 <- read.table("../data/cs_imm22_db-pass_clone-pass_germ-pass.tab", header = T, sep = "\t")
cs22$SAMPLE <- "CS_IMM22"
cs23 <- read.table("../data/cs_imm23_db-pass_clone-pass_germ-pass.tab", header = T, sep = "\t")
cs23$SAMPLE <- "CS_IMM23"
cs24 <- read.table("../data/cs_imm24_db-pass_clone-pass_germ-pass.tab", header = T, sep = "\t")
cs24$SAMPLE <- "CS_IMM24"
cs25 <- read.table("../data/cs_imm25_db-pass_clone-pass_germ-pass.tab", header = T, sep = "\t")
cs25$SAMPLE <- "CS_IMM25"
cs26 <- read.table("../data/cs_imm26_db-pass_clone-pass_germ-pass.tab", header = T, sep = "\t")
cs26$SAMPLE <- "CS_IMM26"



ra1$GROUP <- "RA_CTRL"
ra2$GROUP <- "RA_CTRL"
ra3$GROUP <- "RA_CTRL"
ra4$GROUP <- "RA_IMM"
ra5$GROUP <- "RA_IMM"
ra6$GROUP <- "RA_IMM"
ra7$GROUP <- "RA_IMM"
ra8$GROUP <- "RA_IMM"
ra9$GROUP <- "RA_IMM"
ra11$GROUP <- "RA_IMM"
ra12$GROUP <- "RA_IMM"
cs14$GROUP <- "CS_CTRL"
cs15$GROUP <- "CS_CTRL"
cs16$GROUP <- "CS_CTRL"
cs17$GROUP <- "CS_IMM"
cs18$GROUP <- "CS_IMM"
cs19$GROUP <- "CS_IMM"
cs21$GROUP <- "CS_IMM"
cs22$GROUP <- "CS_IMM"
cs23$GROUP <- "CS_IMM"
cs24$GROUP <- "CS_IMM"
cs25$GROUP <- "CS_IMM"
cs26$GROUP <- "CS_IMM"
all <- rbind(ra1, ra2, ra3, ra4, ra5, ra6, ra7, ra8, ra9, ra12, ra11, cs14, cs15, cs16,
             cs17, cs18, cs19, cs21, cs22, cs23, cs24, cs25, cs26)


Isotype_proportions <- all %>% select("CREGION", "SAMPLE") %>% filter(CREGION != "")
positions <- c("RA_CTRL1", "RA_CTRL2",  "RA_CTRL3", "RA_IMM4", "RA_IMM5", "RA_IMM6","RA_IMM7", "RA_IMM8",
               "RA_IMM9","RA_IMM11","RA_IMM12","CS_CTRL14","CS_CTRL15","CS_CTRL16","CS_IMM17",
               "CS_IMM18","CS_IMM19","CS_IMM21","CS_IMM22","CS_IMM23","CS_IMM24","CS_IMM25","CS_IMM26")

g <- ggplot(Isotype_proportions, aes(x = SAMPLE)) + geom_bar(aes(fill=CREGION), width = 0.5) 
# +
#   theme(axis.text.x = element_text(angle=65, vjust=0.6)) + ylab("Number of Individual Clones") +
#   labs(title="Percentage of Isotype") + scale_x_discrete(limits = positions)
save_plot("./isotype_proportions.png", g, ncol=3, base_height = 10, base_width = 15)

all_gene <- countGenes(all, gene="V_CALL", groups=c("SAMPLE", "GROUP"),
                       mode="gene", copy="DUPCOUNT")
all_gene
# all_gene <- countGenes(all, gene="V_CALL", groups=c("SAMPLE"), 
#                        mode="gene", copy="DUPCOUNT")

head(all_gene)

write.csv(all_gene, file = "./all_gene.csv")
# Assign sorted levels and subset to IGHV1
# ighv1 <- gene %>%
#   mutate(GENE=factor(GENE, levels=sortGenes(unique(GENE), method="name"))) %>%
#   filter(getFamily(GENE) == "IGHV1")

# Plot V gene usage in the IGHV1 family by sample
g1 <- ggplot(all_gene, aes(x=GENE, y=SEQ_FREQ)) +
  theme_bw() +
  ggtitle("Gene Usage") +
  theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1)) +
  ylab("Percent of repertoire") +
  xlab("") +
  scale_color_brewer(palette="Set1") +
  geom_point(aes(color=GROUP), size=5, alpha=0.8) # + facet_grid(.~CREGION)
g3 <- ggplot(all_gene, aes(x=GENE, y=COPY_FREQ)) +
  theme_bw() +
  ggtitle("Copy Number") +
  theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1)) +
  ylab("Percent of repertoire") +
  xlab("") +
  scale_color_brewer(palette="Set1") +
  geom_point(aes(color=GROUP), size=5, alpha=0.8)# + facet_grid(. ~ CREGION)

g1

g_combined <- plot_grid(g1, g3,  labels = c("A", "B"))
g_combined
save_plot("./20190527_gene_seq_freq.png", g1, ncol=1, base_height = 10, base_width = 30, limitsize = FALSE)
save_plot("./20190527_gene_copy_freq.png", g3, ncol=1, base_height = 10, base_width = 30, limitsize = FALSE)
# CairoPNG(filename = "gene_usage.png",
#          width = 15000, height = 7000, pointsize = 12,
#          fallback_resolution = 600,res = 600)
# plot(g3)
# dev.off()


##ploting only IgM, and IgG
getPalette <- colorRampPalette(brewer.pal("Set1"))

colourCount <- length(unique(all$SAMPLE))
colourCount

n <- 23
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
col_vector
g3 <- ggplot(filter(all_gene, CREGION %in% c("Mouse-IGHA")), aes(x=GENE, y=COPY_FREQ)) +
  ggtitle("Copy Number") +
  theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1)) +
  ylab("Percent of repertoire") +
  xlab("") + scale_color_manual(values=col_vector)+ 
#  scale_fill_manual(values = colorRampPalette(brewer.pal(12, "Accent"))(colourCount)) +
  geom_point(aes(colour=SAMPLE), size=5, alpha=0.8) +
  facet_grid(. ~ CREGION)
g3 
all_gene
CairoPNG(filename = "gene_usage_igm_igg_iga.png",
         width = 20000, height = 7000, pointsize = 12,
         fallback_resolution = 600,res = 600)
plot(g3)
dev.off()



# Calculate V family copy numbers by sample and isotype
family <- countGenes(all, gene="V_CALL", groups=c("SAMPLE","GROUP"), 
                     mode="family", copy="DUPCOUNT")
head(family, n=4)
# Plot V family copy abundance by sample and isotype
g4 <- ggplot(family, aes(x=GENE, y=COPY_FREQ)) +
  theme_bw() +
  ggtitle("Copy Number") +
  theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1)) +
  ylab("Percent of repertoire") +
  xlab("") + scale_color_manual(values=col_vector)+
  geom_point(aes(color=GROUP), size=5, alpha=0.8) 

# + facet_grid(. ~ CREGION)
CairoPNG(filename = "family_usage.png",
         width = 10000, height = 7000, pointsize = 12,
         fallback_resolution = 600,res = 600)
plot(g4)
dev.off()
