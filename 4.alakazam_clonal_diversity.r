library(shazam)
library(alakazam)
library(Cairo)
library(dplyr)
library(cowplot)
library("stringr")


setwd("/home/zyfniu//analysis/IgA_IgG/")
clones <- countClones(all, groups="SAMPLE")
head(clones, 5)

head(all)
# Partitions the data based on both the SAMPLE and ISOTYPE columns
# Weights the clone sizes by the DUPCOUNT column
clones <- countClones(all,groups=c("SAMPLE", "CREGION"), copy="DUPCOUNT")

##IgA Specific 
clones <- countClones(all %>% filter(CREGION %in% c("Mouse-IGHA")), groups=c("SAMPLE", "CREGION"), copy="DUPCOUNT")
head(clones, 5)

#g <- ggplot(as.data.frame(clones), aes(SAMPLE)) + geom_bar(aes(fill="CLONE"), width = 0.5)
clones$SAMPLE <- factor(clones$SAMPLE, levels = c("MM86_WT_NAIVE","MM87_KO_NAIVE","MM85-M3_WT","MM85-M4_WT","MM86-F1_WT","MM86-M6_KO","MM89-F1_KO","MM89-M3_KO"))

g <- ggplot(clones, aes(x = SAMPLE, y= SEQ_FREQ, fill= factor(CLONE))) + geom_bar(width = 0.5, stat="identity") + theme(legend.position = "none")

CairoPNG(filename = "IgA_clonal_abundance.png", width = 7500, height = 5000, pointsize = 12,
         fallback_resolution = 600,res = 600)
g 
dev.off()

# ggplot(clones, aes(SAMPLE)) + geom_bar(aes(fill=CLONE))

# Partitions the data on the SAMPLE column
# Calculates a 95% confidence interval via 200 bootstrap realizations



clones <- estimateAbundance(all, group="SAMPLE", ci=0.95, nboot=200, copy="DUPCOUNT")

# Plots a rank abundance curve of the relative clonal abundances
#sample_colors <- c("-1h"="seagreen", "+7d"="steelblue")
CairoPNG(filename = "clonal_abundance.png",
         width = 5000, height = 5000, pointsize = 12,
         fallback_resolution = 600,res = 600)
plot(clones, legend_title="Sample")
dev.off()

# Compare diversity curve across values in the "SAMPLE" column
# q ranges from 0 (min_q=0) to 4 (max_q=4) in 0.05 incriments (step_q=0.05)
# A 95% confidence interval will be calculated (ci=0.95)
# 200 resampling realizations are performed (nboot=200)
sample_div <- rarefyDiversity(all, "SAMPLE", min_q=0, max_q=4, step_q=0.05,
                              ci=0.95, nboot=200)

# Compare diversity curve across values in the "ISOTYPE" column
# Analyse is restricted to ISOTYPE values with at least 30 sequences by min_n=30
# Excluded groups are indicated by a warning message
isotype_div <- rarefyDiversity(all, "CREGION", min_n=30, min_q=0, max_q=4, 
                               step_q=0.05, ci=0.95, nboot=200)

# Plot a log-log (log_q=TRUE, log_d=TRUE) plot of sample diversity
# Indicate number of sequences resampled from each group in the title
sample_main <- paste0("Sample diversity (n=", sample_div@n, ")")
#sample_colors <- c("-1h"="seagreen", "+7d"="steelblue")
CairoPNG(filename = "sample_diversity.png",
         width = 5000, height = 5000, pointsize = 12,
         fallback_resolution = 600,res = 600)
plot(sample_div, main_title=sample_main, 
     legend_title="Sample")
dev.off()


CairoPNG(filename = "isotype_diversity.png",
         width = 5000, height = 5000, pointsize = 12,
         fallback_resolution = 600,res = 600)
isotype_main <- paste0("Isotype diversity (n=", isotype_div@n, ")")
plot(isotype_div, main_title=isotype_main, 
     legend_title="Isotype")
dev.off()

sample_test <- testDiversity(all, 2, "SAMPLE", nboot=200)
CairoPNG(filename = "sample_diversity_fixed_order.png",
         width = 5000, height = 5000, pointsize = 12,
         fallback_resolution = 600,res = 600)
plot(sample_test, main_title=sample_main, 
     legend_title="Sample Diversity")
dev.off()


##IgA Only
#
#
#
#

all_IgA <- all %>% filter(CREGION %in% c("Mouse-IGHA"))
all_IgG <- all %>% filter(CREGION %in% c("Mouse-IGHG12","IGHG3"))

clones_IgA <- estimateAbundance(all_IgA, group="SAMPLE", ci=0.95, nboot=200, copy="DUPCOUNT")
clones_IgA$SAMPLE <- factor(clones_IgA$SAMPLE, levels = c("MM86_WT_NAIVE","MM87_KO_NAIVE","MM85-M3_WT","MM85-M4_WT","MM86-F1_WT","MM86-M6_KO","MM89-F1_KO","MM89-M3_KO"))
clones_IgG <- estimateAbundance(all_IgG, group="SAMPLE", ci=0.95, nboot=200, copy="DUPCOUNT")
clones_IgG$SAMPLE <- factor(clones_IgG$SAMPLE, levels = c("MM86_WT_NAIVE","MM87_KO_NAIVE","MM85-M3_WT","MM85-M4_WT","MM86-F1_WT","MM86-M6_KO","MM89-F1_KO","MM89-M3_KO"))

# Plots a rank abundance curve of the relative clonal abundances
#sample_colors <- c("-1h"="seagreen", "+7d"="steelblue")
setwd("../IgA_IgG/")
typeof(clones_IgA)
CairoPNG(filename = "clonal_abundance.png",
         width = 10000, height = 5000, pointsize = 12,
         fallback_resolution = 600,res = 600)
p1 <- plot(clones_IgG, legend_title="IgG")
p2 <- plot(clones_IgA, legend_title="IgA")
plot_grid(p1,p2,labels = c('IgG', 'IgA'))
dev.off()

# Compare diversity curve across values in the "SAMPLE" column
# q ranges from 0 (min_q=0) to 4 (max_q=4) in 0.05 incriments (step_q=0.05)
# A 95% confidence interval will be calculated (ci=0.95)
# 200 resampling realizations are performed (nboot=200)
sample_div_IgA <- rarefyDiversity(all_IgA, "SAMPLE", min_q=0, max_q=4, step_q=0.05,
                              ci=0.95, nboot=200)


all_IgG

# Compare diversity curve across values in the "ISOTYPE" column
# Analyse is restricted to ISOTYPE values with at least 30 sequences by min_n=30
# Excluded groups are indicated by a warning message
isotype_div <- rarefyDiversity(all_IgA, "CREGION", min_n=30, min_q=0, max_q=4, 
                               step_q=0.05, ci=0.95, nboot=200)
isotype_div$SAMPLE <- factor(isotype_div$SAMPLE, levels = c("MM86_WT_NAIVE","MM87_KO_NAIVE","MM85-M3_WT","MM85-M4_WT","MM86-F1_WT","MM86-M6_KO","MM89-F1_KO","MM89-M3_KO"))

# Plot a log-log (log_q=TRUE, log_d=TRUE) plot of sample diversity
# Indicate number of sequences resampled from each group in the title
sample_main <- paste0("Sample diversity (n=", sample_div@n, ")")
#sample_colors <- c("-1h"="seagreen", "+7d"="steelblue")
CairoPNG(filename = "sample_diversity.png",
         width = 5000, height = 5000, pointsize = 12,
         fallback_resolution = 600,res = 600)
plot(sample_div, main_title=sample_main, 
     legend_title="Sample")
dev.off()


CairoPNG(filename = "isotype_diversity.png",
         width = 5000, height = 5000, pointsize = 12,
         fallback_resolution = 600,res = 600)
isotype_main <- paste0("Isotype diversity (n=", isotype_div@n, ")")
plot(isotype_div, main_title=isotype_main, 
     legend_title="Isotype")
dev.off()
head(all_IgA)
sample_test_IgA <- alphaDiversity(all_IgA, group=c("SAMPLE","GROUP"), min_q=0, max_q=2, step_q=1, nboot=200)

head(all)

sample_test_IgG <- alphaDiversity(all, group=c("SAMPLE"), min_q=0, max_q=2, step_q=1, nboot=200)
sample_test_IgG
plot_temp_IgA <- as.data.frame((sample_test_IgA@diversity))
plot_temp_IgA
plot_temp_IgA[,1] <- factor(plot_temp_IgA[,1], levels = c("MM86_WT_NAIVE","MM87_KO_NAIVE", "MM85-M3_WT","MM85-M4_WT","MM86-F1_WT","MM86-M6_KO","MM89-F1_KO","MM89-M3_KO"))
plot_temp_IgG <- as.data.frame((sample_test_IgG@diversity))
plot_temp_IgG[,1] <- factor(plot_temp_IgG[,1], levels = positions)
plot_temp_IgG

write.csv(plot_temp_IgG,"./hill_number_igg.csv")
write.csv(plot_temp_IgA,"./hill_number_iga.csv")
CairoPNG(filename = "sample_diversity_fixed_order.png",
         width = 20000, height = 10000, pointsize = 12,
         fallback_resolution = 600,res = 600)
g3 <-ggplot(plot_temp_IgG, aes(x=SAMPLE, y=D,color=SAMPLE)) + 
        theme_bw() + ggtitle("Hills Diversity") +
        xlab("Sample") + ylab("Mean QD +_") +
 #       scale_fill_brewer(palette="Dark2") +
        geom_boxplot() +scale_color_manual(values=c("black","black", "black",
                                                    "blue", "blue", "blue", "blue", "blue", "blue", "blue","blue",
                                                    "green","green","green",
                                                    "red","red","red","red","red","red","red","red","red"
                                                    )) +
        geom_errorbar(aes(ymin=D_LOWER, ymax=D_UPPER), width=.2,
                      position=position_dodge(.9)) + facet_grid(. ~ Q)
plot_grid(g3, ncol = 1)
dev.off()
#setwd("../")

?rcolorbrewer

grand.mean <- function(M, N) {weighted.mean(M, N)}
grand.sd   <- function(S, M, N) {sqrt(weighted.mean(S^2 + M^2, N) -
                                              weighted.mean(M, N)^2)}
IgG

library(dplyr)
as.vector(plot_temp_IgA %>%
        filter(str_detect(SAMPLE, 'WT')) %>%
        filter(!str_detect(SAMPLE, 'NAIVE')) %>%
        filter(Q == 0))$D
        

grand.mean((as.vector(plot_temp_IgA %>%
                              filter(str_detect(SAMPLE, 'WT')) %>%
                              filter(!str_detect(SAMPLE, 'NAIVE')) %>%
                              filter(Q == 0))$D3))
