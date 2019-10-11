library(shazam)

# Calculate R and S mutation counts
db_obs_count <- observedMutations(all, sequenceColumn="SEQUENCE_IMGT",
                                  germlineColumn="GERMLINE_IMGT_D_MASK",
                                  regionDefinition=NULL,
                                  frequency=FALSE, 
                                  nproc=4)
db_obs_count_total <- observedMutations(all, sequenceColumn="SEQUENCE_IMGT",
                                  germlineColumn="GERMLINE_IMGT_D_MASK",
                                  regionDefinition=NULL,
                                  frequency=FALSE,
                                  combine = T,
                                  nproc=4)
# Show new mutation count columns
db_obs_count %>% 
  select(SEQUENCE_ID, starts_with("MU_COUNT_")) %>%
  head(n=4)
db_obs_count
g1 <- ggplot(db_obs_count_total, aes(x=GROUP, y=MU_COUNT, fill=SAMPLE)) +
  theme_bw() + ggtitle("Total mutations") +
  xlab("Group") + ylab("Total Mutation Count")  +
  scale_color_manual(values=col_vector) +
  geom_boxplot()
g2 <- ggplot(db_obs_count, aes(x=GROUP, y=MU_COUNT_SEQ_S, fill=SAMPLE)) +
  theme_bw() + ggtitle("Total silent mutations") +
  xlab("Group") + ylab("Mutation count") +
  scale_color_manual(values=col_vector) +
  geom_boxplot()
g3 <- ggplot(db_obs_count, aes(x=GROUP, y=MU_COUNT_SEQ_R, fill=SAMPLE)) +
  theme_bw() + ggtitle("Total replacement mutations") +
  xlab("Group") + ylab("Mutation count") +
  scale_color_manual(values=col_vector) +
  geom_boxplot()
CairoPNG(filename = "mutation_count_silent_replacement.png",
         width = 15000, height = 5000, pointsize = 12,
         fallback_resolution = 600,res = 600)
alakazam::gridPlot(g1, g2, g3, ncol=3)
dev.off()




# Calculate R and S mutation frequencies
db_obs_freq_total <- observedMutations(all, sequenceColumn="SEQUENCE_IMGT",
                                 germlineColumn="GERMLINE_IMGT_D_MASK",
                                 regionDefinition=NULL,
                                 frequency=TRUE, 
                                 nproc=4)
# Show new mutation frequency columns
db_obs_freq %>% 
  select(SEQUENCE_ID, starts_with("MU_FREQ_")) %>%
  head(n=4)

head(db_obs_freq)

db_obs_freq_total <- observedMutations(all, sequenceColumn="SEQUENCE_IMGT",
                                       germlineColumn="GERMLINE_IMGT_D_MASK",
                                       regionDefinition=NULL,
                                       frequency=TRUE,
                                       combine = T,
                                       nproc=4)

g1 <- ggplot(db_obs_freq_total, aes(x=GROUP, y=MU_FREQ, fill=SAMPLE)) +
   theme_bw() + ggtitle("Total mutations") +
   xlab("Group") + ylab("Total Mutation frequency")  +
  scale_color_manual(values=col_vector) +
   geom_boxplot()
g2 <- ggplot(db_obs_freq, aes(x=GROUP, y=MU_FREQ_SEQ_S, fill=SAMPLE)) +
  theme_bw() + ggtitle("Total silent mutations") +
  xlab("Group") + ylab("Mutation frequency") +
  scale_color_manual(values=col_vector) +
  geom_boxplot()
g3 <- ggplot(db_obs_freq, aes(x=GROUP, y=MU_FREQ_SEQ_R, fill=SAMPLE)) +
  theme_bw() + ggtitle("Total replacement mutations") +
  xlab("Group") + ylab("Mutation frequency") +
  scale_color_manual(values=col_vector) +
  geom_boxplot()
CairoPNG(filename = "mutation_frequency_silent_replacement.png",
         width = 15000, height = 5000, pointsize = 12,
         fallback_resolution = 600,res = 600)
alakazam::gridPlot(g1, g2, g3, ncol=3)
dev.off()

# Calculate  mutation counts for individual CDRs and FWRs
db_obs_v <- observedMutations(all, sequenceColumn="SEQUENCE_IMGT",
                              germlineColumn="GERMLINE_IMGT_D_MASK",
                              regionDefinition=IMGT_V_BY_REGIONS,
                              frequency=FALSE, combine = T,
                              nproc=4)
head(db_obs_v_IgA)
## IgA Count
db_obs_v_IgA <- filter(db_obs_v, CREGION %in% c("Mouse-IGHA"))
db_obs_v_IgA$SAMPLE <- factor(db_obs_v_IgA$SAMPLE, levels = c("MM86_WT_NAIVE","MM87_KO_NAIVE","MM85-M3_WT","MM85-M4_WT","MM86-F1_WT","MM86-M6_KO","MM89-F1_KO","MM89-M3_KO"))
df <- db_obs_v_IgA %>% group_by(SAMPLE) %>% summarize(n = n())
g2 <- ggplot(db_obs_v_IgA, aes(x=SAMPLE, y=MU_COUNT, fill=SAMPLE)) +
  theme_bw() + ggtitle("Total  mutations") +
  xlab("Isotype") + ylab("Mutation frequency") +
  scale_fill_brewer(palette="Dark2") +
  geom_boxplot() + geom_text(data = df, aes(y = -2, label=n),
             position=position_dodge(width=1.0))
CairoPNG(filename = "IgA_CDR_mutation_count_replacement.png",
         width = 7500, height = 5000, pointsize = 12,
         fallback_resolution = 600,res = 600)
alakazam::gridPlot(g2, ncol=2)
dev.off()


# g3 <- ggplot(db_obs_v_IgA, aes(x=SAMPLE, y=MU_COUNT, fill=SAMPLE)) +
#   theme_bw() + ggtitle("Total  mutations") +
#   xlab("Isotype") + ylab("Mutation frequency") +
#   scale_fill_brewer(palette="Dark2") +geom_violin()+
#   stat_summary(fun.data=mean_sdl, mult=1, 
#                  geom="pointrange", color="black")
# 
# CairoPNG(filename = "IgA_CDR_mutation_count_replacement_violin.png",
#          width = 7500, height = 5000, pointsize = 12,
#          fallback_resolution = 600,res = 600)
# alakazam::gridPlot(g3, ncol=2)
# dev.off()




# Show new FWR mutation columns
db_obs_v %>% 
  select(SEQUENCE_ID, starts_with("MU_COUNT_FWR")) %>%
  head(n=4)


# Calculate aggregate CDR and FWR V-segment R and S mutation frequencies
db_obs_v <- observedMutations(db_obs_v, sequenceColumn="SEQUENCE_IMGT",
                              germlineColumn="GERMLINE_IMGT_D_MASK",
                              regionDefinition=IMGT_V,
                              frequency=TRUE, 
                              nproc=4)
# Show new CDR and FWR mutation frequency columns
db_obs_v %>% 
  select(SEQUENCE_ID, starts_with("MU_FREQ_")) %>%
  head(n=4)
g2 <- ggplot(db_obs_v, aes(x=SAMPLE, y=MU_FREQ_CDR_S, fill=CREGION)) +
  theme_bw() + ggtitle("CDR silent mutations") +
  xlab("Isotype") + ylab("Mutation frequency") +
  scale_fill_brewer(palette="Dark2") +
  geom_boxplot()
g3 <- ggplot(db_obs_v, aes(x=SAMPLE, y=MU_FREQ_CDR_R, fill=CREGION)) +
  theme_bw() + ggtitle("CDR replacement mutations") +
  xlab("Isotype") + ylab("Mutation frequency") +
  scale_fill_brewer(palette="Dark2") +
  geom_boxplot()
CairoPNG(filename = "CDR_mutation_frequency_silent_replacement.png",
         width = 10000, height = 5000, pointsize = 12,
         fallback_resolution = 600,res = 600)
alakazam::gridPlot(g2, g3, ncol=2)
dev.off()

## IgA 
db_obs_v_IgA <- filter(db_obs_v, CREGION %in% c("Mouse-IGHA"))
db_obs_v_IgA$SAMPLE <- factor(db_obs_v_IgA$SAMPLE, levels = c("MM86_WT_NAIVE","MM87_KO_NAIVE","MM85-M3_WT","MM85-M4_WT","MM86-F1_WT","MM86-M6_KO","MM89-F1_KO","MM89-M3_KO"))
df <- db_obs_v_IgA %>% group_by(SAMPLE) %>% summarize(n = n())

g2 <- ggplot(db_obs_v_IgA, aes(x=SAMPLE, y=MU_FREQ_CDR_S, fill=SAMPLE)) +
  theme_bw() + ggtitle("CDR silent mutations") +
  xlab("Isotype") + ylab("Mutation frequency") +
  scale_fill_brewer(palette="Dark2") +
  geom_boxplot() + geom_text(data = df, aes(y = -0.2, label=n),
                             position=position_dodge(width=1.0))
g3 <- ggplot(db_obs_v_IgA, aes(x=SAMPLE, y=MU_FREQ_CDR_R, fill=SAMPLE)) +
  theme_bw() + ggtitle("CDR replacement mutations") +
  xlab("Isotype") + ylab("Mutation frequency") +
  scale_fill_brewer(palette="Dark2") +
  geom_boxplot() + geom_text(data = df, aes(y = -0.2, label=n),
                             position=position_dodge(width=1.0))
CairoPNG(filename = "IgA_CDR_mutation_frequency_silent_replacement.png",
         width = 12500, height = 5000, pointsize = 12,
         fallback_resolution = 600,res = 600)
alakazam::gridPlot(g2, g3, ncol=2)
dev.off()




# Create substitution model using silent mutations
sub_matrix <- createSubstitutionMatrix(all, model="S")
# Create mutability model using silent mutations
mut_matrix <- createMutabilityMatrix(all, sub_matrix, model="S")

# Extend models to include ambiguous 5-mers
sub_matrix <- extendSubstitutionMatrix(sub_matrix)
mut_matrix <- extendMutabilityMatrix(mut_matrix)

# Create targeting model matrix from substitution and mutability matrices
tar_matrix <- createTargetingMatrix(sub_matrix, mut_matrix)

# Collapse sequences into clonal consensus
clone_db <- collapseClones(all, nproc=4)
# Create targeting model in one step using only silent mutations
# Use consensus sequence input and germline columns
model <- createTargetingModel(clone_db, model="S", sequenceColumn="CLONAL_SEQUENCE", 
                              germlineColumn="CLONAL_GERMLINE")

# Generate hedgehog plot of mutability model

CairoPNG(filename = "mutability_all_samples_A.png",
         width = 5000, height = 5000, pointsize = 12,
         fallback_resolution = 600,res = 600)
plotMutability(model, nucleotides="A", style="hedgehog") 
dev.off()

plotMutability(model, nucleotides="C", style="hedgehog")

# Generate bar plot of mutability model
plotMutability(model, nucleotides="G", style="bar") + facet_wrap(~SAMPLE)
plotMutability(model, nucleotides="T", style="bar")