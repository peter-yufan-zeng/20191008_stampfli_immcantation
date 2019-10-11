
db_props <- aminoAcidProperties(all, seq="JUNCTION", nt=TRUE, trim=TRUE, 
                                label="CDR3")
#db_props <- filter(db_props, CREGION %in% c("Mouse-IGHA"))
# The full set of properties are calculated by default
dplyr::select(db_props[1:3, ], starts_with("CDR3"))

tmp_theme <- theme_bw() + theme(legend.position="bottom")
g1
# Generate plots for a four of the properties
g1 <- ggplot(db_props, aes(x=GROUP, y=CDR3_AA_LENGTH, fill=SAMPLE)) + tmp_theme +
  ggtitle("CDR3 length") + 
  xlab("Group") + ylab("Amino acids") +
  scale_color_manual(values=col_vector) +
  geom_boxplot(aes(fill=SAMPLE))
g2 <- ggplot(db_props, aes(x=GROUP, y=CDR3_AA_GRAVY, fill=SAMPLE)) + tmp_theme + 
  ggtitle("CDR3 hydrophobicity") + 
  xlab("Group") + ylab("GRAVY") +
  scale_color_manual(values=col_vector) +
  geom_boxplot(aes(fill=SAMPLE))
g3 <- ggplot(db_props, aes(x=GROUP, y=CDR3_AA_BASIC, fill=SAMPLE)) + tmp_theme +
  ggtitle("CDR3 basic residues") + 
  xlab("Group") + ylab("Basic residues") +
  scale_y_continuous(labels=scales::percent) +
  scale_color_manual(values=col_vector) +
  geom_boxplot(aes(fill=SAMPLE))
g4 <- ggplot(db_props, aes(x=GROUP, y=CDR3_AA_ACIDIC, fill=SAMPLE)) + tmp_theme +
  ggtitle("CDR3 acidic residues") + 
  xlab("Group") + ylab("Acidic residues") +
  scale_y_continuous(labels=scales::percent) +
  scale_color_manual(values=col_vector) +
  geom_boxplot(aes(fill=SAMPLE))

# Plot in a 2x2 grid
CairoPNG(filename = "AA_Properties.png",
         width = 10000, height = 5000, pointsize = 12,
         fallback_resolution = 600,res = 600)
gridPlot(g1, g2, g3, g4, ncol=2)
dev.off()







db_props <- aminoAcidProperties(WT, seq="JUNCTION", nt=TRUE, trim=TRUE, 
                                label="CDR3")
#db_props <- filter(db_props, CREGION %in% c("Mouse-IGHA"))
# The full set of properties are calculated by default
dplyr::select(db_props[1:3, ], starts_with("CDR3"))

tmp_theme <- theme_bw() + theme(legend.position="bottom")

# Generate plots for a four of the properties
g1 <- ggplot(db_props, aes(x=CREGION, y=CDR3_AA_LENGTH)) + tmp_theme +
  ggtitle("CDR3 length") + 
  xlab("Isotype") + ylab("Amino acids") +
  scale_fill_brewer(name = "Isotype", palette="Dark2") +
  geom_boxplot(aes(fill=CREGION))
g2 <- ggplot(db_props, aes(x=CREGION, y=CDR3_AA_GRAVY)) + tmp_theme + 
  ggtitle("CDR3 hydrophobicity") + 
  xlab("Isotype") + ylab("GRAVY") +
  scale_fill_brewer(name = "Isotype", palette="Dark2") +
  geom_boxplot(aes(fill=CREGION))
g3 <- ggplot(db_props, aes(x=CREGION, y=CDR3_AA_BASIC)) + tmp_theme +
  ggtitle("CDR3 basic residues") + 
  xlab("Isotype") + ylab("Basic residues") +
  scale_y_continuous(labels=scales::percent) +
  scale_fill_brewer(name = "Isotype", palette="Dark2") +
  geom_boxplot(aes(fill=CREGION))
g4 <- ggplot(db_props, aes(x=CREGION, y=CDR3_AA_ACIDIC)) + tmp_theme +
  ggtitle("CDR3 acidic residues") + 
  xlab("Isotype") + ylab("Acidic residues") +
  scale_y_continuous(labels=scales::percent) +
  scale_fill_brewer(name = "Isotype", palette="Dark2") +
  geom_boxplot(aes(fill=CREGION))

# Plot in a 2x2 grid
CairoPNG(filename = "WT_IgA_only_AA_Properties.png",
         width = 10000, height = 5000, pointsize = 12,
         fallback_resolution = 600,res = 600)
gridPlot(g1, g2, g3, g4, ncol=2)
dev.off()


db_props <- aminoAcidProperties(all, seq="JUNCTION", nt=TRUE, trim=TRUE, 
                                label="CDR3")
db_props_IgA <- filter(db_props, CREGION %in% c("Mouse-IGHA"))
tmp_theme <- theme_bw() + theme(legend.position="bottom")

# Generate plots for a four of the properties
g1 <- ggplot(db_props_IgA, aes(x=CREGION, y=CDR3_AA_LENGTH)) + tmp_theme +
  ggtitle("CDR3 length") + 
  xlab("Isotype") + ylab("Amino acids") +
  scale_fill_brewer(name = "Isotype", palette="Dark2") +
  geom_boxplot(aes(fill=SAMPLE))
g2 <- ggplot(db_props_IgA, aes(x=CREGION, y=CDR3_AA_GRAVY)) + tmp_theme + 
  ggtitle("CDR3 hydrophobicity") + 
  xlab("Isotype") + ylab("GRAVY") +
  scale_fill_brewer(name = "Isotype", palette="Dark2") +
  geom_boxplot(aes(fill=SAMPLE))
g3 <- ggplot(db_props_IgA, aes(x=CREGION, y=CDR3_AA_BASIC)) + tmp_theme +
  ggtitle("CDR3 basic residues") + 
  xlab("Isotype") + ylab("Basic residues") +
  scale_y_continuous(labels=scales::percent) +
  scale_fill_brewer(name = "Isotype", palette="Dark2") +
  geom_boxplot(aes(fill=SAMPLE))
g4 <- ggplot(db_props_IgA, aes(x=CREGION, y=CDR3_AA_ACIDIC)) + tmp_theme +
  ggtitle("CDR3 acidic residues") + 
  xlab("Isotype") + ylab("Acidic residues") +
  scale_y_continuous(labels=scales::percent) +
  scale_fill_brewer(name = "Isotype", palette="Dark2") +
  geom_boxplot(aes(fill=SAMPLE))

# Plot in a 2x2 grid
CairoPNG(filename = "IgA_only_AA_Properties.png",
         width = 10000, height = 5000, pointsize = 12,
         fallback_resolution = 600,res = 300)
gridPlot(g1, g2, g3, g4, ncol=2)
dev.off()



