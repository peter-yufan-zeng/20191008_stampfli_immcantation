library(shazam)
library(Cairo)

outputMutabilityFigures <- function(db, name){
clone_db <- collapseClones(db, nproc=4)
# Create targeting model in one step using only silent mutations
# Use consensus sequence input and germline columns
model <- createTargetingModel(clone_db, model="S", sequenceColumn="CLONAL_SEQUENCE", 
                              germlineColumn="CLONAL_GERMLINE")

# Generate hedgehog plot of mutability model
setwd("c:/Users/zyfni/Desktop/20190527_immcantation_400_200/mouse/20190527_jon_mouse_400_200/")
CairoPNG(filename = paste("./mutability_figures/",name, "_A_gseaplot.png",sep = ""),
         width = 5000, height = 5000, pointsize = 12,
         fallback_resolution = 600,res = 600)
plotMutability(model, nucleotides="A", style="hedgehog") 
dev.off()

CairoPNG(filename = paste("./mutability_figures/",name, "_C_gseaplot.png",sep = ""),
         width = 5000, height = 5000, pointsize = 12,
         fallback_resolution = 600,res = 600)
plotMutability(model, nucleotides="C", style="hedgehog") 
dev.off()

CairoPNG(filename = paste("./mutability_figures/",name, "_T_gseaplot.png",sep = ""),
         width = 5000, height = 5000, pointsize = 12,
         fallback_resolution = 600,res = 600)
plotMutability(model, nucleotides="T", style="hedgehog") 
dev.off()

CairoPNG(filename = paste("./mutability_figures/",name, "_G_gseaplot.png",sep = ""),
         width = 5000, height = 5000, pointsize = 12,
         fallback_resolution = 600,res = 600)
plotMutability(model, nucleotides="G", style="hedgehog") 
dev.off()
}

outputMutabilityFigures(WT1, "MM85-M3_WT")
outputMutabilityFigures(WT2, "MM85-M4_WT")
outputMutabilityFigures(KO1, "MM89-F1_KO")
outputMutabilityFigures(KO2, "MM89-M3_KO")
outputMutabilityFigures(WT_Naive, "MM86_WT_NAIVE")
