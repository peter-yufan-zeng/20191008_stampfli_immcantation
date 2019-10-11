library(circlize)
library(dplyr)
#Sys.setenv(R_GSCMD="C:/Program Files/gs/gs9.25/bin/gswin64c.exe")



  WTCircos = WT_Naive[, c("CLONE","CREGION")] 
  KOCircos = KO_Naive[, c("CLONE","CREGION")]
  head(WTCircos)
  
  colnames(WTCircos) <- c("from", "to")
  colnames(KOCircos) <- c("from", "to")
  
  WTCircos$value = 1
  KOCircos$value = 1
  
  WTCircos <- WTCircos[rowSums(is.na(WTCircos)) == 0,]
  KOCircos <- KOCircos[rowSums(is.na(KOCircos)) == 0,]
  
  
  circos.clear()
  circos.par(gap.after = c(rep(0.025, length(unique(WTCircos[[1]]))-1), 3, 
                           rep(3, length(unique(WTCircos[[2]]))-1), 3))
  pdf("./WT_Naive_Circos.pdf")
  par(cex = 0.25)
  chordDiagram(WTCircos)
  dev.off()
  
  circos.clear()
  circos.par(gap.after = c(rep(0.025, length(unique(KOCircos[[1]]))-1), 3, 
                           rep(3, length(unique(KOCircos[[2]]))-1), 3))
  pdf("./KO_Naive_Circos.pdf")
  par(cex = 0.25)
  chordDiagram(KOCircos)
  dev.off()


