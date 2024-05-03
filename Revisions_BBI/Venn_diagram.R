library(ggVennDiagram)
Polio_markers <- openxlsx::read.xlsx("/data1/meaneylab/eamon/Microglia_ePRS/Polioudakis cluster enriched genes.xlsx", 
                                     sheet = "Cluster enriched genes")

Polio_markers <- Polio_markers[Polio_markers$Cluster == "Mic",]

Tran_markers <- read.csv("/data1/meaneylab/eamon/Microglia_ePRS/sACC_markers_Tran.csv")


Genesets <- list(Adult_microglia_genes=Tran_markers$Gene,
                 Fetal_microglia_genes=Polio_markers$Gene)


plotme <- ggvenn(Genesets, show_percentage = F)+ 
  theme(text = element_text(size = 12)) +  
  scale_fill_brewer(palette = "Pastel1")  

ggsave("/data1/meaneylab/eamon/Microglia_ePRS/Revisions_BBI/Venn_diagram.pdf", 
       plot = plotme, width = 17, height = 12, units = "cm", bg="white")




