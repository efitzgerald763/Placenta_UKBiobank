Full_cyan <- openxlsx::read.xlsx ("/data1/meaneylab/eamon/Placenta study/DGIdb full cyan module.xlsx")
ePRS <- openxlsx::read.xlsx ("/data1/meaneylab/eamon/Placenta study/DGIdb full cyan ePRS.xlsx")
gene_interactions <- fread("/data1/meaneylab/eamon/Placenta_study/DGIdb 3+ genes full cyan.csv")



plotme <- ggplot(gene_interactions) +
  geom_col(aes(x = reorder(factor(Drug), Number), y = Number), fill = "darkolivegreen") +
  theme_minimal() +
  coord_flip() + 
  theme(panel.grid = element_blank(),
        text = element_text(size = 16)) +
  xlab("") +
  ylab("Number of interacting genes") +
  geom_text(data = gene_interactions, color = "white",
            aes(x = Drug, y = Number, label = Genes), 
            position = position_stack(0.5), size=3)
ggsave("/data1/meaneylab/eamon/Placenta_study/Revisions/Drug.pdf", 
       plot = plotme, width = 20, height = 15, units = "cm", bg="white")
