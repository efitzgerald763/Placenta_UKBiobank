so <- readRDS("/data1/meaneylab/eamon/Microglia_ePRS/Revisions_BBI/e9e38c8d-fd21-475d-acf6-f8994b003d70.rds")

Polio_markers <- openxlsx::read.xlsx("/data1/meaneylab/eamon/Microglia_ePRS/Polioudakis cluster enriched genes.xlsx", 
                                     sheet = "Cluster enriched genes")

micro <- Polio_markers[Polio_markers$Cluster == "Mic",]
rownames(so)

# Human Ensembl database
micro_ENSEMBL <- micro$Ensembl %>%
  as.data.frame()

# Create score
Idents(so) <- "assay"

# Remove cells with unclear ID
so <- subset(x = so, 
             idents = c("10x 3' v3"))


Idents(so) <- "cell_type"

# Remove cells with unclear ID
so <- subset(x = so, 
             idents = c("microglial cell"))

so <- AddModuleScore(object = so, features = micro_ENSEMBL, name = "micro_enriched")
Idents(so) <- "development_stage"

RidgePlot(so, features = 'micro_enriched1', group.by = 'region') +
  ggtitle("Enrichment in Olah et al microglia clusters") +
  ylab('') + xlab('Module score')
  theme(legend.position = 'none',
        plot.title = element_text(hjust = 0.5, face = 'plain'),
        panel.grid = element_blank())



#------------------------- GO -------------------------
Polio_back <- openxlsx::read.xlsx("/data1/meaneylab/eamon/Microglia_ePRS/Polioudakis cluster enriched genes.xlsx", 
                                      sheet = "Sheet1")
list_of_columns <- lapply(Polio_back, function(x) list(x))
backgound <- unique(unlist(list_of_columns))


Micro_GO <- enrichGO(gene         = micro$Gene,
                           OrgDb         = 'org.Hs.eg.db',
                           keyType       = 'SYMBOL',
                           ont           = "ALL",
                           pAdjustMethod = "BH",
                           pvalueCutoff  = 1,
                           qvalueCutoff  = 1)

Micro_GO_res <- Micro_GO@result

Micro_GO_res <- Micro_GO_res[with(Micro_GO_res,order(p.adjust)),]
Micro_GO_res <- Micro_GO_res[1:10,]



# Capatilize the first letter of each GO term
Capitalize <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}


Micro_GO_res$Description <- Capitalize(Micro_GO_res$Description)
Micro_GO_res$logPvalue <- -log10(Micro_GO_res$p.adjust)


# Use reorder to order by logPvalue
Micro_GO_res$Description <- reorder(Micro_GO_res$Description, Micro_GO_res$logPvalue)

GO <- ggplot(data = Micro_GO_res, aes(x = Description, y = logPvalue)) +
  geom_segment(aes(xend = Description, yend = 0), color = "black") +
  geom_point(size = 5) + 
  coord_flip() +
  theme_bw() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = 'none',
        text = element_text(size = 16)) +
  geom_hline(yintercept = 1.3, linetype = "dashed", color = "grey44") +
  xlab(NULL) +
  ylab(expression(paste("-log"[10]," (adj P-value)")))

ggsave("/data1/meaneylab/eamon/Microglia_ePRS/Revisions_BBI/GO_plot_Fetal.pdf", 
       plot = GO, width = 18, height = 7.5, units = "cm", bg="white")

