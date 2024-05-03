so <- readRDS("/data1/meaneylab/eamon/Microglia_ePRS/Revisions_BBI/3730d395-82a4-4ed0-9d29-3307c20c44d7.rds")

Tran_markers <- read.csv("/data1/meaneylab/eamon/Microglia_ePRS/sACC_markers_Tran.csv")

# Human Ensembl database
ensembl_human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
Genes <- getBM(filters= "external_gene_name", attributes= c("ensembl_gene_id",
                                                         "external_gene_name"),
               values=Tran_markers$Gene,mart= ensembl_human)

micro_merged <- merge(Tran_markers,Genes, by.x = "Gene", by.y = "external_gene_name")

micro_ENSEMBL <- unique(micro_merged$ensembl_gene_id) %>%
  as.data.frame()

# Create score
so <- AddModuleScore(object = so, features = micro_ENSEMBL, name = "micro_enriched")
Idents(so) <- "author_cell_type"

# Remove cells with unclear ID
so <- subset(x = so, 
             idents = c('Microglia cluster 1', 
                        'Microglia cluster 2', 
                        'Microglia cluster 3',
                        'Microglia cluster 4', 
                        'Microglia cluster 5', 
                        'Microglia cluster 6',
                        'Microglia cluster 7', 
                        'Microglia cluster 8', 
                        'Microglia cluster 9'))
table(so$author_cell_type)
singcell_plot <- RidgePlot(so, features = 'micro_enriched1', group.by = 'author_cell_type') +
  ggtitle("Enrichment in Olah et al microglia clusters") +
  ylab('') + xlab('Module score') +
  theme(legend.position = 'none',
        plot.title = element_text(hjust = 0.5, face = 'plain'),
        panel.grid = element_blank(),
        axis.title.x = element_text(hjust = 0.5))

ggsave("/data1/meaneylab/eamon/Microglia_ePRS/Revisions_BBI/Single_Cell_plot.pdf", 
       plot = singcell_plot, width = 15, height = 7.5, units = "cm", bg="white")



#------------------------- GO -------------------------
PFC_genes <- read.csv("/data1/meaneylab/eamon/Microglia_ePRS/dlPFC_markers_Tran.csv")
background <- PFC_genes$gene
package.version("clusterProfiler")
Micro_GO <- enrichGO(gene         = unique(micro_merged$Gene),
                           OrgDb         = 'org.Hs.eg.db',
                           keyType       = 'SYMBOL',
                           ont           = "ALL",
                           pAdjustMethod = "BH",
                           pvalueCutoff  = 1,
                           qvalueCutoff  = 1,
                     universe = background)

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

GO_plot <- ggplot(data = Micro_GO_res, aes(x = Description, y = logPvalue)) +
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


ggsave("/data1/meaneylab/eamon/Microglia_ePRS/Revisions_BBI/GO_plot.pdf", 
       plot = GO_plot, width = 15, height = 7.5, units = "cm", bg="white")
