load("/data1/meaneylab/eamon/Placenta_study/Single_cell_for_revision/Term_placenta.RData")
#-------------------------------------

meta_full <- term_placenta@meta.data
table(meta_full$annotation_merged)
meta_full$celltype <-ifelse(grepl("APC", meta_full$annotation_merged, ignore.case=TRUE), "Antigen presenting cell",
                            ifelse(grepl("Bcell", meta_full$annotation_merged, ignore.case=TRUE), "B-cell",
                                   ifelse(grepl("DSC", meta_full$annotation_merged, ignore.case=TRUE), "Stromal",
                                          ifelse(grepl("vil.FB", meta_full$annotation_merged, ignore.case=TRUE), "Fibroblast",
                                                 ifelse(grepl("dec.FB", meta_full$annotation_merged, ignore.case=TRUE), "Fibroblast",
                                                        ifelse(grepl("SMC", meta_full$annotation_merged, ignore.case=TRUE), "Smooth muscle cell",
                                                               ifelse(grepl("Gran", meta_full$annotation_merged, ignore.case=TRUE), "Granulocyte",
                                                                      ifelse(grepl("Mono", meta_full$annotation_merged, ignore.case=TRUE), "Monocyte",
                                                                             ifelse(grepl("NK_", meta_full$annotation_merged, ignore.case=TRUE), "Natural killer cell",
                                                                                    ifelse(grepl("Tcell", meta_full$annotation_merged, ignore.case=TRUE), "T-cell",
                                                                                           ifelse(grepl("Ery", meta_full$annotation_merged, ignore.case=TRUE), "Erythrocyte",
                                                                                                  ifelse(grepl("CT", meta_full$annotation_merged, ignore.case=TRUE), "Trophoblast",
                                                                                                         ifelse(grepl("VT", meta_full$annotation_merged, ignore.case=TRUE), "Trophoblast",
                                                                                                                ifelse(grepl("vil.Hofb", meta_full$annotation_merged, ignore.case=TRUE), "Hofbauer cell", NA))))))))))))))

table(meta_full$celltype)
term_placenta@meta.data <- meta_full
DefaultAssay(term_placenta) <- 'RNA'

# Module Score
Modules <- read.xlsx("/data1/meaneylab/eamon/Placenta_study/Modules.xlsx")
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))

Genes <- getBM(filters= "ensembl_gene_id", 
              attributes= c("ensembl_gene_id","external_gene_name"),
              values=Modules$Ensembl.Gene.ID,mart= mart) %>%
  dplyr::rename(Ensembl.Gene.ID=ensembl_gene_id)

Modules <- merge(Genes,Modules, by = "Ensembl.Gene.ID")

# Get unique module names
unique_module_names <- unique(Modules$Module.Name)

# Iterate over each unique module name
for (module_name in unique_module_names) {
  module_data <- Modules[Modules$Module.Name == module_name,] %>%
    dplyr::select(external_gene_name) %>%
    as.data.frame()
  term_placenta <- AddModuleScore(object = term_placenta, features = module_data, name = module_name)
}
str(df_num)
Meta <- term_placenta@meta.data
df_num <- as.data.frame(Meta[,16:44])
summary_df <- df_num %>%
  group_by(celltype) %>%
  summarise(across(.cols = turquoise1:darkorange1, .fns = mean)) %>%
  na.omit()

colnames(summary_df) <- gsub("1","", colnames(summary_df))
summary_df <- column_to_rownames(summary_df, var = "celltype")
df_num_scale = scale(summary_df)

pheatmap(
  mat               = df_num_scale,
  color             = cividis(50),
  border_color      = NA,
  cluster_cols      = T,
  cluster_rows      = FALSE,
  show_colnames     = TRUE,
  show_rownames     = TRUE,
  drop_levels       = F,
  filename = '/data1/meaneylab/eamon/Placenta_study/Revisons_2.0/All_modules_single_cell_heatmap.pdf',
  width = 15, height = 8,
  fontsize          = 15,
  main              = "",
  angle_col         = 90)

