library(openxlsx)
term_placenta <- readRDS("/data1/meaneylab/eamon/Placenta_study/Revisons_2.0/GSE182381_processed_single_cell_Seurat_object.rda")
table(Modules$Module.Name)

DefaultAssay(term_placenta) <- 'RNA'

# Module Score
Modules <- read.xlsx("/data1/meaneylab/eamon/Placenta_study/Modules.xlsx")
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))

Cyan <- Modules[Modules$Module.Name == 'cyan',]
Cyan <- getBM(filters= "ensembl_gene_id", 
              attributes= c("ensembl_gene_id","external_gene_name"),
              values=Cyan$Ensembl.Gene.ID,mart= mart) %>%
  dplyr::select(external_gene_name)
term_placenta <- AddModuleScore(object = term_placenta, features = Cyan, name = "Cyan")

darkorange <- Modules[Modules$Module.Name == 'darkorange',]
darkorange <- getBM(filters= "ensembl_gene_id", 
                    attributes= c("ensembl_gene_id","external_gene_name"),
                    values=darkorange$Ensembl.Gene.ID,mart= mart) %>%
  dplyr::select(external_gene_name)
term_placenta <- AddModuleScore(object = term_placenta, features = darkorange, name = "darkorange")

darkturquoise <- Modules[Modules$Module.Name == 'darkturquoise',]
darkturquoise <- getBM(filters= "ensembl_gene_id", 
                       attributes= c("ensembl_gene_id","external_gene_name"),
                       values=darkturquoise$Ensembl.Gene.ID,mart= mart) %>%
  dplyr::select(external_gene_name)
term_placenta <- AddModuleScore(object = term_placenta, features = darkturquoise, name = "darkturquoise")

grey60 <- Modules[Modules$Module.Name == 'grey60',]
grey60 <- getBM(filters= "ensembl_gene_id", 
                attributes= c("ensembl_gene_id","external_gene_name"),
                values=grey60$Ensembl.Gene.ID,mart= mart) %>%
  dplyr::select(external_gene_name)
term_placenta <- AddModuleScore(object = term_placenta, features = grey60, name = "grey60")

royalblue <- Modules[Modules$Module.Name == 'royalblue',]
royalblue <- getBM(filters= "ensembl_gene_id", 
                   attributes= c("ensembl_gene_id","external_gene_name"),
                   values=royalblue$Ensembl.Gene.ID,mart= mart) %>%
  dplyr::select(external_gene_name)
term_placenta <- AddModuleScore(object = term_placenta, features = royalblue, name = "royalblue")

table(term_placenta$cell.type)

Meta <- term_placenta@meta.data

df <- aggregate( Cyan1 ~ cell.type, Meta, mean )

df1 <- aggregate( darkorange1 ~ cell.type, Meta, mean )
df<- full_join(df,df1, by = 'cell.type')

df1 <- aggregate( darkturquoise1 ~ cell.type, Meta, mean )
df<- full_join(df,df1, by = 'cell.type')

df1 <- aggregate( grey601 ~ cell.type, Meta, mean )
df<- full_join(df,df1, by = 'cell.type')

df1 <- aggregate( royalblue1 ~ cell.type, Meta, mean )
df<- full_join(df,df1, by = 'cell.type')

df <- df%>%
  dplyr::rename(Cyan = Cyan1)%>%
  dplyr::rename(Darkorange = darkorange1)%>%
  dplyr::rename(Darkturquoise = darkturquoise1)%>%
  dplyr::rename(Grey60 = grey601)%>%
  dplyr::rename(Royalblue = royalblue1)


df_num = as.matrix(df[,2:6])
rownames(df_num) <- df$cell.type
df_num_scale = scale(df_num)

pheatmap(
  mat               = df_num_scale,
  color             = cividis(50),
  border_color      = NA,
  cluster_cols      = T,
  cluster_rows      = FALSE,
  show_colnames     = TRUE,
  show_rownames     = TRUE,
  drop_levels       = F,
  filename = '/data1/meaneylab/eamon/Placenta_study/Revisons_2.0/single_cell_heatmap.pdf',
  width = 7, height = 7,
  fontsize          = 15,
  main              = "",
  angle_col         = 90)
