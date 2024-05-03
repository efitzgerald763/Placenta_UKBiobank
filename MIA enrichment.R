

Modules <- read.xlsx("/data1/meaneylab/eamon/Placenta_study/Modules.xlsx")
Cyan <- Modules[Modules$Module.Name == 'cyan',]
Cyan <- Cyan$Ensembl.Gene.ID

Leukens <- read.csv("/data1/meaneylab/eamon/Placenta_study/Revisions/Leukens 2022 DEGs.csv")
Leukens <- Leukens$ortholog_ensg %>%
  as.data.frame()
Leukens$. <- gsub("N/A","",Leukens$.)
Leukens <- Leukens[!(is.na(Leukens$.) | Leukens$.==""), ]


Conor <- read.csv("/data1/meaneylab/eamon/Placenta_study/Revisions/Connor et al 2022 DEGs.csv")
Conor <- Conor$ortholog_ensg %>%
  as.data.frame()
Conor$. <- gsub("N/A","",Conor$.)
Conor <- Conor[!(is.na(Conor$.) | Conor$.==""), ]

go.obj <- newGeneOverlap(Cyan, Leukens, genome.size =31097)

go.obj <- testGeneOverlap(go.obj)

go.obj
