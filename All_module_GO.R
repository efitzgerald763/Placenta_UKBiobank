library(gProfileR)
library(gprofiler2)

Modules <- read.xlsx("/data1/meaneylab/eamon/Placenta_study/Modules.xlsx")
Mods <- unstack(Modules)
#background <- Modules$Ensembl.Gene.ID

gostres <- gost(query =  Mods, 
                organism = "hsapiens", ordered_query = FALSE, significant = TRUE, exclude_iea = TRUE, 
                measure_underrepresentation = FALSE, evcodes = FALSE, 
                user_threshold = 0.05, correction_method = "g_SCS", 
                domain_scope = "annotated", custom_bg = NULL, 
                numeric_ns = "", sources = NULL, as_short_link = FALSE)
go <- gostres$result

go <- go[go$source == 'GO:BP',]
go <- go %>%
  group_by(query) %>%
  slice_min(order_by = p_value, n = 10)

go$logPvalue <- -log(go$p_value)

Capitalize <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

go$term_name <- Capitalize(go$term_name)
go$query <- Capitalize(go$query)

cyan <- go[go$query == 'cyan',]
darkorange <- go[go$query == 'darkorange',]
darkturquoise <- go[go$query == 'darkturquoise',]
grey60 <- go[go$query == 'grey60',]
royalblue <- go[go$query == 'royalblue',]

total <- rbind(cyan,darkorange, darkturquoise, grey60, royalblue)



plotme <- ggplot(total, aes(x=term_name, y=query, color = logPvalue, size = logPvalue )) + 
  geom_point() +
  scale_color_paletteer_c(`"ggthemes::Classic Area-Brown"`,name = "-log P-value") +
  theme_minimal() + labs(size="")+
  coord_flip()+
  theme(panel.grid.major = element_blank(),
        legend.position = 'bottom',
        text = element_text(size=16),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+
  scale_size(guide = 'none') +
  xlab('')+ ylab('') 
ggsave("/data1/meaneylab/eamon/Placenta_study/Revisions/Inflam_GO_dotplot.pdf", 
       plot = plotme, width = 20, height = 30, units = "cm", bg="white")

go <- apply(go,2,as.character)
write.csv(go,file = "/data1/meaneylab/eamon/Placenta_study/Revisions/t10_GO_all_mods.csv")

