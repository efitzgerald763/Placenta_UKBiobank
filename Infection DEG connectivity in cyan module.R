library(awtools)
library(paletteer)


values <- fread ("/data1/meaneylab/eamon/Placenta_study/Cyan module DEG connectivity.csv")
values$Connectivity
ggplot(values, aes(x = Outcome, y = Connectivity, fill = Outcome)) +
  geom_col(width = 0.5) +
  scale_fill_manual(values=c("darkslategray3", "darkslategray")) +
  theme_minimal() +
  theme(axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, color = "grey20", angle = 45, hjust = 0.5, vjust = 0.5),
        legend.position = 'none',
        panel.grid = element_blank(),
        axis.line = element_line(color='grey90')) +
  xlab("") +
  ylab("Connectivity in cyan module") 

dev.off()

