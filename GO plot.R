library("viridis")    
library(tidyverse)

GO <-openxlsx::read.xlsx("/data1/meaneylab/eamon/Placenta_study/GO plot data.xlsx")



#To order axis correctly
#Turn your 'Cluster' column into a character vector
GO$term_name <- as.character(GO$term_name)

#Then turn it back into a factor with the levels in the correct order
GO$term_name <- factor(GO$term_name, levels=unique(GO$term_name))

plotB <- ggplot(GO, aes(x = term_name, y = `negative_log10_of_adjusted_p_value`)) +
  geom_col(width = 0.6, fill = "cyan3") +
  theme_minimal() +
  theme(axis.text.x = element_text(color = "grey20", angle = 45, hjust = 1, vjust = 1),
        panel.grid = element_blank(),
        text = element_text(size=14),
        axis.line = element_line(color='grey90'),
        legend.position="none") +
  geom_hline(yintercept=1.3, linetype="dashed", color = "grey44") +
  xlab("") +
  ylab("-log Adj p-value") +
  coord_flip()

