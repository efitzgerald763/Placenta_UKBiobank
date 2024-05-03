library(tidyverse)

Total <- read_csv ("/data1/meaneylab/eamon/Placenta_study/NegControl_ssGSEA_GUSTO_reg_PCs_sex_26Jan2022.csv")

# To order as in file
Total$Predictor <- as.character(Total$Predictor)
Total$Predictor  <- factor(Total$Predictor, levels=unique(Total$Predictor))

plotD <- ggplot(Total, aes(x=Predictor, y=Beta, ymin=Beta-1.96*SE_Beta, ymax=Beta+1.96*SE_Beta, colour=PGS)) +
  geom_pointrange() +
  scale_color_manual(values=c("#E69F00","grey69")) +
  geom_hline(yintercept=0, lty=2, color = "grey60") +  
  coord_flip() +  # flip coordinates (puts labels on y axis)
  ylab("Coefficient for PGS (95% CI)")  +
  theme_minimal() + # use a white background
  labs(x = "") + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title=element_blank(),
        legend.position = 'bottom',
        text = element_text(size = 16),
        plot.margin = margin(1,1,1,0, "cm")) 
plotD + guides(color=guide_legend(nrow=2))




ggsave(path = path, width = width, height = height, device='tiff', dpi=700)
