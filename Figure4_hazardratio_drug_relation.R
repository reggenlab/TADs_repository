#########################################################################################################################################################################

######## ####################        FIGURE 4A			#####################

library(tidyverse)
library(hrbrthemes)
library(viridis)

mat_new=read.table('boxplot_hazard_all_table.txt',header = T,sep='\t')
set.seed(18)
cols = rainbow(6, s=.5, v=.6)[sample(1:6,6)]
mat_new$category<- factor(mat_new$category, levels = unique(mat_new$category))


# Plot
p=mat_new %>%
  ggplot( aes(x=category, y=hazard_ratio,fill=category)) +
  geom_boxplot(notch=TRUE)+
  theme_ipsum() +
  theme(
    legend.position="right",
    plot.title = element_text(size=25)
  ) +
  ggtitle("hazard ratio for different TADs ")

p2=p +theme_classic()+ 
  theme(plot.title = element_text(size=25, face="bold", margin = margin(10, 10, 10, 10)))+
  theme(
    axis.title.x = element_text(color="black", vjust=-0.55,size = 20),
    axis.title.y = element_text(color="black" , vjust=0.55,size = 20)   
  )+
  theme(axis.text.x=element_text(color="black",size=18, vjust=0.5),
        axis.text.y=element_text(color="black",size=18, vjust=0.5))+
  labs(x='category',y='hazard ratio')+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  theme(legend.position="none",legend.text=element_text(size=18))+
  scale_fill_manual(values=cols)

ggsave("all_topbox_hazard.pdf",height = 10,width = 18)

#########################################################################################################################################################################
