#########################################################################################################################################################################
#########################################################################################################################################################################

##############################        ###############      FIGURE 3A (TCGA MEDIAN PLOT)				##############################################

library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)
library(digest)
library(cluster)
library(BBmisc)

options(scipen = 999)

abc=read.table('final_exp_tcga_mat.txt',"\t",header=T,check.names=F,row.names = 1)
a=abc
abc<- apply(a,2, function(x) x/var(x))


Label<-rep(c('BLCA','BRCA','CESC','COAD','GBM','HNSC','KIRC','KIRP','LIHC','LUAD','LUSC','PRAD','READ','STAD','THCA','UCEC'), c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1))

meta<-as.factor(Label)

myCol <- colorRampPalette(c('black', 'white', 'red'))(100)
myBreaks <- seq(min(abc), max(abc), length.out = 100)


ann <- data.frame(
  Sample = meta,
  stringsAsFactors = FALSE)


colours <- list(
  Sample = c('BLCA' = 'red','BRCA' = '#faf569','CESC' = '#eb4034','COAD' = '#eb9334','GBM' = '#ebd634','HNSC' = '#b7eb34','KIRC' = '#5feb34','KIRP' = '#f54295','LIHC' = '#34ebd6','LUAD' = '#34a8eb','LUSC'= '#3465eb','PRAD'='#8934eb','READ'='#c334eb','STAD'='#426640','THCA' = '#5febb4','UCEC' = '#18d3e3'))


Sample =list(
  title = '',
  title_position = 'topcenter',
  legend_direction = 'vertical',
  title_gp = gpar(fontsize = 12, fontface = 'bold'),
  labels_gp = gpar(fontsize = 12, fontface = 'bold'))




colAnn <- HeatmapAnnotation(
  df = ann,
  col = colours,
  annotation_height = 0.6,
  annotation_width = unit(1, 'cm'),
  gap = unit(1, 'mm'),
  annotation_legend_param =Sample)

kclus <- kmeans(abc, 12)
split <- paste0("Cluster\n", kclus$cluster)


abc<-as.matrix(abc)

hmap <- Heatmap(abc,
                split = split,
                cluster_row_slices = F,
                cluster_column_slices = F,
                
                name = 'Differential \n TAD-score',
                
                col = colorRamp2(myBreaks, myCol),
                
                heatmap_legend_param = list(
                  color_bar = 'continuous',
                  legend_direction = 'vertical',
                  legend_width = unit(8, 'cm'),
                  legend_height = unit(5.0, 'cm'),
                  title_position = 'topcenter',
                  title_gp=gpar(fontsize = 12, fontface = 'bold'),
                  labels_gp=gpar(fontsize = 12, fontface = 'bold')),
                
                row_split=split,
                show_row_dend =T,
                row_title_side = 'left',
                row_title_gp = gpar(fontsize = 12,  fontface = 'bold'),
                row_title_rot = 90,
                show_row_names = T,
                row_names_gp = gpar(fontsize = 0.5),
                row_names_side = 'right',
                row_dend_width = unit(15,'mm'),
                
                show_column_dend = T,
                column_km = 1,
                column_title = 'TADs expression across PAN-cancer',
                column_title_side = 'top',
                column_title_gp = gpar(fontsize = 15, fontface = 'bold'),
                column_title_rot = 0,
                show_column_names = T,
                column_names_gp = gpar(fontsize = 10, fontface = 'bold'),
                column_names_max_height = unit(10, 'cm'),
                column_dend_height = unit(8,'mm'),
                

                top_annotation = colAnn)

pdf(file="saving_plot_median.pdf",width=14,height = 9)

ht=draw(hmap,
        padding = unit(c(5, 10, 5, 10), "mm"),
        heatmap_legend_side = 'left',
        annotation_legend_side = 'right')
dev.off()



#########################################################################################################################################################################

##############          SURVIVAL SIGNIFICANT TAD FREQUENCY ACROSS PAN-CANCER              ##############  

#########################################################################################################################################################################


temp = list.files(pattern="*.csv")
myfiles = lapply(temp, read.csv)
tads = as.list(NA)
tad_p_value <-  as.list(NA)
tad_hazard_ratio <- as.list(NA)
naming = data.frame(matrix(NA,nrow = length(temp), ncol = 1))                                      
i = 1; b = 0; 
for ( i in (1:length(temp)))
{
  a <- as.data.frame(myfiles[i], sep = ',')
  a <- a[,-1]
  a <- a[order(a$p.value, decreasing = FALSE),]
  a <- dplyr::filter(a,a$p.value < 0.05) #decides cutt off
  rownames(a) = NULL
  tads[[i]] <- a$TAD
  tad_p_value[[i]] <- a$p.value
  tad_hazard_ratio[[i]] <- a$hazard.ratio
  if(i == 1)
  {naming[i,1] = dim(a)[1]; e = dim(a)[1]}
  if(i > 1)
  {naming[i,1] = dim(a)[1] + naming[(i-1),1]}
  b = b + dim(a)[1]
}
i = i + 1

df <- plyr::ldply(tads, cbind)
df_p_value <- plyr::ldply(tad_p_value, cbind)
df_hazard_value <- plyr::ldply(tad_hazard_ratio, cbind)


df <- cbind(df,df_p_value,df_hazard_value)
colnames(df) = c("TAD","p-value","hazard ratio")


substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

cancer_name = as.data.frame(NA)
df = cbind(df$TAD,cancer_name,df$`p-value`,df$`hazard ratio`)
colnames(df) = c("TAD","Cancer in which p <0.05 ","p-value","hazard ratio")


i=1 ;
for( i in (1:length(temp)))
{
  t = substrRight(temp[i],9)
  t= substr(t,1,4)
  if(i == 1)
  {df[1:naming[i,1],2] = t}
  if(i > 1)
  {df[c(naming[i-1,1]):c(naming[i,1]),2] = t}  
}
c <- as.data.frame(table(df$TAD))
c = c[order(c$Freq, decreasing = TRUE),]
rownames(c) = NULL




i =1 ;for(i in (1:dim(c)[1]))
{
  t = c[i,1]
  m  = dplyr::filter(df,df$TAD == t)
  freq = dim(m)[1]
  m = cbind(m$TAD,freq,m[,2:4])
  if( i == 1) { m1 = m}
  if( i >1) { m1 = rbind(m1,m)}
  }

colnames(m1)[2] = c("frequency")
write.csv(c,"p less than 0.05 TAD-overlapping frequency in 20 cancers from survival 51-49 rule.csv", na = "")
write.csv(m1,'p less than 0.05 TAD-overlapping frequency + p value + hazard ratio in 20 cancers from survival 51-49 rule.csv', na = "")


#########################################################################################################################################################################
#########################################################################################################################################################################
#########################################################################################################################################################################

####################        FIGURE 3D (Survival across multiple cancers)    #####################
library(data.table)
library(reshape2)
library(ggplot2)
library(ggthemes)

abc=read.table('p less than 0.05 TAD-overlapping frequency + p value + hazard ratio in 20 cancers from survival 51-49 rule.csv',sep="\t",header=T,check.names=F)
data=subset(abc,abc$frequency>7)

data2=abc[pos_1,]
data3=rbind(data,data2)
data3=data3[,c(1,3,4)]
aqm2=acast(data3, data3$TAD ~ data3$Cancer)
aqm2[is.na(aqm2)]=1
aqm2=melt(aqm2)
colnames(aqm2)=c('TAD','CANCER','p_value')
#hazard
data3=rbind(data,data2)
data4=data3[,c(1,3,5)]
aqm3=acast(data4, data4$TAD ~ data4$Cancer)
aqm3[is.na(aqm3)]=0
aqm3=melt(aqm3)
colnames(aqm3)=c('TAD','CANCER','hazard_ratio')


dot_data=merge(aqm2,aqm3,by=c('TAD','CANCER'),all=T)
dot_data2=dot_data
dot_data2$p_value=(-log10(dot_data2$p_value))
dot_data2$hazard_ratio=log10(dot_data2$hazard_ratio)
dot_data2$hazard_ratio[is.infinite(dot_data2$hazard_ratio)]=0


p <- ggplot(data = dot_data2, 
            aes(x=dot_data2$CANCER, y=dot_data2$TAD, color=dot_data2$hazard_ratio,size=dot_data2$p_value))+
  scale_colour_gradientn(colours=c("red","#f26c18","#38588a","#02000a"),  breaks = c(-0.8,-0.4,0,0.4)) +
  ggtitle("TAD survival across multiple cancers")+
  geom_point()+
  scale_size(range = c(5, 15))+
  theme(plot.background=element_rect(fill="white"),
        plot.margin = unit(c(5, 10, 10, 10), "cm")) 

p2=p +
  theme_bw()+
  theme(plot.title = element_text(size=30, face="bold", margin = margin(0, 0, 10, 0)))+
  theme(
    axis.title.x = element_text(color="black", vjust=-0.55,size = 22,face="bold"),
    axis.title.y = element_text(color="black" , vjust=0.55,size = 22,face="bold")   
  )+
  theme(axis.text.x=element_text(color="black",size=16, angle = 30, hjust=1),
        axis.text.y=element_text(color="black",size=16, vjust=0.5))+
  labs(x='CANCER',y='TADs')+
  theme(legend.title = element_text(size=20, face="bold"))+
  theme(legend.position="right")+ theme(legend.text=element_text(size=18))

ggsave('dotplot_survival.pdf',width = 18, height = 10)

#########################################################################################################################################################################
#########################################################################################################################################################################
#########################################################################################################################################################################
#########################################################################################################################################################################
#########################################################################################################################################################################
