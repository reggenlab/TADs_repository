
#########################################################################################################################################################################
#########################################################################################################################################################################
#########################################################################################################################################################################
###############      FIGURE 2A 2B				##############################################

ccle_matrix=read.table('TADscore-GSVA-unionTADs-CCLE-CTRP_without_white_space_4dec22_random.csv',header=T,sep='\t',row.names=1)
ccle_colname_tcgabased=read.table("final_colname_ccle_to_TCGA.txt")
a=ccle_matrix
colnames(a)=ccle_colname_tcgabased[,2]



######
######   Avoid for figure 2B   ############
abc_blca=a[,grepl("HNSC",colnames(a),ignore.case = T)]
a=abc_blca
######   Avoid for figure 2B   ############
###########################################


b=read.table('ctrp2-ccle-drugMat.txt',sep="\t")
b1=-log10(b*10^-6)
b1[is.na(b1)]=0
b1[is.infinite(as.matrix(b1))]=0
colnames(b1)=ccle_colname_tcgabased[,2]


######   Avoid for figure 2B   ############
abc_blca=b1[,grepl("HNSC",colnames(b1),ignore.case = T)]
b1=abc_blca
######   Avoid for figure 2B   ############
###########################################

cordd=cor(t(a),t(b1),method='spearman')

abc<-as.matrix(cordd)

############### cutoff based on correlation ###################
library(reshape2)
data2=melt(as.matrix(abc))
data3=subset(data2,abs(data2$value)>0.3)
data4=acast(data3,data3$Var1~data3$Var2)
dim(data4)

data4[is.na(data4)]=0
sum(is.na(data4))

abc=as.matrix(data4)


#################     correlation heatmap

myCol <- colorRampPalette(c('Red', '#ced7d8','Blue'))(1000)
myBreaks <- seq(min(abc), max(abc), length.out = 1000)

abc<-as.matrix(abc)

library(dendextend)
set.seed(100)
col_dend = as.dendrogram(hclust(as.dist(1 - cor(abc,method = "spearman"))))
col_dend = color_branches(col_dend, k = 8)


set.seed(800)
kclus <- kmeans(abc, 8)



#set.seed(800)
hmap <- Heatmap(abc,

                split = paste0("Cluster\n", kclus$cluster),
                cluster_row_slices = F,
                cluster_column_slices = F,
                
                name = 'Correlation',
                
                col = colorRamp2(myBreaks, myCol),
                
                heatmap_legend_param = list(
                  color_bar = 'continuous',
                  legend_direction = 'vertical',
                  legend_width = unit(8, 'cm'),
                  legend_height = unit(5.0, 'cm'),
                  title_position = 'topcenter',
                  title_gp=gpar(fontsize = 12, fontface = 'bold'),
                  labels_gp=gpar(fontsize = 12, fontface = 'bold')),
                

                cluster_rows =TRUE,
                show_row_dend =T,

                row_title_side = 'left',
                row_title_gp = gpar(fontsize = 12,  fontface = 'bold'),
                row_title_rot = 90,
                show_row_names = T,
                row_names_gp = gpar(fontsize = 0.5),

                row_names_side = 'right',
                row_dend_width = unit(15,'mm'),
                
                cluster_columns = col_dend,
                show_column_dend = T,

                column_title = 'TADs correlation in BRCA cancer with drugs',
                column_title_side = 'top',
                column_title_gp = gpar(fontsize = 15, fontface = 'bold'),
                column_title_rot = 0,
                show_column_names = T,
                column_names_gp = gpar(fontsize = 1.5, fontface = 'bold'),
                column_names_max_height = unit(10, 'cm'),
                column_dend_height = unit(15,'mm'),

                row_dend_reorder =TRUE,
                column_dend_reorder=FALSE,
                clustering_distance_rows = function(x) as.dist(1 - cor(t(x),method = "spearman")),)

pdf(file="heatmap.pdf",width=15)

ht=draw(hmap,
        padding = unit(c(5, 10, 5, 10), "mm"),
        heatmap_legend_side = 'left')

dev.off()


#########################################################################################################################################################################
#########################################################################################################################################################################


