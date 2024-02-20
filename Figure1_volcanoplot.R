### HN120MCR Vs HN120M ###
geneset = read.csv('Non_redundant_UnionTADs_for_data_genes.csv',header=F,strip.white=T,stringsAsFactor=F)
mat=read.table('unipath_pathwayscore_with_UnionTADs_5dec.txt',header=T,row.names = 1)
col=read.table('GSE117872_good_Data_cellinfo1.txt')
col1=col[,1]
wilk= matrix(1 , nrow(mat) , 1) ;
FC= matrix(1 , nrow(mat) , 1) ;
pos = which( col1 == 'HN120MCR') ;
pos1 = which(col1 == 'HN120M'); 
for(k in 1:nrow(mat)){
  tri = wilcox.test( as.numeric(mat[k,pos]) , as.numeric(mat[k,pos1]) )
  wilk[k] = tri$p.value 
  FC[k] = median(as.numeric(mat[k,pos])) - median(as.numeric(mat[k,pos1])) ;
}
HN120MCR_HN120M_TAD=geneset[,1:1]
HN120MCR_HN120M_TAD=data.frame(HN120MCR_HN120M_TAD)
HN120MCR_HN120M_TAD1=paste0("TAD_",1:nrow(HN120MCR_HN120M_TAD))
HN120MCR_HN120M_geneset=geneset[,2:ncol(geneset)]
HN120MCR_HN120M=cbind(HN120MCR_HN120M_TAD1,HN120MCR_HN120M_TAD,wilk,FC,HN120MCR_HN120M_geneset)
colnames(HN120MCR_HN120M)[1]="TADs_number"
colnames(HN120MCR_HN120M)[2]="HN120MCR_HN120M_TADs"
write.table(HN120MCR_HN120M, file = "HN120MCR_HN120M_union.csv", sep=",", row.names = F , col.names = T, quote=F) ;

###### volcano #######
library(EnhancedVolcano)
HN120MCR_HN120M=read.table('HN120MCR_HN120M_union_06dec_med.csv',sep=",",header=T)
abc=data.frame(HN120MCR_HN120M)
bn = abc[,2:ncol(abc)]
rownames(bn)= abc[,1]
bn$FC[is.na(bn$FC)]=0
bn$FC[is.infinite(bn$FC)]=0

bn$wilk[is.na((bn$wilk))]=1

pdf('HN120MCR_HN120M.pdf',width = 7.5,height=6)
EnhancedVolcano(bn, lab = rownames(bn), x = 'FC', y = 'wilk',xlab = bquote('log2 Fold Change'),ylab = '-log10 p-Value',pCutoff = 0.05,FCcutoff = 1,legendLabels=c('Not Signi.','Fold Change','p-value','p-value & Fold Change'),
                col=c('#ab34eb', 'green', 'orange1', 'red'),colAlpha = 1,pointSize = 4.0,
                labSize = 6.0,title = 'HN120MCR Vs. HN120M Volcano Plot',gridlines.major = FALSE,
                gridlines.minor = FALSE)
dev.off()
##########################################################################################################################################################
##########################################################################################################################################################
##########################################################################################################################################################

### HN137MCR Vs HN137M ###
geneset = read.csv('Non_redundant_UnionTADs_for_data_genes.csv',header=F,strip.white=T,stringsAsFactor=F)
mat=read.table('unipath_pathwayscore_with_UnionTADs_5dec.txt',header=T,row.names = 1)
col=read.table('GSE117872_good_Data_cellinfo1.txt')
col1=col[,1]
wilk= matrix(1 , nrow(mat) , 1) ;
FC= matrix(1 , nrow(mat) , 1) ;
pos = which( col1 == 'HN137MCR') ;
pos1 = which(col1 == 'HN137M'); 
for(k in 1:nrow(mat)){
  tri = wilcox.test( as.numeric(mat[k,pos]) , as.numeric(mat[k,pos1]) )
  wilk[k] = tri$p.value 
  FC[k] = median(as.numeric(mat[k,pos])) - median(as.numeric(mat[k,pos1])) ;
}
HN137MCR_HN137M_TAD=geneset[,1:1]
HN137MCR_HN137M_TAD=data.frame(HN137MCR_HN137M_TAD)
HN137MCR_HN137M_TAD1=paste0("TAD_",1:nrow(HN137MCR_HN137M_TAD))
HN137MCR_HN137M_geneset=geneset[,2:ncol(geneset)]
HN137MCR_HN137M=cbind(HN137MCR_HN137M_TAD1,HN137MCR_HN137M_TAD,wilk,FC,HN137MCR_HN137M_geneset)
colnames(HN137MCR_HN137M)[1]="TADs_number"
colnames(HN137MCR_HN137M)[2]="HN137P_EP_HN137PCR_TADs"
write.table(HN137MCR_HN137M, file = "HN137MCR_HN137M_union.csv", sep=",", row.names = F , col.names = T, quote=F) ;

###### volcano #######
library(EnhancedVolcano)
HN137MCR_HN137M=read.table('HN137P_EP_HN137PCR_union_06dec_med.csv',sep=",",header=T)
abc=data.frame(HN137MCR_HN137M)
bn = abc[,2:ncol(abc)]
rownames(bn)= abc[,1]
bn$FC[is.na(bn$FC)]=0
bn$FC[is.infinite(bn$FC)]=0

bn$wilk[is.na((bn$wilk))]=1

pdf('HN137MCR_HN137M.pdf',width = 7.5,height=6)
EnhancedVolcano(bn, lab = rownames(bn), x = 'FC', y = 'wilk',xlab = bquote('log2 Fold Change'),ylab = '-log10 p-Value',pCutoff = 0.05,FCcutoff = 1,legendLabels=c('Not Signi.','Fold Change','p-value','p-value & Fold Change'),
                col=c('#ab34eb', 'green', 'orange1', 'red'),colAlpha = 1,pointSize = 4.0,
                labSize = 6.0,title = 'HN137MCR Vs. HN137M Volcano Plot',gridlines.major = FALSE,
                gridlines.minor = FALSE)
dev.off()


##########################################################################################################################################################
##########################################################################################################################################################
##########################################################################################################################################################

### HN137PCR Vs HN137P ###
geneset = read.csv('Non_redundant_UnionTADs_for_data_genes.csv',header=F,strip.white=T,stringsAsFactor=F)
mat=read.table('unipath_pathwayscore_with_UnionTADs_5dec.txt',header=T,row.names = 1)
col=read.table('GSE117872_good_Data_cellinfo1.txt')
col1=col[,1]
wilk= matrix(1 , nrow(mat) , 1) ;
FC= matrix(1 , nrow(mat) , 1) ;
pos = which( col1 == 'HN137PCR') ;
pos1 = which(col1 == 'HN137P_EP'); 
for(k in 1:nrow(mat)){
  tri = wilcox.test( as.numeric(mat[k,pos]) , as.numeric(mat[k,pos1]) )
  wilk[k] = tri$p.value 
  FC[k] = median(as.numeric(mat[k,pos])) - median(as.numeric(mat[k,pos1])) ;
}
HN137PCR_HN137P_TAD=geneset[,1:1]
HN137PCR_HN137P_TAD=data.frame(HN137PCR_HN137P_TAD)
HN137PCR_HN137P_TAD1=paste0("TAD_",1:nrow(HN137PCR_HN137P_TAD))
HN137PCR_HN137P_geneset=geneset[,2:ncol(geneset)]
HN137PCR_HN137P=cbind(HN137PCR_HN137P_TAD1,HN137PCR_HN137P_TAD,wilk,FC,HN137PCR_HN137P_geneset)
colnames(HN137PCR_HN137P)[1]="TADs_number"
colnames(HN137PCR_HN137P)[2]="HN137PCR_HN137P_TADs"
write.table(HN137PCR_HN137P, file = "HN137PCR_HN137P_union.csv", sep=",", row.names = F , col.names = T, quote=F) ;

###### volcano #######
library(EnhancedVolcano)
HN137PCR_HN137P=read.table('HN137PCR_HN137P_union_06dec_med.csv',sep=",",header=T)
abc=data.frame(HN137PCR_HN137P)
bn = abc[,2:ncol(abc)]
rownames(bn)= abc[,1]
bn$FC[is.na(bn$FC)]=0
bn$FC[is.infinite(bn$FC)]=0

bn$wilk[is.na((bn$wilk))]=1

pdf('HN137PCR_HN137P.pdf',width = 7.5,height=6)
EnhancedVolcano(bn, lab = rownames(bn), x = 'FC', y = 'wilk',xlab = bquote('log2 Fold Change'),ylab = '-log10 p-Value',pCutoff = 0.05,FCcutoff = 1,legendLabels=c('Not Signi.','Fold Change','p-value','p-value & Fold Change'),
                col=c('#ab34eb', 'green', 'orange1', 'red'),colAlpha = 1,pointSize = 4.0,
                labSize = 6.0,title = 'HN137PCR Vs. HN137P Volcano Plot',gridlines.major = FALSE,
                gridlines.minor = FALSE)
dev.off()

##########################################################################################################################################################
##########################################################################################################################################################
##########################################################################################################################################################
