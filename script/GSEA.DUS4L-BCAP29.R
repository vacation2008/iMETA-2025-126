


#���ð�
library(limma)
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)

expFile="symbol.txt"         

gmtkegg="c2.all.v2023.1.Hs.symbols.gmt"    

setwd("C:\\Users\\cikeli\\Desktop\\0213\\gsea\\DUS4L-BCAP29.GSEA")     


gmtkegg=read.gmt(gmtkegg)


rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)



risk <- read.delim("C:/Users/cikeli/Desktop/0213/gsea/DUS4L-BCAP29.GSEA/risklandscap.DUS4L-BCAP29.tsv")
rownames(risk)=risk[,1]
risk=risk[,2:ncol(risk)]


data=data[rowMeans(data)>0.5,]
dataL=data[,row.names(risk[risk[,"risk"]=="low",])]
dataH=data[,row.names(risk[risk[,"risk"]=="high",])]
meanL=rowMeans(dataL)
meanH=rowMeans(dataH)
meanL[meanL<0.00001]=0.00001
meanH[meanH<0.00001]=0.00001
logFC=log2(meanH)-log2(meanL)
logFC=sort(logFC,decreasing=T)
genes=names(logFC)


kk=GSEA(logFC, TERM2GENE=gmtkegg, pvalueCutoff=1, minGSSize=15, maxGSSize=500)
kkTab=as.data.frame(kk)
kkTab=kkTab[kkTab$pvalue<0.05,]    

write.table(kkTab,file="DUS4L-BCAP29.txt",sep="\t",quote=F,row.names = F)


termNum=1     
kkUp=kkTab[kkTab$NES>=1,]
if(nrow(kkUp)>=termNum){
  rows <- c("NELSON_RESPONSE_TO_ANDROGEN_UP")
  selected_rows <- kkUp[row.names(kkUp) %in% rows, ]
  showTerm=row.names(selected_rows)[1:1]     
  gseaplot=gseaplot2(kk, showTerm, base_size=8, title="NELSON_RESPONSE_TO_ANDROGEN_UP")
  pdf(file="GSEA.DUS4L-BCAP29.pdf", width=7, height=5.5)
  print(gseaplot)
  dev.off()
}


