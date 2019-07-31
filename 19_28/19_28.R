library(tidyverse)
dir<-file.path("C:/Users/My computer/Dropbox/CSHL/R/qpcr_analysis")
fc_path<-file.path(dir,"19_28","19_28_fc.csv")
data<-read_csv(fc_path) # read fc data into tibble
# remove some outliers
which(data$FC<0.3)
data<-data[-c(21,93),]
gm_rpl32 <- 
# calculate group means
gm<-aggregate(x= data$FC, by = list(tx = data$Tx, ref_gene = data$Ref_gene), mean)

#separate by reference gene and exclude control groups
gm_rpl32<-filter(gm,ref_gene == "RPL32", !tx %in%c("NT","scramble","488 507 526 120 nM (60 pmol)"))
gm_hprt<-filter(gm,ref_gene == "HPRT", !tx %in%c("NT","scramble","488 507 526 120 nM (60 pmol)"))

# order the tx by ascending order of fc
txorder1<-gm_rpl32[order(gm_rpl32$x),"tx"]
txorder2<-gm_rpl32[order(gm_hprt$x),"tx"]

# sort data by ref_gene
data1<-filter(data, Ref_gene == "RPL32")
data2<-filter(data, Ref_gene == "HPRT")
data3<-filter(data, Ref_gene == "RPL32")
gm_rpl32
# plot by the rank of the tx
data1$Tx<-factor(data1$Tx, levels = c("NT","scramble","488 507 526 120 nM (60 pmol)",txorder1))
data2$Tx<-factor(data2$Tx, levels = c("NT","scramble","488 507 526 120 nM (60 pmol)",txorder2))
data3$Tx<-factor(data1$Tx, levels = c(txorder1[21:1], "488 507 526 120 nM (60 pmol)","scramble","NT"))
library(ggpubr)

p1<-ggplot(data1,aes(x=data1$Tx, y = data1$FC)) + geom_boxplot(color = "red")+ geom_point() + 
  theme(axis.text.x = element_text(hjust = 1, angle = 45)) +
  xlab("treatment") + ylab("CFTR fold change") +  ylim(0.1,4.5) +
  ggtitle("CFTR/RPL32 vs NT")

p2<-ggplot(data2,aes(x=data2$Tx, y = data2$FC)) + geom_boxplot(color = "blue")+ geom_point() + ylim(0.1,4.5)+
  theme(axis.text.x = element_text(hjust = 1, angle = 45)) +
  xlab("treatment") + 
  ylab("CFTR fold change") + 
  ggtitle("CFTR/HPRT vs NT")

ggarrange(p1,p2)
# vs. NT
summary(lm(FC ~ Tx, data1))
summary(lm(FC ~ Tx, data2))
# vs. 20-20-20 pmol
summary(lm(FC ~ Tx, data3))


lookup<-tibble(txorder1, index1 = c(1:21), txorder2, index2 = match(txorder2,txorder1))
cor(lookup$index1,lookup$index2, method = "spearman")
# [1] 0.9467532, this means the rank of fc calucated with RPL32 or HPRT as ref gene correlates very well. 

bin_path<-file.path(dir,"19_28","bin19_28.csv")
bindata<-read_csv(bin_path)
#clustering and heatmaps https://earlglynn.github.io/RNotes/package/gplots/heatmap2.html
library(gplots)
binlabels<-apply(binmat[,1:3],1,paste,collapse="-")
binmat<-matrix(unlist(bindata[,4:7]),36,4, dimnames=list(binlabels, c("aso24","aso25","aso26","fc")))
binmat[,1:3]<-log(binmat[,1:3],3)
binmat2<-binmat[order(binmat[,"fc"], decreasing = TRUE),]
heatmap.2(binmat, trace="none", cexRow = 0.8, cexCol = 1.5, col=bluered, Colv=FALSE, reorderfun = function (d, w) reorder(d, w, agglo.FUN = median) )
heatmap.2(binmat, trace="none", cexRow = 0.8, cexCol = 1.5, col=bluered, Colv=FALSE)
heatmap.2(binmat2, trace="none", cexRow = 0.5, cexCol = 1.5, col=bluered, Colv=FALSE, Rowv=FALSE)
heatmap.2(binmat2, trace="none", cexRow = 0.5, cexCol = 1.5, col=bluered)
# doseage of which ASO correlates the most with the fold change?
for(n in 1:3){
  print(cor(binmat2[,n],binmat2[,4]))
}

# analyze similar experiment: exp18-65
fc_path2<-file.path(dir,"18_65","18_65_fc.csv")
fc65<-read_csv(fc_path2)

fc65<-aggregate(x = fc65$fc, 
          by = list(Tx = fc65$Tx, 
                    aso24 = fc65$aso24, 
                    aso25 = fc65$aso25,
                    aso26 = fc65$aso26), mean, na.rm =TRUE)

fc65<-matrix(unlist(fc65[,2:5]), ncol = 4, 
             dimnames = list(
               fc65$Tx,
               c(colnames(fc65)[2:4],"fc")
               )
)

fc65<-fc65[order(fc65[,"fc"], decreasing = TRUE),]
fc65[,1:3]<-log(fc65[,1:3]+1,3)
heatmap.2(fc65, trace="none", cexRow = 0.5, cexCol = 1.5, Colv=FALSE, Rowv=FALSE)
heatmap.2(fc65, trace="none", cexRow = 0.5, col = bluered, cexCol = 1.5, Colv=FALSE)

# doseage of which ASO correlates the most with the fold change?
for(n in 1:3){
  print(cor(fc65[,n],fc65[,4]))
}
