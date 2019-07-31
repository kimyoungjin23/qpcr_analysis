setwd("C:/Users/My computer/Dropbox/CSHL/qPCR/Data")


qpcr.df <- X20190514_exp_19_20_genes
expression<-qpcr.df['expression_vs_rpl32']
expression<-sapply(expression, as.numeric)
qpcr.df['expression_vs_rpl32']<-expression

#How do you specifically order ggplot2 x axis instead of alphabetical order? 
#data$Treatment <- factor(data$Treatment, levels=unique(data$Treatment))

qpcr.df$treatment <- factor(qpcr.df$treatment, level = unique(qpcr.df$treatment))
plot<-ggplot(qpcr.df)+geom_boxplot(aes(x= qpcr.df$treatment, y= qpcr.df$expression_vs_rpl32, color =qpcr.df$Target_gene))
plot + theme(axis.text.x = element_text(face="bold", color="#993333", size=7, angle=45)) + scale_y_continuous(breaks=c(1:6))

plot<-ggplot(qpcr.df)+geom_boxplot(aes(x= qpcr.df$treatment, y= qpcr.df$expression_vs_rpl32, color =qpcr.df$Treatment_length
                                       ))

qpcr_plot <- function (qpcr.data) {
  ggplot(qpcr.data)+geom_boxplot(aes(x= qpcr.data$treatment, y= qpcr.data$expression_vs_rpl32, color =qpcr.data$Treatment_length))
}


qpcr.df.cftr<-filter(qpcr.df, qpcr.df$Target_gene == "CFTR")
cftr.plot<- ggplot(qpcr.df.cftr)+geom_boxplot(aes(x= qpcr.df.cftr$treatment, y= qpcr.df.cftr$expression_vs_rpl32, color =qpcr.df.cftr$Treatment_length))
cftr.plot+theme(axis.text.x = element_text(face="bold", color="#993333", size=7, angle=45)) + scale_y_continuous(breaks=c(1:6))

qpcr.df.c1orf37<-filter(qpcr.df, qpcr.df$Target_gene == "C1orf37")
c1orf37.plot <- ggplot(qpcr.df.c1orf37)+geom_boxplot(aes(x= qpcr.df.c1orf37$treatment, y= qpcr.df.c1orf37$expression_vs_rpl32, color =qpcr.df.c1orf37$Treatment_length))
c1orf37.plot + theme(axis.text.x = element_text(face="bold", color="#993333", size=7, angle=45)) + scale_y_continuous(breaks=c(1:8))

qpcr.df.eif4a2<-filter(qpcr.df, qpcr.df$Target_gene == "eIF4A2")
eif4a2.plot<-ggplot(qpcr.df.eif4a2)+ geom_boxplot(aes(x= qpcr.df.eif4a2$treatment, y= qpcr.df.eif4a2$expression_vs_rpl32, color =qpcr.df.eif4a2$Treatment_length))
eif4a2.plot + theme(axis.text.x = element_text(face="bold", color="#993333", size=7, angle=45)) + scale_y_continuous(breaks=c(1:8))

qpcr.df.gadd45a<-filter(qpcr.df, qpcr.df$Target_gene == "GADD45A")
gadd42a.plot<-qpcr_plot(qpcr.df.gadd45a)
gadd42a.plot + theme(axis.text.x = element_text(face="bold", color="#993333", size=7, angle=45)) + scale_y_continuous(breaks=c(1:8))

qpcr.df.rpl12<-filter(qpcr.df, qpcr.df$Target_gene == "RPL12")
rpl12.plot <- qpcr_plot(qpcr.df.rpl12)
rpl12.plot + theme(axis.text.x = element_text(face="bold", color="#993333", size=7, angle=45)) + scale_y_continuous(breaks=c(1:8))

qpcr.df.sf3b1<-filter(qpcr.df, qpcr.df$Target_gene == "SF3B1")
sf3b1.plot <- qpcr_plot(qpcr.df.sf3b1)
sf3b1.plot +theme(axis.text.x = element_text(face="bold", color="#993333", size=7, angle=45)) + scale_y_continuous(breaks=c(1:8))

qpcr.df.srsf2<-filter(qpcr.df, qpcr.df$Target_gene == "SRSF2")
srsf2.plot <- qpcr_plot(qpcr.df.srsf2)
srsf2.plot +theme(axis.text.x = element_text(face="bold", color="#993333", size=7, angle=45)) + scale_y_continuous(breaks=c(1:8))

qpcr.df.srsf4<-filter(qpcr.df, qpcr.df$Target_gene == "SRSF4")
srsf4.plot <- qpcr_plot(qpcr.df.srsf4)
srsf4.plot +theme(axis.text.x = element_text(face="bold", color="#993333", size=7, angle=45)) + scale_y_continuous(breaks=c(1:8))

