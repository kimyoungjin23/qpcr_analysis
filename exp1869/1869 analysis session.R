melt(exp1859)

lm(exp1859_CFTR~1+2+3+4+5+6+7, data.frame(exp18))

qpcr1869<-read.csv("1869ddCt.csv",header=F) #read qpcr data and save
qpcr1869<-na.omit(qpcr1869) #remove NA data
qpcr1869.rx<-lapply(qpcr1869[1], unlist) # unlist treatment groups
qpcr1869.foldchange<-lapply(qpcr1869[2], unlist) #makes fold change as numeric 
qpcr1869.foldchange<-lapply(qpcr1869[2], as.numeric) #makes fold change as numeric 

## generate data

x = factor(sample(letters[1:5],100, replace=TRUE))

print(levels(x))  ## This will show the levels of x are "Levels: a b c d e"

## To reorder the levels:
## note, if x is not a factor use levels(factor(x))
x = factor(x,levels(x)[c(4,5,1:3)])

print(levels(x))  ## Now "Levels: d e a b c"
