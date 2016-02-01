# deseq.R

# Invoke % R --slave --args [deseqFile] [query1] [query2] < ipd.R

# retrieve the args
Args<-commandArgs();

library( "DESeq" )
countTable <- read.table( Args[4], header=TRUE, row.names=1)
	   
#####HARD=CODED:#####
conds <- factor( c( "S613", "S613", "S613", "R712", "R712", "R712", "TM", "TM", "TM", "LIAR", "LIAR", "LIAR" ) )

cds <- newCountDataSet( countTable, conds )
cds <- estimateSizeFactors( cds )
sizeFactors( cds )
cds <- estimateDispersions( cds )

res <- nbinomTest( cds, Args[5], Args[6] )
pdf(file="out.pdf") #rename it in the bash loop
plotMA(res)
hist(res$pval, breaks=100, col="skyblue", border="slateblue", main="")
dev.off()
write.csv( res, file="out.csv" ) #rename it in the bash loop


