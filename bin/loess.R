# loess.R
# does a smooth curve fitting using local regression
# Invoke % R --slave --args [span] [infile] [outflie] < loess.R
#	where infile is the summary file 
#	from the cfi analysis

# retrieve the args
Args<-commandArgs();

# read the table specified in the args
loe<-read.table(Args[5]);	

# data frame the columns we want
df<-data.frame(x=loe$V1, y=loe$V2);

# fit the loess model
sp<-as.numeric(Args[4]);
cfiloess<-loess(y ~ x, span=sp, data=df);

# predict with df xvals by default
cfi<-predict(cfiloess);

# create dataframe to output
dfout<-data.frame(x=df$x, y=cfi);

# print it out
write.table(dfout, file =Args[6], sep ="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
#print(dfout);