# mm.R
# does a nonlinear regression using the
# m-menton model to calc Km and Vmax
# Invoke % R --slave --args [infile] < mm.R
#	where infile is the summary file 
#	from the ddh analysis 

# retrieve the args
Args<-commandArgs();

# read the table specified in the args
ddh<-read.table(Args[4]);	

# data frame the columns we want
df<-data.frame(x=ddh$V2, y=ddh$V4);

# fit the mm model using nonlin least sqrs
fit<-nls(y~SSmicmen(x, Vm, K), df);

# summarize
summary(fit);
