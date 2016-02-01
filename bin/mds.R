# kmeans_and_mclust.R
# does multi dimmensional scaling of a distance matrix (dissimilarity matrix)
# Invoke % R --slave --args [infile] < mds.R

# retrieve the args
Args<-commandArgs();

# read the table specified in the args
mds<-read.table(Args[4]);	
mdsprod<-cmdscale(mds, k=2);
write.table(mdsprod, quote=FALSE, sep="\t", col.names=FALSE);
plot(mdsprod);

#png(filename="your/file/location/name.png")
#plot(fit)
#dev.off()				

