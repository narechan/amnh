# kmedoids.R
# does kmedoids clustering given an assumed number of clusters and a distance matrix
# Invoke % R --slave --args [dist_matrix] [clusters_expected] < kmedoids.R

# load libraries
library (cluster);

# retrieve the args
Args<-commandArgs();

# read the table specified in the args
k<-read.table(Args[4]);	
km<-pam(k, k=Args[5], diss=TRUE);
clusplot(pam(k, k=Args[5]));
plot (km, which.plots=2);
write.table (km$clustering, quote=FALSE, sep="\t", col.names=FALSE);
