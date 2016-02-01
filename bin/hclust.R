# hclust.R
# does hierarchical clustering given an assumed number of clusters, a distance matrix, and a method
# Invoke % R --slave --args [dist_matrix] [clusters_expected] [method] < hclust.R

# retrieve the args
Args<-commandArgs();

# read the table specified in the args
h<-read.table(Args[4]);	
hd = dist(h, method='euclidian');
hm<-hclust(hd, method=Args[6]);
hmcut<-cutree(hm, k=Args[5]);
write.table (hmcut, quote=FALSE, sep="\t", col.names=FALSE);
plot (hm);
#summary(hm);
