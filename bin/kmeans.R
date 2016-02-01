# kmeans.R
# does kmeans clustering on a series of x,y points
# Invoke % R --slave --args [infile] [clusters] < kmeans.R

# retrieve the args
Args<-commandArgs();

# read the table specified in the args
ipdr<-read.table(Args[4]);	

# do the kmeans clustering using the specifed number of clusters
kc<-kmeans(ipdr, Args[5]);
#kc$totss;
kc$withinss;
#kc$tot.withinss;
