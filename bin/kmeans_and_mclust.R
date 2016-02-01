# kmeans_and_mclust.R
# does kmeans clustering and EM clustering on a series of x,y points
# Invoke % R --slave --args [infile] [clusters] < kmeans_and_mclust.R

# load libraries
library(mclust);
library (fpc);

# retrieve the args
Args<-commandArgs();

# read the table specified in the args
ipdr<-read.table(Args[4]);	

# do the kmeans clustering using the specifed number of clusters
(kc<-kmeans(ipdr, Args[5]));

x<-ipdr$V1;
y<-ipdr$V2;

#matplot(x, y, type="p", col=kc$cluster);
#plotcluster (ipdr, kc$cluster);
plot(x, y, col=kc$cluster);

# do the EM clustering using mclust
#ipdrclust<-Mclust(ipdr);
ipdrclust<-Mclust(ipdr, modelNames=c("VII","VEI","VVI","VEE","VVE","VEV","VVV"));
#ipdrclust<-Mclust(ipdr, modelNames=Args[6]);			 
summary(ipdrclust, parameters = TRUE);
plot (ipdrclust, what = "BIC");
plot (ipdrclust, what = "classification");
plot (ipdrclust, what = "uncertainty");
plot (ipdrclust, what = "density");
