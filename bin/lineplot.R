# lineplot.R
# creates a line graph of all columns in your input table
# Invoke % R --slave --args [infile] < lineplot.R

# retrieve the args
Args<-commandArgs();

# read the table specified in the args
sum<-read.table(Args[4]);	

# plot all the data (ddh)
#x<-sum$V1;
#y<-sum$V5;
x<-sum$V1;
y1<-sum$V2;
y2<-sum$V3;
y3<-sum$V4;
y4<-sum$V5;
#y5<-sum$V6;
#y6<-sum$V7;
all = data.frame(y1, y2, y3, y4); #, y5, y6);
#matplot(x, y, type="l", lty=1, ylim=c(0,1.2), xlab= "Reads", ylab= "PRC");
#matplot(x, y, type="l", xlab= "Frame", ylab= "Average Number of Flocks");
matplot(x, all, type="l", lty = 1, lwd = 4, xlab="", ylab="", ylim=c(0,1));
#legend("topright", c("zero", "ten", "twenty", "thirty", "forty", "fifty"), cex=1, lty=1, lwd=4, col=1:6);
#legend("topright", c("none", "3X", "7X", "19X"), cex=1, lty=1, lwd=4, col=1:6);

# plot all the data cfi
#x<-sum$V1;
#sum$V1<-NULL;
#matplot(x, sum, type="l", col=1:8, lty=1, ylim=c(0,1.2), xlab= "Partition Fraction", ylab= "CFI Fraction");
#legend("topright", c("CFI", "Node1", "Node2", "Node3", "Node4", "Node5", "Node6", "Node7", "Node8", "Node9", "Node10", "Node11", "Node12"), cex=0.35, lty=1, col=1:13);
#legend("topright", names(sum), col=1:13, cex=0.35, lty=1);
#dev.off();


#plot(sum$V1, sum$V2, type="l", ylim=c(0,1.2));
#lines(sum$V1, sum$V3, type="l", ylim=c(0,1.2));
#dev.off(); 



#df<-data.frame(x=sum$V1, y=sum$V2);

