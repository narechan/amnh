# ipd.R
# gets the ipd values for all reads overlapping a given position (forward and reverse)
# from a pacbio cmp.h5 file.
# Invoke % R --slave --args [position] [strand {0 or 1}] [in *.cmp.h5 file] [outflie] < ipd.R

# requires the pbh5 package
require (pbh5);

# retrieve the args
Args<-commandArgs();

# load data
cmpH5 <- PacBioCmpH5(Args[6])

# get the IPDs
tposIPD <- getByTemplatePosition(cmpH5, idx = getReadsInRange(cmpH5, 1, Args[4], Args[4]), f = getIPD)

# filter the IPDs
tposIPD <- subset(tposIPD, read == ref)
tposIPD <- subset(tposIPD, strand == Args[5] & position == Args[4])	

# print it out
write.table(tposIPD, file =Args[7], sep ="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
