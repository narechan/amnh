import argparse
#from pbcore.io.CmpH5IO import *
from pbcore.io.CmpH5Reader import *
import numpy as np

def getIPD( cmpH5, refIdx, refPos, fwdStrand, k=1, minMapQV=10 ):
    refPos = refPos - 1  # motif.gff is 1-indexed, cmp.h5 0-indexed
    outIPD = []
    for aln in cmpH5.readsInRange( refIdx, refPos - k, refPos + k ):
        if aln.isForwardStrand != fwdStrand: # want complement strand of modified motif
            if aln.MapQV >= minMapQV:
                idx = np.where(aln.referencePositions()==refPos)[0][0]
                slc = slice(idx-k,idx+k+1)
                if '-' not in (aln.reference()[slc] + aln.read()[slc]): # match over window
                    outIPD.append(aln.IPD()[idx])
    return outIPD


parser = argparse.ArgumentParser(description='get IPD ratios for reference position')
parser.add_argument('cmpH5', type=str, help='aligned_reads.cmp.h5 corresponding to motif.gff')
parser.add_argument('refIdx', type=int, help='index of reference contig (1 if single contig)')
parser.add_argument('refPos', type=int, help='position of modified cognate, 4th column in motifs.gff')
parser.add_argument('-f',  dest='fwdStrand', action='store_true', help='us -f flag if + strand in motifs.gff')
parser.add_argument('-r',  dest='fwdStrand', action='store_false', help='us -r flag if - strand in motifs.gff')
parser.add_argument('-k', type=int, default=1, help='min number of bases on each side of modified base which must align in read')
parser.add_argument('-q,--minMapQV', dest='minMapQv',type=int, default=10, help='minimum mapping QV of read')

args = parser.parse_args()
print( getIPD( CmpH5Reader(args.cmpH5),args.refIdx,args.refPos,args.fwdStrand,args.k,args.minMapQv ) )
