#!/usr/bin/perl -w

# This program accepts a list of snps from snp_annotator.pl
# and creates a matrix from that data for subsequent phylogenetic analysis.

# NOTE: the number of taxa in the matrix is the number of queries in 
#       your list plus 1, for the reference.

# -i is the input snp table
# -l is the query list
# -o is the outdir
# -q is the snp quality above which snps are accepted
# -r is the reference (string as in the taxon name)

#####SETUP#####

use strict;
use Getopt::Long;
use Bio::AlignIO;

my ($infile, $list, $outdir, $quality, $ref);
GetOptions(
	   'i|infile=s'  => \$infile,
	   'l|list=s'    => \$list,
	   'o|outdir=s'  => \$outdir,
	   'q|quality=s' => \$quality,
	   'r|reference=s' => \$ref,
	   );

`mkdir -p $outdir`;

#####MAIN#####

# readin the quey list
my @queries;
open (L, "$list");
while (my $line = <L>){
    chomp $line;
    push (@queries, $line);
}
close (L);

# readin the file and store data
my $snps = {};
my $refbases = {};
my $annots = {};
#my $tvrsa = 0;
#my $tvssa = 0;
open (I, "$infile");
while (my $line = <I>){
    chomp $line;
    next if ($line =~m/Chrom\tQuery/);

    my ($refctg, $query, $refpos, $refbase, $snp, $sq, $stuff) = 
	split (/\t/, $line, 7);

    # plug in the reference chars
    $snps->{$refctg}->{$refpos}->{$ref}  = $refbase;
    $refbases->{$refctg}->{$refpos} = $refbase;

    # store the annotation
    $annots->{$refctg}->{$refpos} = $stuff;
    
    # filters
    next if ($refbase eq "N"); #ambiguous ref
    if ($sq < $quality){ #sq
	$snps->{$refctg}->{$refpos}->{$query} = "?";
	next;
    }
    if ($snp =~m/\,/){ #hets
	$snps->{$refctg}->{$refpos}->{$query} = "?";
	next;
    }
#    if ($query eq "VRSA"){
#	$tvrsa++;
#	print "$line\n";
#    }
#    if ($query eq "VSSA"){
#	$tvssa++;
#	print "$line\n";
#    }

    # plug in the snp chars
    $snps->{$refctg}->{$refpos}->{$query} = $snp;
}
close (I);

#print STDERR "$tvrsa\t$tvssa\n";

# post process datastruc to fill in ref alleles
# for queries with snp gaps
foreach my $ctg (keys %$snps){
    foreach my $pos (keys %{$snps->{$ctg}}){
	foreach my $query (@queries){
	    if (exists($snps->{$ctg}->{$pos}->{$query})){
		next;
	    }
	    else {
		$snps->{$ctg}->{$pos}->{$query} = $refbases->{$ctg}->{$pos};
	    }
	}
    }
}

=head
# vrsa / vssa specific code
foreach my $pos (sort {$a <=> $b} keys %{$snps->{'saureus_USA300_TCH1516'}}){
#    if ( ($snps->{'saureus_USA300_TCH1516'}->{$pos}->{'VRSA'} eq "?") or ($snps->{'saureus_USA300_TCH1516'}->{$pos}->{'VSSA'} eq "?") ){
#	next;
#    }
    
#    if ($snps->{'saureus_USA300_TCH1516'}->{$pos}->{'VRSA'} eq "?"){
#	next;
#    }
#    print "$pos\t$snps->{'saureus_USA300_TCH1516'}->{$pos}->{'VRSA'}\n";

#    if ($snps->{'saureus_USA300_TCH1516'}->{$pos}->{'VSSA'} eq "?"){                                         
#       next;                                                                                              
#    }  
#    print "$pos\t$snps->{'saureus_USA300_TCH1516'}->{$pos}->{'VSSA'}\n";

#    print "$pos\t$snps->{'saureus_USA300_TCH1516'}->{$pos}->{'VRSA'}\t$snps->{'saureus_USA300_TCH1516'}->{$pos}->{'VSSA'}\t$snps->{'saureus_USA300_TCH1516'}->{$pos}->{'saureus_USA300_TCH1516'}\n";
#    if ( ($snps->{'saureus_USA300_TCH1516'}->{$pos}->{'VRSA'} ne $snps->{'saureus_USA300_TCH1516'}->{$pos}->{'saureus_USA300_TCH1516'}) and ($snps->{'saureus_USA300_TCH1516'}->{$pos}->{'VSSA'} eq $snps->{'saureus_USA300_TCH1516'}->{$pos}->{'saureus_USA300_TCH1516'}) ){
#	print "$pos\t$snps->{'saureus_USA300_TCH1516'}->{$pos}->{'VRSA'}\t$snps->{'saureus_USA300_TCH1516'}->{$pos}->{'VSSA'}\t$snps->{'saureus_USA300_TCH1516'}->{$pos}->{'saureus_USA300_TCH1516'}\n";
#    }
#    if ( ($snps->{'saureus_USA300_TCH1516'}->{$pos}->{'VRSA'} ne $snps->{'saureus_USA300_TCH1516'}->{$pos}->{'saureus_USA300_TCH1516'}) and ($snps->{'saureus_USA300_TCH1516'}->{$pos}->{'VSSA'} eq "?") ){
#        print "$pos\t$snps->{'saureus_USA300_TCH1516'}->{$pos}->{'VRSA'}\t$snps->{'saureus_USA300_TCH1516'}->{$pos}->{'VSSA'}\t$snps->{'saureus_USA300_TCH1516'}->{$pos}->{'saureus_USA300_TCH1516'}\n";
#    }
#    if ( ($snps->{'saureus_USA300_TCH1516'}->{$pos}->{'VSSA'} ne $snps->{'saureus_USA300_TCH1516'}->{$pos}->{'saureus_USA300_TCH1516'}) and ($snps->{'saureus_USA300_TCH1516'}->{$pos}->{'VRSA'} eq $snps->{'saureus_USA300_TCH1516'}->{$pos}->{'saureus_USA300_TCH1516'}) ){
#        print "$pos\t$snps->{'saureus_USA300_TCH1516'}->{$pos}->{'VRSA'}\t$snps->{'saureus_USA300_TCH1516'}->{$pos}->{'VSSA'}\t$snps->{'saureus_USA300_TCH1516'}->{$pos}->{'saureus_USA300_TCH1516'}\n";
#    }
    if ( ($snps->{'saureus_USA300_TCH1516'}->{$pos}->{'VSSA'} ne $snps->{'saureus_USA300_TCH1516'}->{$pos}->{'saureus_USA300_TCH1516'}) and ($snps->{'saureus_USA300_TCH1516'}->{$pos}->{'VRSA'} eq "?") ){
        print "$pos\t$snps->{'saureus_USA300_TCH1516'}->{$pos}->{'VRSA'}\t$snps->{'saureus_USA300_TCH1516'}->{$pos}->{'VSSA'}\t$snps->{'saureus_USA300_TCH1516'}->{$pos}->{'saureus_USA300_TCH1516'}\n";
    }
}
=cut
# get the pairwise counts of diffs between queries
# and print to stderr
my $pairwise = {};
my $qrycopy = {};
foreach my $qry (@queries){
    $qrycopy->{$qry} = 1;
}

foreach my $ctg (keys %$snps){
    foreach my $q1 (sort keys %$qrycopy){
	foreach my $q2 (sort keys %$qrycopy){
	    next if ($q1 eq $q2);
	    
	    my $pair = $q1 . "-" . $q2;
	    foreach my $pos (keys %{$snps->{$ctg}}){
		if ( ($snps->{$ctg}->{$pos}->{$q1} eq "?") or ($snps->{$ctg}->{$pos}->{$q2} eq "?") ){
		    next;
		}
		elsif ($snps->{$ctg}->{$pos}->{$q1} ne $snps->{$ctg}->{$pos}->{$q2}){
                    $pairwise->{$pair}++;
#		    print STDERR "$ctg\t$pos\n";
                }
                else {
                    next;
                }
	    }
	}
	delete $qrycopy->{$q1};
    }
}

foreach my $pair (sort keys %$pairwise){
    print STDERR "$pair\t$pairwise->{$pair}\n";
}


# get some matrix data
my $taxa = @queries;
$taxa++;
my $alnlentot = 0;
foreach my $c (keys %$snps){
    foreach my $s (keys %{$snps->{$c}}){
	$alnlentot++;
    }
}

# sort print the matrix and charpars                                   
open (CAT, ">$outdir/matrix");
open (PRT, ">$outdir/partitions");
open (TBL, ">$outdir/table");
print CAT "#NEXUS\n";
print CAT "BEGIN DATA;\n";
print CAT "DIMENSIONS NTAX=$taxa NCHAR=$alnlentot;\n";
print CAT "FORMAT INTERLEAVE SYMBOLS=\"ABCDEFGHIKLMNPQRSTUVWXYZ\" DATATYPE=PROTEIN MISSING=? GAP=-;\n";
print CAT "MATRIX\n";
print PRT "BEGIN SETS;\n";

# print table header
foreach my $contig (keys %$snps){
    print TBL "Contig\tSNP\t";
    foreach my $count (sort {$a <=> $b} keys %{$snps->{$contig}}){
	my @sp;
	foreach my $sp (sort keys %{ $snps->{$contig}->{$count} }){
	    push (@sp, $sp);
	}
	my $spstring = join "\t", @sp;
	print TBL "$spstring\n";
	last;
    }
    last;
}

my $start = 1;
my $end;
#my $tts = 0;
#my $vrsa = 0;
#my $vssa = 0;
foreach my $contig (keys %$snps){
    foreach my $count (sort {$a <=> $b} keys %{$snps->{$contig}}){
	$end = $start;
	print CAT "[Partition ctg$contig snp$count]\n";

	my @table;
#	print TBL "ctg$contig\tsnp$count\t";
#	print TBL "$contig\t$count\t"; 
	push (@table, $contig, $count);
	
	my $charset = "ctg" . $contig . "snp" . $count;
	print PRT "CHARSET $charset=$start-$end;\n";
    
	my @qrysnps;
	foreach my $sp (sort keys %{ $snps->{$contig}->{$count} }){
	    print CAT "$sp\t$snps->{$contig}->{$count}->{$sp}\n"; 
	    push (@table, $snps->{$contig}->{$count}->{$sp});
	    push (@qrysnps, $snps->{$contig}->{$count}->{$sp}) unless ($sp eq $ref);
	}
	
	# check to see if all the qry snps are the same
	if (keys %{{ map {$_, 1} @qrysnps }} == 1){
	    push (@table, "N");
	}
	else {
	    my $signal = 0;
	    foreach my $qsnp (@qrysnps){
		($signal++) if ($qsnp eq "?");
	    }
	    if ($signal > 0){
		push (@table, "*");
	    }
	    else {
		push (@table, "+");
#		($vrsa++) if ($snps->{$contig}->{$count}->{'VSSA'} eq $snps->{$contig}->{$count}->{'saureus_USA300_TCH1516'});
#		($vssa++) if ($snps->{$contig}->{$count}->{'VRSA'} eq $snps->{$contig}->{$count}->{'saureus_USA300_TCH1516'});
#		$tts++;
	    }
	}
	my $table = join "\t", @table;
	print TBL "$table\t$annots->{$contig}->{$count}\n";

	print CAT "\n";
	$start = $end + 1;
    }
}
#print STDERR "$tts\t$vrsa\t$vssa\n";
print CAT ";\n";
print CAT "END;\n\n";
print PRT "END;\n";
close (CAT);
close (PRT);
close (TBL);

`cat $outdir/matrix $outdir/partitions > $outdir/nexus`;
