#!/usr/bin/perl -w

#####SETUP#####

use strict;
use Getopt::Long;
use Bio::AlignIO;

my ($infile, $finfile, $quality);
GetOptions(
	   'i|infile=s'  => \$infile,
	   'f|finfile=s' => \$finfile,
	   'q|quality=s' => \$quality,
	   );

#`mkdir -p $outdir`;

#####MAIN#####

# readin the file and store data
my $snps = {};
my $refbases = {};
my $annots = {};

my $counter = 0;
open (I, "$infile");
while (my $line = <I>){
    chomp $line;
    next if ($line =~m/Chrom\tQuery/);

    my ($refctg, $query, $refpos, $refbase, $snp, $sq, $stuff) = 
	split (/\t/, $line, 7);

    # plug in the reference chars
#    $snps->{$refpos} = $refbase;
#    $refbases->{$refpos} = $refbase;

    # filters
    next if ($refbase eq "N"); #ambiguous ref
    if ($sq < $quality){ #sq
	$snps->{$refpos} = "?";
	next;
    }
    if ($snp =~m/\,/){ #hets
	$snps->{$refpos} = "?";
	next;
    }
    $counter++;
    
    print "$line\n";
    # plug in the snp chars
    $snps->{$refpos} = $snp;

}
close (I);
print STDERR "$counter\n";

my $lines = 0;
my $question = 0;
my $allqs = 0;
my $refqs = 0;
my $exists = 0;
my $sames = 0;
my $notexists = 0;
my $qryqs = 0;
my $diffs = 0;
open (F, "$finfile");
while (my $line = <F>){
    chomp $line;
    next if ($line =~m/Chrom\tQuery/);

    my ($refctg, $query, $refpos, $refbase, $snp, $sq, $stuff) =
        split (/\t/, $line, 7);
    
    # filters                                                                                             
    my $char;
    next if ($refbase eq "N"); #ambiguous ref                                                               
    if ($sq < $quality){ #sq                                                                              
        $char = "?";
        next;
    }
    if ($snp =~m/\,/){ #hets                                                                              
        $char = "?";
        next;
    }
    $lines++;
    
    if (exists ($snps->{$refpos})){
	$exists++;
	if ($snps->{$refpos} eq $snp){
	    $sames++;
	    if (($snps->{$refpos} eq "?") and ($snp eq "?")){
		$allqs++;
	    }
	    else{
		next;
	    }
	}
	elsif (($snps->{$refpos} eq "?") and ($snp ne "?")){
	    $refqs++;
	}
	elsif (($snps->{$refpos} ne "?") and ($snp eq "?")){
	    $qryqs++
	}
	else {
	    $diffs++;
	}
    }
    else {
	$notexists++;
    }
}
close (F);

print STDERR "$lines\t$exists\t$sames\t$diffs\t$allqs\t$refqs\t$qryqs\t$notexists\n";

=head
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
my $tts = 0;
my $vrsa = 0;
my $vssa = 0;
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
		($vrsa++) if ($snps->{$contig}->{$count}->{'VSSA'} eq $snps->{$contig}->{$count}->{'saureus_USA300_TCH1516'});
		($vssa++) if ($snps->{$contig}->{$count}->{'VRSA'} eq $snps->{$contig}->{$count}->{'saureus_USA300_TCH1516'});
		$tts++;
	    }
	}
	my $table = join "\t", @table;
	print TBL "$table\t$annots->{$contig}->{$count}\n";

	print CAT "\n";
	$start = $end + 1;
    }
}
print STDERR "$tts\t$vrsa\t$vssa\n";
print CAT ";\n";
print CAT "END;\n\n";
print PRT "END;\n";
close (CAT);
close (PRT);
close (TBL);

`cat $outdir/matrix $outdir/partitions > $outdir/nexus`;
=cut
