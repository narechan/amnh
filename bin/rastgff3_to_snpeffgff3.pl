#!/usr/bin/perl -w

#####SETUP#####

use strict;
use Getopt::Long;

my ($infile);
GetOptions(
	   'i|infile=s'  => \$infile,
	   );

#####MAIN#####

my $counter = 0;
open (I, "$infile");
while (my $line = <I>){
    chomp $line;
    if ($line =~m/^#/){
	print "$line\n";
	next;
    }
    else {
	$counter++;
	my ($refctg, $class, $type, $start, $end, $hold1, $strand, $hold2, $info) = 
	    split (/\t/, $line);
	if ($type eq "CDS"){
	    my ($idst, $descst) = split (/\;/, $info);
	    my ($crap1, $id) = split (/\=/, $idst);
	    my ($crap2, $desc) = split (/\=/, $descst);
	    my ($id1, $id2) = split (/\|/, $id);
	    
	    my $cdsid = $id1 . "|" . "CDS." . $id2;
	    my $mrnaid = $id1 . "|" . "mRNA." . $id2;
	    my $geneid = $id1 . "|" . "gene." . $id2;
	    my $exonid = $id1 . "|" . "exon." . $id2;
	    
	    print "$refctg\t$class\tmRNA\t$start\t$end\t$hold1\t$strand\t$hold2\tID=$mrnaid;Parent=$geneid\n";
	    print "$refctg\t$class\tgene\t$start\t$end\t$hold1\t$strand\t$hold2\tID=$geneid;Name=gene$counter\n";
	    print "$refctg\t$class\texon\t$start\t$end\t$hold1\t$strand\t$hold2\tID=$exonid;Name=exon$counter\n";
	    print "$refctg\t$class\tCDS\t$start\t$end\t$hold1\t$strand\t$hold2\tID=$cdsid;Name=CDS$counter\n";
	}
	else{
	    my ($idst, $descst) = split (/\;/, $info);
            my ($crap1, $id) = split (/\=/, $idst);
            my ($crap2, $desc) = split (/\=/, $descst);
            my ($id1, $id2) = split (/\|/, $id);

#            my $cdsid = $id1 . "|" . "CDS" . $id2;
#            my $mrnaid = $id1 . "|" . "mRNA" . $id2;
            my $geneid = $id1 . "|" . "gene." . $id2;
            my $exonid = $id1 . "|" . "exon." . $id2;

            print "$refctg\t$class\t$type\t$start\t$end\t$hold1\t$strand\t$hold2\tID=$id;Name=$desc\n";
            print "$refctg\t$class\tgene\t$start\t$end\t$hold1\t$strand\t$hold2\tID=$geneid;Name=$desc\n";
            print "$refctg\t$class\texon\t$start\t$end\t$hold1\t$strand\t$hold2\tID=$exonid;Name=$desc\n";
#            print "$refctg\t$class\tCDS\t$start\t$end\t$hold1\t$strand\t$hold2\tID=$cdsid;Name=CDS$counter\n";
	}
    }
}
close (I);

   
