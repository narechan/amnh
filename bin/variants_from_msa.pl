#!/usr/bin/perl -w

# -i is the alignment infile
# -f is the input format

#####SETUP#####

use strict;
use Getopt::Long;
use Bio::AlignIO;

my ($infile, $format, $refacc);
GetOptions(
	   'i|infile=s'      => \$infile,
	   'f|format=s'      => \$format,
           'r|refacc=s'      => \$refacc,
	   );

#####MAIN#####

my $alnin = Bio::AlignIO->new(-file   => "$infile",
			      -format => "$format");

# get aln data
my $aln = {};
my $ref = {};
my $alnobj = $alnin->next_aln();
foreach my $seq ($alnobj->each_seq){
    my $id        = $seq->display_id;
    my $sequence = $seq->seq;
    my @sequence = split (//, $sequence);
    my $alncounter = 0;
    foreach my $base (@sequence){
	$alncounter++;
	if ($id eq $refacc){
	    $ref->{$alncounter} = $base;
	}
	else{
	    $aln->{$alncounter}->{$id} = $base;
	}
    }
}

# print the header
my @header;
push (@header, $refacc);
foreach my $alnpos (sort {$a <=> $b} keys %$aln){
    foreach my $id (sort keys %{$aln->{$alnpos}}){
	push (@header, $id);
    }
    last;
}
my $header = join "\t", @header;
print "AlnPos\tRefPos\t$header\n";


my $refpos = 0;
foreach my $alnpos (sort {$a <=> $b} keys %$aln){
    my @bases;
    my $refbase = $ref->{$alnpos};

    # iterate refpos unless gap in ref                                                                        
    unless ($refbase eq "-"){
	$refpos++;
    }

    push (@bases, $refbase);

    foreach my $id (sort keys %{$aln->{$alnpos}}){
	push (@bases, $aln->{$alnpos}->{$id});
    }
    my $bases = join "\t", @bases;
    
    # check to see if all snps are the same
    if (keys %{{ map {$_, 1} @bases }} == 1){
	next;
    }
    else{
	my $ref = shift @bases;
	print "$alnpos\t$refpos\t$ref\t";
	
	my @string;
	foreach my $base (@bases){
	    if ($base eq $ref){
		push (@string, ".");
	    }
	    else{
		push (@string, $base);
	    }
	}
	my $string = join "\t", @string;
	print "$string\n";
    }
}
