#!/usr/bin/perl -w

# -m is the matrix
# -a is comma delimited string of first group (ingroup)
# -b is comma delimited string of the second group (outgroup)
# -o is the outdir

#####SETUP#####

use strict;
use Getopt::Long;
use Bio::AlignIO;
use Bio::TreeIO;
use TreeSupports;

my ($matrixfile, $group1, $group2, $outdir);
GetOptions(
	   'm|matrix=s'    => \$matrixfile,
	   'a|group1=s'    => \$group1,
	   'b|group2=s'    => \$group2,
           'o|outdir=s'    => \$outdir,
	   );

`mkdir -p $outdir`;

#####MAIN######

# store the intaxa and outtaxa
my @intax  = split (/\,/, $group1);
my @outtax = split (/\,/, $group2);

my $intax = {};
foreach my $intaxi (@intax){
    $intax->{$intaxi} = 1;
}

my $outtax = {};
foreach my $outtaxi (@outtax){
    $outtax->{$outtaxi} = 1;
}

# store the matrix data
my $supportobj = TreeSupports->new;
$supportobj->load_aln ($matrixfile);
my $charsets = $supportobj->get_charsets;
my $chars    = $supportobj->get_nchar;

my $alnin = Bio::AlignIO->new(-file   => $matrixfile,
                              -format => 'nexus');
my $aln = $alnin->next_aln();

my $alndata = {};
my $taxa = {};
foreach my $seq ($aln->each_seq){
    my $id        = $seq->display_id;
    $taxa->{$id} = 1;
#    print STDERR "$id\n";
    
    foreach my $charset (sort keys %$charsets){

        my $coords = $charsets->{$charset};
        my ($start, $end) = split (/\-/, $coords);

        my $partition = $seq->subseq($start, $end);
	$alndata->{$charset}->{$id} = $partition;

    }
}

# cycle through the snps and determine which are 
# diagnostic for distingushing the taxa groups
# also create the all.snp report
open (R, ">$outdir/diag.snps");
open (A, ">$outdir/all.snps");
my @taxaheader = sort keys %$taxa;
my $taxaheader = join "\t", @taxaheader;
print A "SNP\t$taxaheader\n";

foreach my $snp (keys %$alndata){
#    print STDERR "$snp\n";
    $snp =~m/(\d+)/;
    my $refpos = $1;

    # print the all.snps info
    my @bases;
    foreach my $taxon (sort keys %{$alndata->{$snp}}){
	push (@bases, $alndata->{$snp}->{$taxon});
    }
    my $bases = join "\t", @bases;
    print A "$snp\t$bases\n";
    
    my @nodein;
    my @nodeout;
    my @nodeouttaxa;
    my @nodeintaxa;
    foreach my $taxon (keys %$taxa){
	
	# populate the outgroup
	if (exists ($outtax->{$taxon})){
	    push (@nodeout, $alndata->{$snp}->{$taxon});
	    push (@nodeouttaxa, $taxon);
	}
	
	# populate the ingroup
	elsif (exists ($intax->{$taxon})){
	    push (@nodein, $alndata->{$snp}->{$taxon});
            push (@nodeintaxa, $taxon);
	}
	
	# die if taxon is unassigned
	else{
	    print STDERR "Unknown taxon\n";
	    die;
	}
    }
    my $nodeoutstring = join ",", sort (@nodeouttaxa);
    my $nodeinstring = join ",", sort (@nodeintaxa);

    # check to see if the ingroup is internally consistent
    if (keys %{{map {$_, 1} @nodein}} == 1){
	my $basein = $nodein[0];
	    
	# check to see if the outgroup bases are all diff
	# from the ingroup -- don't have to be internally consistent
	my $compres = 0;
	my $baseouts = {};
	foreach my $baseout (@nodeout){
	    if ($baseout eq $basein){
		$compres++;
	    }
	    else{
		$baseouts->{$baseout}++;
		next;
	    }
	}
	
	# print if the conditions are satisfied
	if ($compres > 0){
	    next;
	}
	else {
	    my $baseoutstr = join (",", sort (keys %$baseouts));
	    print R "$snp\t$basein\t$baseoutstr\n";
	}
    }
    else {
	next;
    }
}
close (R);
close (A);
