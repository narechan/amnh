#!/usr/bin/perl -w

=head1 NAME

ddh_postprocess.pl

=head1 SYNOPSIS

ddh_postprocess.pl

Options:

--indir is the dir that contains your results
--outdir is your output dir
--cutoff is the slope cutoff

Requires the bioperl libs. 

=head1 DESCRIPTION

processes the data from the ddh system

=head1 AUTHOR

Apurva Narechania
anarechania *a|t* amnh.org

=head1 COPYRIGHT

Copyright (c) 2008 American Museum of Natural History

This library is free software;  
you can redistribute it and/or modify 
it under the same terms as Perl itself.

=cut

# ----------------------------------------------------

#####SETUP#####

use strict;
use Getopt::Long;
use Pod::Usage;

my ($help, $indir, $outdir, $cutoff);
GetOptions(
    'h|help'          => \$help,
    'i|indir=s'       => \$indir,
    'c|cutoff=s'      => \$cutoff,
    'o|outdir=s'      => \$outdir,
    ) or pod2usage;
pod2usage if $help;

for my $option ($indir, $outdir, $cutoff){
    (warn ("Missing a required option\n") and pod2usage)
        unless ($option);
}

`mkdir -p $outdir`;

#####MAIN#####

# get the results
opendir (R, "$indir");
my @expts = sort (readdir(R));
shift @expts;
shift @expts;
closedir (R);

# parse and combine the summary files
open (SUM, ">$outdir/summary.slopes");
open (PTS, ">$outdir/summary.points");
foreach my $expt (sort @expts){

    my @expt = split (/\-v\-/, $expt);
    my $qry = $expt[0];
    my $ref = $expt[1];
    
    my @ref = split (/\-/, $ref);
    my $tech = pop @ref;
    $ref = join "-", @ref;
    
    # store the expt's data
    my @data;
    open (S, "$indir/$expt/summary");
    while (my $line = <S>){
	chomp $line;
	
	my @array = split (/\t/, $line);
	push (@data, [@array]);
    }
    close (S);

    # sort the array by reads (coverage)
    my @sorted = sort {$a->[0] <=> $b->[0]} @data;
    
    # calculate the slopes
    # and take out first data points under the cutoff
    my $counter = 0;
    my %slopes;
    my $inflect;
    foreach my $point (@sorted){
	my $currentcov = @$point[2];
	my $currentref = @$point[7];
	
	if ($sorted[$counter + 1]){
	    my $nextcov    = $sorted[$counter + 1][2];
	    my $nextref    = $sorted[$counter + 1][7];
	
	    my $slope = ($nextref - $currentref) / ($nextcov - $currentcov);
	    
	    if (($slope < $cutoff) and (! defined($inflect))){
		print PTS "$qry\t$ref\t$tech\t$sorted[$counter + 1][4]\t$nextcov\t$nextref\t";
		$inflect = $nextcov;
	    }
	    
	    $slopes{$currentcov} = $slope;
	    $counter++;
	}
	
	else {
	    print PTS "$currentcov\t$currentref\n";
	    last;
	}
    }
    
    # print out the SUM data
    print SUM "$qry\t$ref\t$tech\n";
    foreach my $point (@sorted){
	my $joined = join "\t", @$point;
	print SUM "$joined";
	
	if (exists ($slopes{@$point[2]})){
	    print SUM "\t$slopes{@$point[2]}\n";
	}
	else {
	    print SUM "\n";
	}
    }

}
close (SUM);
close (PTS);
