#!/usr/bin/perl -w

=head1 NAME

blast_parse_self.pl

=head1 SYNOPSIS

  blast_parse_self.pl -- 
              

Options:

 --help        Show brief help and exit
 --blastreport Is your raw blast report from legacy blast of blast against self
 --length      Is the min length of the hsps considered (all hsps since a single blast to self was done)
 --identity    Is the percent identity above which regions are screened   
 --outdir      Is your output dir

=head1 DESCRIPTION

Given a blast to self result, output regions or signifncant hits that are not to self

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

#####
#####     

#####SETUP#####

use strict;
use Getopt::Long;
use Pod::Usage;
use Parallel::ForkManager;
use Bio::SearchIO;

my ($help, $blastreport, $ident, $length);
GetOptions(
    'h|help'          => \$help,
    'b|blastreport=s' => \$blastreport,
    'i|identity=s'    => \$ident,
    'l|length=s'      => \$length,
    ) or pod2usage;

pod2usage if $help;

for my $option ($blastreport, $ident, $length){
    (warn ("Missing a required option\n") and pod2usage)
	unless ($option);
}

#####MAIN#####

my $in = Bio::SearchIO->new(-file   => "$blastreport",
			    -format => 'blast');
my $result = $in->next_result;
my $hit_count = 0;
while( my $hit = $result->next_hit ) {
    $hit_count++;
	
    my $hsp_count = 0;
    while (my $hsp = $hit->next_hsp ) {
	$hsp_count++;
	my $frac   = $hsp->frac_identical;
	my $hstart = $hsp->start('hit');
	my $hend   = $hsp->end('hit');
	my $qstart = $hsp->start('query');
	my $qend   = $hsp->end('query');
	my $qlen   = $hsp->length('query');
	
	if (($hstart == $qstart) and ($hend == $qend)){
	    print STDERR "$hsp_count\tSELF\t$qstart\t$qend\t$hstart\t$hend\t$qlen\t$frac\n";
	    next;
	}
	elsif ($qlen < $length){
	    print STDERR "$hsp_count\tSHORT\t$qstart\t$qend\t$hstart\t$hend\t$qlen\t$frac\n";
	    next;
	}
	elsif ($frac < $ident){
	    print STDERR "$hsp_count\tDIV\t$qstart\t$qend\t$hstart\t$hend\t$qlen\t$frac\n";
	    next;
	}
	else{
	    print STDERR "$hsp_count\tYES\t$qstart\t$qend\t$hstart\t$hend\t$qlen\t$frac\n";
	    print "$qstart-$qend\n";
	}
    }
}
