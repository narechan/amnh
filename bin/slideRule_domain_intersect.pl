#!/usr/bin/perl -w

=head1 NAME

slideRule_domain_intersect.pl

=head1 SYNOPSIS

  slideRule_domain_intersect.pl -- this program searches for overlap between protfam annotations
    and slideRule ILD annotations.

Options:

 --help        Show brief help and exit
 --srfile   The slide rule summary file
 --pffile     The protfam summary file
 --evalcutoff The eval cutoff for counting a domain

=head1 DESCRIPTION

Search for annotation overlap

Usage examp:

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
use Blast;
use Statistics::Descriptive;

my ($help, $srfile, $pffile, $evalcutoff);
GetOptions(
    'h|help'          => \$help,
    's|srfile=s'      => \$srfile,
    'p|pffile=s'      => \$pffile,
    'e|evalcutoff=s'  => \$evalcutoff,
    ) or pod2usage;

pod2usage if $help;

for my $option ($srfile, $pffile, $evalcutoff){
    (warn ("Missing a required option\n") and pod2usage)
	unless ($option);
}

#####MAIN#####

# parse the pffile
my @domainregions;
open (PF, "$pffile");
my $counter = 0;
my $name_cache = 0;
while (my $line = <PF>){
    chomp $line;
    $counter++;
    my @pfdata = split (/\t/, $line);
    $pfdata[0] =~s/\s/_/g;
    my $seqname = $pfdata[0];
    my $domname = $pfdata[5];
    my $start   = $pfdata[6];
    my $end     = $pfdata[7];
    my $evalue  = $pfdata[8];
    $evalue =~s/E/e/;
    
    next if ($evalue > $evalcutoff);

    if ($seqname ne $name_cache){
	$name_cache = $seqname;
	$counter = 1;
    }

    my $acc = $seqname . $counter . "_" . $domname;
    
    push @domainregions,
    {
#	'chrom' => $acc,
	'chrom' => $seqname,
	'start' => $start,
	'end'   => $end,
    };
}
close (PF);

# get the uniq domain regions
my $blastobj = Blast->new;
my $coverage = $blastobj->genome_coverage (\@domainregions);

# parse the sr ild file and determine if you're 
# in or out of a domain region
my @in;
my @out;
open (SR, "$srfile");
while (my $line = <SR>){
    chomp $line;
    my ($part, $range, $pval, $score) = split (/\t/, $line);
    my ($start, $end) = split (/\-/, $range);
    my $midpt = ($start + $end) / 2;

    foreach my $domain (keys %$coverage){
	next unless ($srfile =~m/$domain/);

	my $cachein = 0;
	foreach my $domainregion (sort {$a <=> $b} keys %{$coverage->{$domain}}){
            my $start = $coverage->{$domain}->{$domainregion}[0];
            my $end   = $coverage->{$domain}->{$domainregion}[1];
            my $len   = $end - $start + 1;

	    if (($midpt >= $start) and ($midpt <= $end)){
		$cachein++;
	    }
        }
	
	if ($cachein > 0){
	    push (@in, $pval);
	}
	else {
	    push (@out, $pval);
	}
    }

}
close (SR);

# get the stats for all windows in a domain and 
# all windows out of domains
my $statin = Statistics::Descriptive::Full->new();
$statin->add_data(@in);
my $inmean = $statin->mean();
my $incount = $statin->count();

my $statout = Statistics::Descriptive::Full->new();
$statout->add_data(@out);
my $outmean = $statout->mean();
my $outcount = $statout->count();
    
print "IN\t$evalcutoff\t$incount\t$inmean\nOUT\t$evalcutoff\t$outcount\t$outmean\n";
