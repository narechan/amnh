#!/usr/bin/perl -w

=head1 NAME

blast_postprocess.pl

=head1 SYNOPSIS

  blast_postprocess.pl -- this program post-processes output from blast_forked.pl

Options:

 --help        Show brief help and exit
 --in          Is your input file

=head1 DESCRIPTION

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
use Bio::SearchIO;

my ($help, $infile);
GetOptions(
    'h|help'          => \$help,
    'i|infile=s'      => \$infile,
    ) or pod2usage;

pod2usage if $help;

for my $option ($infile){
    (warn ("Missing a required option\n") and pod2usage)
	unless ($option);
}

#####MAIN#####

open (F, "$infile");
my $hitornot = {};
my @queries;
while (my $line = <F>){
    chomp $line;
    
    my ($qname,
	$qlen,
	$hname,
	$hlen,
	$hcount,
	$hsprank,
	$hspeval,
	$hspscore,
	$hspfrac,
	$hspstartq,
	$hspendq,
	$hspgapsq,
	$hspfracq,
	$hspstrandq,
	$hspstarth,
        $hspendh,
        $hspgapsh,
        $hspfrach,
        $hspstrandh,
	$descripq,
	$descriph) = split (/\t/, $line);
    
    (push (@queries, $qname)) unless (exists ($hitornot->{$qname}));
    
    if ($hname eq "No hits found"){
	$hitornot->{$qname} = 0;
    }
    else {
	$hitornot->{$qname}++;
    }
}
close (F);

foreach my $query (@queries){
    print "$query\t$hitornot->{$query}\n";
}
