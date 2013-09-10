#!/usr/bin/perl -w

=head1 NAME

ild_parse.pl

=head1 SYNOPSIS

ild_parse.pl

Options:

--outdir is your output dir
--config is your configfile
--matrix is your alignfile
--window       is your window size if doing                                                                
                LILD sildeRule (optional)                                                                 
--motion       is the number of residues the                                                            
                window moves per calculation (5' -> 3')                                                   
                if doing LILD slideRule (optional)                                                     
--start        is the start position of slideRule; set to 1 by default (optional)                      
--end          is the end position of slideRule; set to nchar by default    (optional)                     

The matrix must be in nexus format.
Requires the bioperl libs. 
Requires Ild.pm

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

#####SETUP#####

use strict;
use Getopt::Long;
use Pod::Usage;
use Ild;

my ($help, $matrixfile, $configfile, $outdir, $window, $motion, $start, $end);
GetOptions(
    'h|help'           => \$help,
    'm|matrix=s'       => \$matrixfile,
    'c|config=s'       => \$configfile,
    'o|outdir=s'       => \$outdir,
    'w|window=s'       => \$window,
    'j|motion=s'       => \$motion,
    's|start=s'        => \$start,
    'e|end=s'          => \$end,
    ) or pod2usage;
pod2usage if $help;

for my $option ($matrixfile, $configfile, $outdir){
    (warn ("Missing a required option\n") and pod2usage)
        unless ($option);
}

#####MAIN#####

# instantiate the object and load stuff we need
my $ildobj = Ild->new;
$ildobj->load_config ($configfile);
$ildobj->load_aln    ($matrixfile);

# define start as position 1 and and end as nchar                                                        
# unless otherwise specified                                                                          
my $nchar = $ildobj->get_nchar;
($start = 1) unless ($start);
($end   = $nchar) unless ($end);

# get the charsets depending on expt (slide rule or defined charsets)                               
my $charsets;
if ( ($window) and ($motion) and ($start) and ($end) ){
    $charsets = $ildobj->generate_partitions ($window, $motion, $start, $end);
}
else {
    $charsets = $ildobj->get_charsets;
}

# cycle through the output files, parse, and
# store the ild
opendir (R, "$outdir/logs");
my @logs = sort (readdir(R));
shift @logs;
shift @logs;
closedir (R);

open (RES, ">$outdir/summary");
foreach my $log (@logs){
    my ($pval, $slen, $ildchars) = $ildobj->parse_ild ("$outdir/logs/$log");
#    my ($partition, $stuff) = split (/\./, $log, 2);
    if ( ($window) and ($motion) and ($start) and ($end) ){
	$ildchars =~s/C//;
#	print RES "$partition\t$charsets->{$partition}\t$pval\t$slen\n";
	print RES "$ildchars\t$charsets->{$ildchars}\t$pval\t$slen\n";
    }
    else {
#	print RES "$partition\t$pval\t$slen\n";
	print RES "$ildchars\t$pval\t$slen\n";
    }
}
close (RES);
