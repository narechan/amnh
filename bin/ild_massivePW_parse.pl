#!/usr/bin/perl -w

=head1 NAME

ild_massivePW_parse.pl

=head1 SYNOPSIS

ild_massivePW_parse.pl

Options:

--outdir       is your output dir for the run

the matrix must be in
nexus format.

Requires the bioperl libs.
Requires Lild.pm

=head1 DESCRIPTION

This program does the parse on ild_massivePW.pl

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

my ($help, $outdir);
GetOptions(
    'h|help'           => \$help,
    'o|outdir=s'       => \$outdir,
    ) or pod2usage;
pod2usage if $help;

for my $option ($outdir){
    (warn ("Missing a required option\n") and pod2usage)
        unless ($option);
}


#####MAIN#####

# instantiate the object and load required data
my $ildobj = Ild->new;

# parse all the sg commands
my $sglens = {};
my $sgnumtrees = {};

opendir (L, "$outdir/sg_logs");
my @sglogs = sort (readdir (L));
shift @sglogs;
shift @sglogs;
closedir (L);

foreach my $sglog (@sglogs){
    my ($lens, $numtrees) = 
	$ildobj->parse_treelens ("$outdir/sg_logs/$sglog");

    foreach my $part (keys %$lens){
	$sglens->{$part} = $lens->{$part};
	$sgnumtrees->{$part} = $numtrees->{$part};
    }
}

# parse all the pw commands
my $pwlens = {};
my $pwnumtrees = {};

opendir (X, "$outdir/pw_logs");
my @pwlogs = sort (readdir (X));
shift @pwlogs;
shift @pwlogs;
closedir (X);

foreach my $pwlog (@pwlogs){
    my ($lens, $numtrees) =
        $ildobj->parse_treelens ("$outdir/pw_logs/$pwlog");

    foreach my $part (keys %$lens){
        $pwlens->{$part} = $lens->{$part};
        $pwnumtrees->{$part} = $numtrees->{$part};
    }
}
da
# calc the ILD stats
foreach my $pwlen (keys %$pwlens){
    my ($part1, $part2) = split (/v/, $pwlen);
    
    my $ild = 
	($pwlens->{$pwlen} - ($sglens->{$part1} + $sglens->{$part2})) / $pwlens->{$pwlen};
     print "$pwlen\t$pwlens->{$pwlen}\t$sglens->{$part1}\t$sglens->{$part2}\t$ild\n";
    
}
