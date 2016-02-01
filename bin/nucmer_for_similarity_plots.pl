#!/usr/bin/perl -w

=head1 NAME

nucmer_for_similarity_plots.pl

=head1 SYNOPSIS

Options:

 --help        Show brief help and exit
 --fasta       Is your query fasta
 --ref         Is your ref fasta
 --config      Is the configuration for your blast
 --outdir      Is your output dir
 --length      Is the length of your ref genome

mummer executables must be in your path

The config must also specify parameters for mummer executables

NUCMER=...
SHOW-COORDS=...

=head1 DESCRIPTION

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
use Statistics::Descriptive;

my ($help, $fasta, $ref, $config, $outdir, $length);
GetOptions(
    'h|help'          => \$help,
    'f|fasta=s'       => \$fasta,
    'o|outdir=s'      => \$outdir,
    'r|ref=s'         => \$ref,
    'c|config=s'      => \$config,
    'l|length=s'       => \$length,
    ) or pod2usage;

pod2usage if $help;

for my $option ($fasta, $outdir, $ref, $config, $length){
    (warn ("Missing a required option\n") and pod2usage)
	unless ($option);
}

`mkdir $outdir`;

#####MAIN#####

# get the config
my $conf = parse_config ($config);

# do the mummer algns
my $fasta_name;
if ($fasta =~/\//g){
    $fasta =~m/.*\/(.*)$/;
    $fasta_name = $1;
}

else {
    $fasta_name = $fasta;
}

my $ref_name;
if ($ref =~/\//g){
    $ref =~m/.*\/(.*)$/;
    $ref_name = $1;
}

else {
    $ref_name = $ref;
}

mummer ($conf, "$fasta", "$ref", "$outdir/$fasta_name.$ref_name"); 

# parse the show-coords file for pids
my $refdata = {};
open (F, "$outdir/$fasta_name.$ref_name.showcoords");
while (my $line = <F>){
    chomp $line;
    my @data  = split (/\t/, $line);
    my $start = $data[0];
    my $end   = $data[1];
    my $score = $data[6];
    
    for (my $i = $start; $i <= $end; $i++){
	push @{$refdata->{$i}}, $score;
    }
}
close (F);

# print out the per coordinate average pid with
# respect to the reference
for (my $j = 1; $j <= $length; $j++){
    if ($refdata->{$j}){
	my @data = @{$refdata->{$j}};
	my $statobj = Statistics::Descriptive::Full->new();
	$statobj->add_data(@data);
	my $mean  = $statobj->mean();
	print "$j\t$mean\n";
    }
    else{
	print "$j\t0\n";
    }
}

#####SUBS#####

sub mummer{
    my $conf   = shift;
    my $fasta  = shift;
    my $ref    = shift;
    my $out    = shift;

    # nucmer
    my $nucmer = "nucmer";
    ($nucmer .= " $conf->{'NUCMER'}") if ($conf->{'NUCMER'});
    $nucmer .= " -p $out";
    $nucmer .= " $ref";
    $nucmer .= " $fasta";
    `$nucmer`;
    
    # show-coords
    my $showcoords = "show-coords";
    ($showcoords .= " $conf->{'SHOW-COORDS'}") if ($conf->{'SHOW-COORDS'});
    $showcoords .= " $out.delta";
    $showcoords .= " > $out.showcoords";
    `$showcoords`;

}

sub parse_config{
    my $file = shift;

    my %config;
    open (F, "$file");
    while (my $line = <F>){
        chomp $line;
        my ($key, $value) = split (/\=/, $line, 2);
        $config{$key} = $value;

    }

    return (\%config);
    close (F);
}
