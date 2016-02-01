#!/usr/bin/perl -w

#####SETUP#####
use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;

my ($help, $posfile, $chipfile, $outdir);
GetOptions(
    'h|help'          => \$help,
    'p|posfile=s'   => \$posfile, #position file you want to interrogate
    'c|chipfile=s'   => \$chipfile, #positions of chipseq peaks and intervals
    'o|outdir=s'    => \$outdir,
	   ) or pod2usage;
pod2usage if $help;

`mkdir -p $outdir`;

#####MAIN#####

# parse the chip seq data
my $chip = {};
open (C, "$chipfile");
while (my $line = <C>){
    chomp $line;
    my ($regu, $gene, $pval, $score, $vpm, $start, $end, $center, $seq) = split (/\t/, $line);
    my $coords = $start . "-" . $end;
#    $chip->{$regu} = $coords;
    $chip->{$coords} = $regu;
}
close (C);

my $tfs = {};
open (X1, "$posfile");
open (A, ">$outdir/posfile.annotated");
while (my $line1 = <X1>){
    chomp $line1;
    next if ($line1 =~m/^coord/);
    
    # note the for strands, 1 is F/+ and 2 is R/-
#    my ($coord, $position, $seq, $strand, $cov, $qv, $ipdr) = split (/\t/, $line1);
    my ($position, $ipdr) = split (/\t/, $line1);

    # bail if the site is not mapped to rv
    if (($position eq "N") or ($position eq "X")){
	print A "$line1\tNO_MAP\n";
	next;
    }

    # brute force binning                      
    my @regs;
    foreach my $coordset (sort keys %$chip){
        my ($start, $end) = split (/\-/, $coordset);
        if (($start <= $position) and ($position <= $end)){
	    push (@regs, $chip->{$coordset});
	    push @{$tfs->{$chip->{$coordset}}}, $ipdr;
        }
        else {
            next;
        }
    }

    my $regnumber = @regs;
    my $regstring = join ";", @regs;
    print A "$line1\t$regstring\t$regnumber\n";
}
close (X1);

open (B, ">$outdir/tfs.ipdrs");
foreach my $tf (keys %$tfs){
    my $ipdrs = join "\t", @{$tfs->{$tf}};
    print B "$tf\t$ipdrs\n";
}
close (B);

