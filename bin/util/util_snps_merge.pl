#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use Pod::Usage;
use Bio::Seq;
use Bio::SeqIO;

my ($help, $infile1, $infile2, $snpqual);
GetOptions(
    'h|help'           => \$help,
    'a|infile1=s'        => \$infile1,
    'b|infile2=s'        => \$infile2,	   
    's|snpqual=s'       => \$snpqual,	   
    ) or pod2usage;

pod2usage if $help;

# parse the vcf
my $pileup = {};
open (PU, "$infile1");
while (my $line = <PU>){
    chomp $line;
    next if ($line =~m/\#/);
    
    my ($refseq, $refposy, $refbase, $snp, $cq, $sq, $rms, $cov, $bases, $quals, $rest) 
	= split (/\t/, $line);
    
    next if ($sq < $snpqual);
    
    my @pileup = split (/\t/, $line);
    shift @pileup;
    my $refpos = shift @pileup;
    
    $pileup->{$refpos} = [@pileup];
}
close (PU);

# parse the bioscope stuff
my $bioscope = {};
open (BS, "$infile2");
while (my $line = <BS>){
    chomp $line;
    
    my @bioscope = split (/\t/, $line);
    shift @bioscope;
    my $refpos = shift @bioscope;
    
    $bioscope->{$refpos} = [@bioscope];
}

foreach my $pos (sort {$a <=> $b} keys %$bioscope){
    my $bstring = join ("\t", @{$bioscope->{$pos}});
    print "$pos\t$bstring";
    
    if (exists ($pileup->{$pos})){
	my $pstring = join ("\t", @{$pileup->{$pos}});
	print "\t$pstring\n";
	delete $pileup->{$pos};
    }
    else{
	print "\n";
    }
}

foreach my $pos (sort {$a <=> $b} keys %$pileup){
    my $pstring = join ("\t", @{$pileup->{$pos}});
    print "$pos\t\t\t\t\t\t\t\t$pstring\n";
}
