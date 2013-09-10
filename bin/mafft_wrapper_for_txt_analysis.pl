#!/usr/bin/perl -w

=head1 NAME

mafft_wrapper_for_text_analysis.pl

=head1 SYNOPSIS

  mafft_wrapper_for_text_analysis.pl -- this program will partition a clean text
    string  into "gene" sizes of arbitrary length,
    randomize the gene's chars and create gene members for additional taxa,
    align genes in order as extracted from text fastas,
    and recombine to give concatenated alignment and partitions

Options:

 --help        Show brief help and exit
 --infile      Contains the text to be analyzed (single string slurp)
 --config      Is the configuration for your run (MAFFT params)
 --chars       The number of chars in each "gene"
 --total       The total number of chars you want
 --members     The number of taxa members to simulate randomly
 --outdir      Is your output dir
 --procs       Is the number of procs you want to use

mafft must be in your path

All text should be on one line in one string (no special chars)
The config must specify parameters for the mafft run

=head1 DESCRIPTION

Run mafft on extracted "genes" from a "genomic" text, and random taxa

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
use Bio::SeqIO;
use Bio::AlignIO;
use Algorithm::Numerical::Shuffle qw /shuffle/;
use Parallel::ForkManager;

my ($help, $infile, $config, $outdir, $chars, $total, $members, $procs);
GetOptions(
    'h|help'          => \$help,
    'i|infile=s'      => \$infile,
    'o|outdir=s'      => \$outdir,
    'c|config=s'      => \$config,
    'x|chars=s'       => \$chars,
    't|total=s'       => \$total,
    'm|members=s'     => \$members,
    'p|procs=s'       => \$procs,
    ) or pod2usage;

pod2usage if $help;

for my $option ($infile, $outdir, $config, $chars, $total, $members, $procs){
    (warn ("Missing a required option\n") and pod2usage)
	unless ($option);
}

`mkdir -p $outdir/alns`;
`mkdir -p $outdir/subfastas`;
`mkdir -p $outdir/stderr`;

#####MAIN#####

# get the config
my $conf = parse_config ($config);

# slurp the text in the file
open (F, "$infile");
my $txt = do {local $/; <F>};
chomp $txt;
close (F);

# generate the subfasta files for each "gene"
print STDERR "Build subfastas\n";
my $counter = 0;
my $totlen  = 0;
while (1){
#    last;
    $counter++;
    
    open (O, ">$outdir/subfastas/part$counter.fasta");
    my $substring = substr ($txt, 0, $chars, '');
    my $length = length $substring;
    $totlen += $length;
    
    print O ">part$counter.tax1\n$substring\n";
    
    my @substring = split (//, $substring);
    for (my $i = 2; $i <= $members; $i++){
	my @shuffsubstr = shuffle @substring;
	my $shuffsubstr = join "", @shuffsubstr;
	print O ">part$counter.tax$i\n$shuffsubstr\n";
    }
    close (O);
    
    last if ($totlen >= $total);
}

# do the mafft alignments
opendir (A, "$outdir/subfastas");
my @fastas = sort readdir (A);
shift @fastas;
shift @fastas;
closedir (A);
my $pm = Parallel::ForkManager->new($procs);
foreach my $fasta (sort @fastas){
#    last;
    $pm->start and next;
    
    my $seqin = Bio::SeqIO->new (-format=>'Fasta', -file=>"$outdir/subfastas/$fasta");
    my $seqcount = 0;
    while (my $sequence_obj = $seqin->next_seq()){
	$seqcount++;
    }
    if ($seqcount > 1){
	print STDERR "Mafft $fasta\n";
	mafft_run ($fasta, $conf, $outdir);
    }
    else {
	print STDERR "Mafft skipped $fasta\n";
    }
    $pm->finish;
}
$pm->wait_all_children;

# concatenate the mafft alignments and 
# define charpars
my $alndata = {};
my $alnlens = {};
opendir (J, "$outdir/alns");
my @alns = sort readdir (J);
shift @alns;
shift @alns;
closedir (J);

my $alnlentot = 0;
foreach my $aln (@alns){
    my @alnname = split (/\./, $aln);
    my $alnin = Bio::AlignIO->new(-file   => "$outdir/alns/$aln",
				  -format => "fasta");
    
    my $alnobj = $alnin->next_aln();
    my $alnlen = $alnobj->length;
    $alnlentot += $alnlen;
    $alnlens->{$alnname[0]} = $alnlen;
    
#    my %acc_cache = %accessions;
    foreach my $seq ($alnobj->each_seq){
	my $id        = $seq->display_id;
	my ($fasta_acc, $aln_acc) = split (/\./, $id);
#	delete $acc_cache{$acc};
	my $sequence = $seq->seq;
	$alndata->{$alnname[0]}->{$aln_acc} = $sequence;
    }
    
    # populate empty taxa with missing chars
#    foreach my $acc (keys %acc_cache){
#	my $sequence = '?' x $alnlen;
#	$alndata->{$alnname[0]}->{$acc} = $sequence;
#    }

}
    
# sort print the concatenation and charpars
open (CAT, ">$outdir/matrix");
open (PRT, ">$outdir/partitions");
print CAT "#NEXUS\n";
print CAT "BEGIN DATA;\n";
print CAT "DIMENSIONS NTAX=$members NCHAR=$alnlentot;\n";
print CAT "FORMAT INTERLEAVE SYMBOLS=\"ABCDEFGHIKLMNPQRSTUVWXYZ\" DATATYPE=PROTEIN MISSING=? GAP=-;\n";
print CAT "MATRIX\n";
print PRT "BEGIN SETS;\n";

my $start = 1;
my $end;
foreach my $count (sort keys %$alndata){
    $end = $start - 1 + $alnlens->{$count};
    print CAT "[Partition $count length $alnlens->{$count} chars $start-$end]\n";
    print PRT "CHARSET $count=$start-$end;\n";
    foreach my $sp (sort keys %{ $alndata->{$count} }){
       print CAT "$sp\t$alndata->{$count}->{$sp}\n";                                                       
    }
    print CAT "\n";
    $start = $end + 1;
}

print CAT ";\n";
print CAT "END;\n\n";
print PRT "END;\n";

close (CAT);
close (PRT);

`cat $outdir/matrix $outdir/partitions > $outdir/nexus`;

#####SUBS#####

sub mafft_run{
    my $fas  = shift;
    my $conf = shift;
    my $out  = shift;

    # mafft
    my $mafft = "mafft";
    ($mafft .= " $conf->{'MAFFT'}") if ($conf->{'MAFFT'});
    $mafft .= " $out/subfastas/$fas";

    `$mafft 1>$out/alns/$fas.aln 2>$out/stderr/$fas.stderr`;
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
    close (F);

    return (\%config);
}

