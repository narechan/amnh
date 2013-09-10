#!/usr/bin/perl -w

# -f is the fasta file of nucs
# -o is the outdir
# -c is your config file

#####SETUP#####

use strict;
use Getopt::Long;
use Bio::AlignIO;
use Bio::Align::Utilities qw(aa_to_dna_aln);
use Bio::SeqIO;

my ($fasta, $outdir, $configfile);
GetOptions(
	   'f|fasta=s'  => \$fasta,
	   'o|outdir=s' => \$outdir,
	   'c|config=s' => \$configfile,
	   );

#####MAIN######

`mkdir -p $outdir/pts`;
`mkdir -p $outdir/alns`;

# parse the config                                                                  
my $conf = parse_config ($configfile);

# get the fasta name
my $fasta_name;
if ($fasta =~m/\//g){
    $fasta =~m/.*\/(.*)$/;
    $fasta_name = $1;
}
else {
    $fasta_name = $fasta;
}

# declare hash reference for storing nuc seq objs
my $nucobjs = {};

# cycle through the fasta file getting id and sequence
my $seqin   = Bio::SeqIO->new (-format => "fasta", 
			       -file   => $fasta);
open (T, ">$outdir/pts/$fasta_name");
while (my $seqobj = $seqin->next_seq()){
    my $id     = $seqobj->display_id();
    my $nucseq = $seqobj->seq();
    
    # populate the nuc objs
    $nucobjs->{$id} = $seqobj;
    
    # translate the sequence and check trans in each frame
    my $ptobj0 = $seqobj->translate(-frame=>0);
    my $ptobj1 = $seqobj->translate(-frame=>1);
    my $ptobj2 = $seqobj->translate(-frame=>2);
    my $bestpt;
    my $min = 10000000000000;
    foreach my $pobj ($ptobj0, $ptobj1, $ptobj2){
	my $pt = $pobj->seq;
	my $ct =()= $pt =~m/\*/g;
	
	if ($ct < $min){
	    $bestpt = $pt;
	    $min = $ct;
	}
	else {
	    next;
	}
    }

    # print out the best translation
    print T ">$id\n$bestpt\n";

    # find the complement                                                                                  
    my $nucseqcomp = $nucseq;
    $nucseqcomp =~ tr/ACGTacgt/TGCAtgca/;

    # reverese the sequence
    my $nucseqrev = reverse ($nucseq);
    
    # reverse complement the sequence
    my $nucseqrevcomp = $nucseqrev;
    $nucseqrevcomp =~ tr/ACGTacgt/TGCAtgca/;
}
close (T);

# aln the translation with mafft
my $alnname = mafft_run ($fasta_name, $conf, "$outdir/pts", "$outdir/alns");

# do the codon alignment
my $aaaln = Bio::AlignIO->new(-file   => "$outdir/alns/$alnname",
			      -format => "fasta");
my $aaobj = $aaaln->next_aln();
my $dnaaln = aa_to_dna_aln($aaobj, $nucobjs);
my $alnout = Bio::AlignIO->new(-file   => ">$outdir/alns/$alnname.codonaln",
			       -format => "fasta");
$alnout->write_aln($dnaaln);


##### SUBS #####

sub mafft_run{
    my $fas  = shift;
    my $conf = shift;
    my $in   = shift;
    my $out  = shift;

    # mafft                                                                       
    my $mafft = "mafft";
    ($mafft .= " $conf->{'MAFFT'}") if ($conf->{'MAFFT'});
    $mafft .= " $in/$fas";

    `$mafft 1>$out/$fas.mafft.aln 2>$out/$fas.mafft.stderr`;
    return ("$fas.mafft.aln");
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
