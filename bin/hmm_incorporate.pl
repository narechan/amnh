#!/usr/bin/perl

# This program will score a fasta file against a library of 
# hmms (oid based, usually), assign sequences to families given 
# a minimum cutoff,  rebuild the alignment (mafft), and 
# reconstruct an ML best tree (raxml) and bootstrap tree (raxml).

use Getopt::Std;
use Bio::Seq;
use Bio::SeqIO;
use Bio::AlignIO;
use Bio::Index::Fasta;

my %opts = ();
getopts ('h:f:c:o:d:b:', \%opts);
my $fasta  = $opts{'f'}; # input fasta file to score
my $oiddir = $opts{'d'}; # directory of oid families as fasta files (can be alignments or sequences)
my $hmmdb = $opts{'h'}; # the hmm database to score against
my $cutoff = $opts{'b'}; # using bit score cutoff here, ie 50
my $outdir = $opts{'o'}; # output
my $config = $opts{'c'}; # contains raxml parameters for best (RAXML_BEST) and boots (RAXML_BOOTS;RAXML_CONBOOTS)
                         # optionally contains parameters for MAFFT and HMMSCAN
`mkdir -p $outdir/alns`;
`mkdir -p $outdir/trees`;

# parse the config
my $conf = {};
open (F, "$config");
while (my $line = <F>){
    chomp $line;
    my ($key, $value) = split (/\=/, $line, 2);
    $conf->{$key} = $value;

}
close (F);

# get the fasta name
my $fasta_name;
if ($fasta =~/\//g){
    $fasta =~m/.*\/(.*)$/;
    $fasta_name = $1;
}
else {
    $fasta_name = $fasta;
}

# index the fasta file
my $index = Bio::Index::Fasta->new(-filename => $fasta . ".idx", -write_flag => 1);
$index->make_index($fasta);

# score using hmmscan
print STDERR "Scanning HMMs\n";
my $hmmcmd = "hmmscan";
$hmmcmd   .= " $conf->{'HMMSCAN'}";
$hmmcmd   .= " -o $outdir/$fasta_name.hmm.out --tblout $outdir/$fasta_name.hmm.tblout --notextw $hmmdb $fasta";
`$hmmcmd`;

# parse the hmmhit
my $fams = {};
my $scores = {};
open (A, "$outdir/$fasta_name.hmm.tblout");
while (my $line = <A>){
    chomp $line;
    next if ($line =~m/\#/);
    my @line = split (/\s+/, $line);
    
    if (exists ($fams->{$line[2]})){
        if ($line[5] > $scores->{$line[2]}){
            $fams->{$line[2]} = $line[0];
            $scores->{$line[2]} = $line[5];
        }
        else {
            next;
        }
    }
    else {
        $fams->{$line[2]} = $line[0];
        $scores->{$line[2]} = $line[5];
    }
}
close (A);

# add the new sequences to their respective oid families and log
# base of query fasta file to track which ones have no hits
open (L, ">$outdir/log");
my $newfams = {};
my $seqo = Bio::SeqIO->new (-format => 'Fasta', -file => "$fasta");
while (my $sobj = $seqo->next_seq()){
    my $seq = $sobj->display_id();
#foreach my $seq (keys %$fams){
    if ($fams->{$seq}){
	my $fam = $fams->{$seq};
	my $sequence = $index->fetch($seq);
	$newfams->{$fam}->{$seq} = $sequence->seq();
	print L "$seq\t$fam\t$scores->{$seq}\n";
    }
    else{
	print L "$seq\tNO_HITS\n";
    }
}

# add back the sequences already present in the families
foreach my $fam (keys %$newfams){
    my $seqin = Bio::SeqIO->new (-format => 'Fasta', -file => "$oiddir/$fam");
    while (my $sequence_obj = $seqin->next_seq()){
	my $id         = $sequence_obj->display_id();
	my $sequence   = $sequence_obj->seq();
	$sequence =~s/\-//g; #get rid of gaps if the families are aligned
	$newfams->{$fam}->{$id} = $sequence;
    }
}

# create new fasta files and realign
print STDERR "Building alignments\n";
foreach my $fam (keys %$newfams){
    open (F, ">$outdir/alns/$fam.fa");
    foreach my $acc (sort keys %{$newfams->{$fam}}){
	print F ">$acc\n$newfams->{$fam}->{$acc}\n";
    }
    close (F);

    my $mafftcmd = "mafft";
    $mafftcmd   .= " $conf->{'MAFFT'}";
    $mafftcmd   .= " $outdir/alns/$fam.fa > $outdir/alns/$fam.aln";
    `$mafftcmd`;

    my $alnin = Bio::AlignIO->new(-file   => "$outdir/alns/$fam.aln",
				  -format => 'fasta');
    my $alnout = Bio::AlignIO->new(-file   => ">$outdir/alns/$fam.phy",
				   -idlength => 40,
				   -format => 'phylip');
    while (my $aln = $alnin->next_aln){
	$alnout->write_aln($aln);
    }
}

# run the tree searches
print STDERR "Running trees\n";
foreach my $fam (keys %$newfams){
    my $raxmlcmd = "raxml728";
    $raxmlcmd   .= " $conf->{'RAXML_BEST'}";
    $raxmlcmd   .= " -n $fam -s $outdir/alns/$fam.phy";
    `$raxmlcmd`;
    
    my $raxmlbootcmd = "raxml728";
    $raxmlbootcmd   .= " $conf->{'RAXML_BOOTS'}";
    $raxmlbootcmd   .= " -n $fam.boots -s $outdir/alns/$fam.phy";
    `$raxmlbootcmd`;
    
    my $raxmlconbootcmd = "raxml728";
    $raxmlconbootcmd   .= " $conf->{'RAXML_CONBOOTS'}";
    $raxmlconbootcmd   .= " -s $outdir/alns/$fam.phy -t RAxML_bestTree.$fam -z RAxML_bootstrap.$fam.boots -n $fam.conboots";
    `$raxmlconbootcmd`;
    `mv RAxML_* $outdir/trees`;
}

