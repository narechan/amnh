#!/usr/bin/perl -w

# -i is the alignment

#####SETUP#####

use strict;
use Getopt::Long;
use Bio::AlignIO;
use Algorithm::Numerical::Shuffle qw /shuffle/;

my ($alnfile);
GetOptions(
	   'i|alnfile=s'       => \$alnfile,
	   );

#####MAIN#####
my $alnin = Bio::AlignIO->new(-file   => "$alnfile",
			      -format => "phylip");
my $alnobj = $alnin->next_aln();
my $alnlen = $alnobj->length;

my $seqcounter = 0;

my $alndata = {};
foreach my $seq ($alnobj->each_seq){
    $seqcounter++;
    my $id        = $seq->display_id;
    my $sequence = $seq->seq;
    my @sequence = split "", $sequence;
#    $alndata->{$id} = $sequence;
    $alndata->{$id} = [@sequence]; 
}

my @columns;
my @tid;
for (my $j = 0; $j <= $alnlen - 1; $j++){
    my @column;
    foreach my $tid (sort keys %$alndata){
	push (@column, $alndata->{$tid}[$j]);
	(push (@tid, $tid)) if ($j == $alnlen - 1);
    }
    my $column = join "", @column;
#    push @columns, [@column];
    push (@columns, $column);
}

my @columns_shuff = shuffle (@columns);

my $tcount = -1;
foreach my $t (@tid){
    print STDERR "$t\n";
    print "$t\t";
    $tcount++;
    foreach my $col (@columns_shuff){
	my @col = split "", $col;
	print "$col[$tcount]";
    }
    print "\n";
}
