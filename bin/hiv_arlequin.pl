#!/usr/bin/perl -w

=head1 NAME

hiv_arelquin.pl

=head1 SYNOPSIS

hiv_arlequin.pl

Options:

--outdir is your output dir
--infile1 is your real data (fasta)
--infile2 is your random data (fasta)
--config is the arlecore settings file (must end in .ars)

Requires the bioperl libs. 

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
use Bio::SeqIO;
use Bio::Seq;
use Data::Compare;

my ($help, $infile1, $infile2, $outdir, $configfile);
GetOptions(
    'h|help'          => \$help,
    'a|infile1=s'     => \$infile1,
    'b|infile2=s'     => \$infile2, 	   
    'o|outdir=s'      => \$outdir,
    'c|config=s'      => \$configfile,
	   ) or pod2usage;
pod2usage if $help;

for my $option ($infile1, $infile2, $outdir, $configfile){
    (warn ("Missing a required option\n") and pod2usage)
        unless ($option);
}

`mkdir -p $outdir`;

#####MAIN#####

# initialize the hash key for aas
my %aakey = (
	     'A' => 1,
	     'R' => 2,
	     'N' => 3,
	     'D' => 4,
	     'C' => 5,
	     'E' => 6,
	     'Q' => 7,
	     'G' => 8,
	     'H' => 9,
	     'I' => 10,
	     'L' => 11,
	     'K' => 12,
	     'M' => 13,
	     'F' => 14,
	     'P' => 15,
	     'S' => 16,
	     'T' => 17,
	     'W' => 18,
	     'Y' => 19,
	     'V' => 20,
	     'X' => 21,
	     );

# get file names
my $infile1name = get_name ($infile1);
my $infile2name = get_name ($infile2);

# read the pair of sequences and translate respecting gaps
my $alnlen = 0;
my $seqcount = 0;
my $realpts = {};
my $realseqin = Bio::SeqIO->new (-format=>'Fasta', -file=>"$infile1");
while (my $sequence_obj = $realseqin->next_seq()){
    $seqcount++;
    my $id       = $sequence_obj->display_id();
    my $seq      = $sequence_obj->seq();
    $alnlen      = $sequence_obj->length();
    
    my $prot_obj = $sequence_obj->translate;
    my $prot_seq = $prot_obj->seq();
    my @prot_seq = split (//, $prot_seq);
    
    my $counter = 0;
    foreach my $aa (@prot_seq){
	$counter++;
	$realpts->{$counter}->{$aa}++;
    }
}

my $randpts = {};
my $randseqin = Bio::SeqIO->new (-format=>'Fasta', -file=>"$infile2");
while (my $sequence_obj = $randseqin->next_seq()){
    my $id       = $sequence_obj->display_id();
    my $seq      = $sequence_obj->seq();

    my $prot_obj = $sequence_obj->translate;
    my $prot_seq = $prot_obj->seq();
    my @prot_seq = split (//, $prot_seq);

    my $counter = 0;
    foreach my $aa (@prot_seq){
	$counter++;
        $randpts->{$counter}->{$aa}++;
    }
}

# calculate the aa frequencies per position
# and write to arlequin input file
foreach my $pos (sort {$a <=> $b} keys (%$realpts)){
    my $compare = Compare ($realpts->{$pos}, $randpts->{$pos});
    
    if ($compare == 1){
	open (AR, ">$outdir/$infile1name.$pos.arp");
	print AR "Aborted: populations are the same!\n";
	close (AR);
	next;
    }
    else {
	open (AR, ">$outdir/$infile1name.$pos.arp");
	print AR "[Profile]\n\n";
	print AR "Title=\"$infile1name.$pos\"\n";
	print AR "NbSamples=2\n";
	print AR "DataType=FREQUENCY\n";
	print AR "GenotypicData=0\n";
	print AR "LocusSeparator=WHITESPACE\n";
	print AR "MissingData='?'\n";
	print AR "Frequency=REL\n";
	print AR "\n";
	print AR "[DATA]\n\n";
	print AR "[[SAMPLES]]\n";
	print AR "SampleName=\"$infile1name.$pos\"\n";
	print AR "SampleSize=$seqcount\n";
	print AR "SampleData={\n";
	foreach my $aa (sort {$aakey{$a} <=> $aakey{$b}} keys %aakey){
	    my $freq = 0;
	    if (exists ($realpts->{$pos}->{$aa})){
		$freq = $realpts->{$pos}->{$aa} / $seqcount;
		print AR "$aakey{$aa} $freq\n";
	    }
	    else {
		print AR "$aakey{$aa} $freq\n";
	    }
	}
	print AR "}\n";
	
	print AR "SampleName=\"$infile2name.$pos\"\n";
	print AR "SampleSize=$seqcount\n";
	print AR "SampleData={\n";
	foreach my $aa (sort {$aakey{$a} <=> $aakey{$b}} keys %aakey){
	    my $freq = 0;
	    if (exists ($randpts->{$pos}->{$aa})){
		$freq = $randpts->{$pos}->{$aa} / $seqcount;
		print AR "$aakey{$aa} $freq\n";
	    }
	    else {
		print AR "$aakey{$aa} $freq\n";
	    }
	}
	print AR "}\n\n";
	
	print AR "[[STRUCTURE]]\n";
	print AR "StructureName=\"$infile1name\"\n";
	print AR "NbGroups=1\n";
	print AR "Group={\n";
	print AR "\"$infile1name.$pos\"\n";
	print AR "\"$infile2name.$pos\"\n";
	print AR "}\n";
	close (AR);
	
	print "$infile1name.$pos\n";
	
	# run arlequin
#	my $cwd = `pwd`;
#	chomp $cwd;
#	`mkdir $outdir/$infile1name.$pos.arp.out`;
#	`cd $outdir/$infile1name.$pos.arp.out`;
#	`~apurva/packages/arlecore_linux/arlecore3512_64bit $outdir/$infile1name.$pos.arp $configfile`;
#	`cd $cwd`;
    }
}



####SUBS####

sub get_name{
    my $fasta = shift;
    my $fasta_name;
    if ($fasta =~/\//g){
        $fasta =~m/.*\/(.*)$/;
        $fasta_name = $1;
    }

    else {
        $fasta_name = $fasta;
    }
    return ($fasta_name);
}
   
    
