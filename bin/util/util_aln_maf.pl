#!/usr/bin/perl -w
use strict;
use Getopt::Long;

my ($infile, $reference);
GetOptions(
	   'i|infile=s'       => \$infile,
	   'r|ref=s'   => \$reference, # SeRP62A.SeRP62A
	   );

# stream the maf file into memory
my $aln = {};
my $header;
my $counter = 0;
open (F, "$infile");
while (my $line = <F>){
    chomp $line;
    if ($line =~m/^\#\#/){
	$header = $line;
    }
    elsif ($line =~m/^a/){
	$counter++;
        $aln->{$counter}->{'data'} = $line;
    }
    elsif ($line =~m/^s/){
	my @splits = split (/\s+/, $line, 3);
	$aln->{$counter}->{'aln'}->{$splits[1]} = $splits[2];
    }
    else {
	next;
    }
}
close (F);

# reorder each alignment block with respect to
# the reference when available
print "$header\n";
foreach my $block (sort {$a <=> $b} keys %$aln){
    
    # print the alignment block data line
    print "$aln->{$block}->{'data'}\n";

    # print the reference alignment first if it exists
    if ($aln->{$block}->{'aln'}->{$reference}){
	print "s $reference $aln->{$block}->{'aln'}->{$reference}\n";
	
	# print all the rest
	foreach my $acc (sort keys %{$aln->{$block}->{'aln'}}){
	    if ($acc eq $reference){
		next;
	    }
	    else {
		print "s $acc $aln->{$block}->{'aln'}->{$acc}\n";
	    }
	}
    }
    # print in arbitrary order if reference doesn't exist
    else {
	foreach my $acc (sort keys %{$aln->{$block}->{'aln'}}){
	    print "s $acc $aln->{$block}->{'aln'}->{$acc}\n";
	}
    }
}

















=head
my $alnin = Bio::AlignIO->new(-file   => "$infile",
			      -format => 'maf');

# get aln data
my $alncounter = 0;
while (my $aln = $alignobj->next_aln()){
    my $label
    my $score = $aln->{'score'};
    $alncounter++;

    foreach my $seq ($alnobj->each_seq){
	$seqcounter++;
	my $id        = $seq->display_id;
	my $sequence = $seq->seq;
	

}
=cut
