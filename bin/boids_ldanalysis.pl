#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use Statistics::Descriptive;
use Bio::Seq;
use Bio::SeqIO;

my ($help, $flockdir, $indir, $burnstart, $burnend, $ldfile); #where indir is the directory from flocking
GetOptions(
    'h|help'          => \$help,
    'f|flockdir=s'    => \$flockdir,
    'i|indir=s'       => \$indir,
    's|startburn=s'   => \$burnstart,
    'e|endburn=s'     => \$burnend,
    'l|ldfile=s'      => \$ldfile,
    ) or pod2usage;

pod2usage if $help;

for my $option ($indir, $flockdir, $burnstart, $burnend, $ldfile){
    (warn ("Missing a required option\n") and pod2usage)                            
        unless ($option);                                                           
}                                                                                   

####MAIN####                                                                         

# create a symbol library (each letter of the alphabet)
my $library = {};
my $c = -1;
my @letters;
foreach my $letter ("A".."Z"){
    $c++;
    $library->{$c} = $letter;
    push (@letters, $letter);
}
my $symbolstring = join " ", @letters;

# parse out the boids from fasta file names                                          
my $taxnum = {};
my $taxacount = {};
opendir (D, "$indir");
my @charsets = sort (readdir (D));
shift @charsets;
shift @charsets;
foreach my $charset (@charsets){
    $taxnum->{$charset} = 0;
    my $seqin = Bio::SeqIO->new (-format=>'Fasta', -file=>"$indir/$charset");
    while (my $sequence_obj = $seqin->next_seq()){
	my $id       = $sequence_obj->display_id();
        my $seq      = $sequence_obj->seq();
	my $len = length ($seq);
	my $matchstring = '\\' . '?' . "{$len}";
	($taxnum->{$charset}++) if ($seq =~m/$matchstring/);
	$taxacount->{$id}++;
    }
}	
closedir (D);

# parse and store the LD table
# if available store the precomputed ILD list                                         
my $ilds = {};
open (L, "$ldfile");
while (my $line = <L>){
    chomp $line;
    my ($pair, $pval) = split (/\t/, $line);
    my ($i1, $i2) = split (/v/, $pair);
    push @{$ilds->{$i1}}, $pval;
    push @{$ilds->{$i2}}, $pval;
}
close (L);

# parse the flocking info
opendir (F, "$flockdir");
my @repdirs = sort (readdir (F));
shift @repdirs;
shift @repdirs;
closedir (F);
my $repdirs = @repdirs;

my $pa = {};
foreach my $repdir (@repdirs){
    for (my $i = $burnstart; $i <= $burnend; $i++){
	print STDERR "$repdir\t$i\n";
	
	# if the optics file exists, parse it,
	# otherwise insert missing data for a replicate
	# that failed at this point
	if (-e "$flockdir/$repdir/flocks/$i.optics"){
	    open (I, "$flockdir/$repdir/flocks/$i.optics");
	    while (my $line = <I>){
		chomp $line;
		my ($flock, $gene) = split (/\t/, $line);
		if ($flock > 26){
		    print STDERR "Too many flocks: not enough characters!\n";
		    die;
		}
		else{
		    push (@{$pa->{$gene}}, $library->{$flock});
		}
	    }
	    close (I);
	}
	else {
	    foreach my $charset (@charsets){
		push (@{$pa->{$charset}}, "?");
	    }
	}
    }
}

my $freqs = {};
foreach my $tax (sort keys %$pa){
    my @paarray = @{$pa->{$tax}};
    
    my $a = 0;
    my $b = 0;
    foreach my $pa (@paarray){
	if ($pa eq "A"){
	    $a++;
	}
	else {
  	    $b++;
	}
    }
    my $freqa = $a / ($a + $b);
    $freqs->{$freqa}++;

    my @ilds = @{$ilds->{$tax}};
    my $zeros = 0;
    my $nonzeros = 0;
    foreach my $il (@ilds){
	if ($il == 0){
	    $zeros++;
	}
	else{
	    $nonzeros++;
	}
    } 
    my $zerofreq = $zeros / ($zeros + $nonzeros);

    my $statobj = Statistics::Descriptive::Full->new();
    $statobj->add_data(@ilds);

    my $mean  = $statobj->mean();
    my $median = $statobj->median();

    my @taxacount = keys %$taxacount;
    my $tcount = @taxacount;
    my $missingtaxa = $taxnum->{$tax};
    my $missingtaxafreq = $taxnum->{$tax} / $tcount;
    if (($freqa <= 0.2) and ($missingtaxafreq <= 0.2)){
	print "$tax\t$a\t$freqa\t$missingtaxa\t$missingtaxafreq\t$zeros\t$zerofreq\t$mean\t$median\n";
    }
}
=head
for (my $i = 0; $i <= 100; $i++){
    my $ith = $i / 100;
    if (exists ($freqs->{$ith})){
#	print STDERR "$ith\t$freqs->{$ith}\n";
    }
    else {
#	print STDERR "$ith\t0\n";
    }
}
=cut
