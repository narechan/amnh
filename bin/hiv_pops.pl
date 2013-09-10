#!/usr/bin/perl -w

=head1 NAME

hiv_pops.pl

=head1 SYNOPSIS

hiv_pops.pl

Options:

--infile is your input data file
--outdir is your output dir

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

my ($help, $infile, $outdir);
GetOptions(
    'h|help'          => \$help,
    'i|infile=s'      => \$infile,
    'o|outdir=s'      => \$outdir,
	   ) or pod2usage;
pod2usage if $help;

for my $option ($infile, $outdir){
    (warn ("Missing a required option\n") and pod2usage)
        unless ($option);
}

`mkdir -p $outdir`;

#####MAIN#####

# parse the infile
my $data = {};
my $seqs = {};
my @accs;
open (I, "$infile");
while (my $line = <I>){
    chomp $line;
    my ($acc,
	$a1,
	$a2,
	$b1,
	$b2,
	$c1,
	$c2,
	$seq) = split (/\t/, $line);

#    $data->{"a1"}->{$a1}->{$acc} = $seq;
#    $data->{"b1"}->{$b1}->{$acc} = $seq;
#    $data->{"c1"}->{$c1}->{$acc} = $seq;
#    $data->{"a2"}->{$a2}->{$acc} = $seq;
#    $data->{"b2"}->{$b2}->{$acc} = $seq;
#    $data->{"c2"}->{$c2}->{$acc} = $seq;

    ($data->{"a"}->{$a1}->{$acc} = $seq) if ($a1);
    ($data->{"b"}->{$b1}->{$acc} = $seq) if ($b1);
    ($data->{"c"}->{$c1}->{$acc} = $seq) if ($c1);
    ($data->{"a"}->{$a2}->{$acc} = $seq) if ($a2);
    ($data->{"b"}->{$b2}->{$acc} = $seq) if ($b2);
    ($data->{"c"}->{$c2}->{$acc} = $seq) if ($c2);

    $seqs->{$acc} = $seq;
    push (@accs, $acc);
}
close (I);

# write out the HLA population files and random samples
foreach my $hla_cat (sort keys %$data){
    `mkdir -p $outdir/$hla_cat`;
    
    foreach my $hla (sort keys %{$data->{$hla_cat}}){
	open (W, ">$outdir/$hla_cat/$hla.fasta");
	open (R, ">$outdir/$hla_cat/$hla-rand.fasta");

	# data structure to track which accs have been randomly sampled
	my $randelims = {};
	foreach my $acc (@accs){
	    $randelims->{$acc} = 1;
	}
	
	my $rands = {};
	foreach my $acc (sort keys %{$data->{$hla_cat}->{$hla}}){
	    print W ">$acc-$hla\n$data->{$hla_cat}->{$hla}->{$acc}\n";
	    
	    while (1) {
		my $randacc = $accs[rand @accs];
		
		# bail if there are no more uniq accs to sample
		delete $randelims->{$randacc};
		my $randleft = keys %$randelims;
		last if ($randleft == 0);

		# disallow random selections in the sister file
		if (exists ($data->{$hla_cat}->{$hla}->{$randacc})){ # or (exists ($rands->{$randacc}))){
		    next;
		}

		# disallow repeated random selections
		elsif (exists ($rands->{$randacc})){
		    next;
		}

		# go for it
		else {
		    print R ">$randacc-$hla\n$seqs->{$randacc}\n";
		    $rands->{$randacc} = 1;
		    last;
		}
	    }
	}
    }
}

#    my $randchrom = int(rand($counter)) + 1;
#    my $strand = int(rand(2)) + 1; # 1 is forward; 2 is reverse
#    my $end = int rand($lengths{$ids{$randchrom}}) + 1;
