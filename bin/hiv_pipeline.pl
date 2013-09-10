#!/usr/bin/perl -w

=head1 NAME

hiv_pipeline.pl

=head1 SYNOPSIS

hiv_pipeline.pl

Options:

--infile is your input data file
--ars is your arlequin config (ars file)
--outdir is your output dir
--procs is the number of cores you want to thread

Requires the bioperl libs. 
SNAP.pl and codons-xyplot.pl must be in your path
arlecore must be in your path (arlecore3513_32bit)

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
use Cwd;
use Parallel::ForkManager;

my ($help, $infile, $outdir, $ars, $procs);
GetOptions(
    'h|help'          => \$help,
    'i|infile=s'      => \$infile,
    'o|outdir=s'      => \$outdir,
    'a|ars=s'         => \$ars,
    'p|procs=s'       => \$procs,	   
	   ) or pod2usage;
pod2usage if $help;

for my $option ($infile, $outdir, $ars, $procs){
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
    
    $seq =~s/[BDEFHIJKLMOPQRSUVWXYZx]/N/g;
    
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

# write out the HLA population files and all others
foreach my $hla_cat (sort keys %$data){
    my $pm = Parallel::ForkManager->new($procs);
    foreach my $hla (sort keys %{$data->{$hla_cat}}){
	$pm->start and next;
	
	# data structure for the final table
	my $datastruc = {};
	
	open (W, ">$outdir/$hla.tabseq");
	open (WF, ">$outdir/$hla.fasta");
	open (R, ">$outdir/$hla-other.tabseq");
	open (RF, ">$outdir/$hla-other.fasta");
	
	# data structure to track which accs have been harvested                                   
        my $randelims = {};
        foreach my $acc (@accs){
            $randelims->{$acc} = 1;
        }
	
	# print the data for this HLA
	foreach my $acc (sort keys %{$data->{$hla_cat}->{$hla}}){
            print W "$acc-$hla\t$data->{$hla_cat}->{$hla}->{$acc}\n";
	    print WF ">$acc-$hla\n$data->{$hla_cat}->{$hla}->{$acc}\n";
	    delete $randelims->{$acc};
	}
	
	# print out the leftover accessions
	foreach my $leftover (sort keys %$randelims){
	    print R "$leftover-$hla\t$seqs->{$leftover}\n";
	    print RF ">$leftover-$hla\n$seqs->{$leftover}\n";
	}
	close (W);
	close (R);
	close (WF);
	close (RF);

	# do the snap analysis for these files
	my $pwd = getcwd;
	
	`mkdir $outdir/$hla`;
	chdir ("$outdir/$hla");
	print STDERR "SNAP $hla\n";
	`SNAP.pl ../$hla.tabseq`;
	`codons-xyplot.pl codons.*`;

	chdir ("$pwd");

	`mkdir $outdir/$hla-other`;
        chdir ("$outdir/$hla-other");
	print STDERR "SNAP $hla-other\n";
	`SNAP.pl ../$hla-other.tabseq`; 
	`codons-xyplot.pl codons.*`;
	chdir ("$pwd");

	# parse and calculate dn/ds
	print STDERR "Calculating dN/dS $hla\n";
	for my $dir ("$outdir/$hla/", "$outdir/$hla-other/"){
	    open (DNDS, ">$dir/dnds.out");
	    open (F, "$dir/codon.data");
	    while (my $line = <F>){
		chomp $line;
		next if ($line =~m/^\#/);
		
		my @line = split (/\s+/, $line);
		shift @line;

		my $dnds;
		if (($line[5] == 0) and ($line[6] > 0)){
		    print DNDS "$line[0]\t100\n";
		    $dnds = 100;
		}
		elsif (($line[5] == 0) and ($line[6] == 0)){
		    print DNDS "$line[0]\t0\n";
		    $dnds = 0;
		}
		else{
		    $dnds = $line[6] / $line[5];
		    print DNDS "$line[0]\t$dnds\n";
		}

		if ($dir =~m/other/){
		    if ($dnds > 1){
			$datastruc->{$line[0]}->{'nonexp'} = 1;
		    }
		    else {
			$datastruc->{$line[0]}->{'nonexp'} = 0;
		    }
		}
		else {
		    if ($dnds > 1){
                        $datastruc->{$line[0]}->{'exp'} = 1;
                    }
                    else {
                        $datastruc->{$line[0]}->{'exp'} = 0;
                    }
		}
	    }
	    close (DNDS);
	    close (F);
	}
	
	# generate the translated data file for arlequin
	print STDERR "Translation $hla\n";
	my $alnlen = 0;
	my $rseqcount = 0;
	my $realpts = {};
	my $realseqin = Bio::SeqIO->new (-format=>'Fasta', -file=>"$outdir/$hla.fasta");
	while (my $sequence_obj = $realseqin->next_seq()){
	    $rseqcount++;
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

	my $oseqcount = 0;
	my $randpts = {};
	my $randseqin = Bio::SeqIO->new (-format=>'Fasta', -file=>"$outdir/$hla-other.fasta");
	while (my $sequence_obj = $randseqin->next_seq()){
	    $oseqcount++;
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
	print STDERR "Arlequin $hla\n";
	`mkdir -p $outdir/$hla-arl`;
	foreach my $pos (sort {$a <=> $b} keys (%$realpts)){
	    print STDERR "Working on Arlequin calc for $hla $pos\n";

	    open (AR, ">$outdir/$hla-arl/$hla.$pos.arp");
	    print AR "[Profile]\n\n";
	    print AR "Title=\"$hla.$pos\"\n";
	    print AR "NbSamples=2\n";
	    print AR "DataType=FREQUENCY\n";
	    print AR "GenotypicData=0\n";
	    print AR "LocusSeparator=WHITESPACE\n";
	    print AR "MissingData='?'\n";
	    print AR "Frequency=REL\n";
	    print AR "\n";
	    print AR "[DATA]\n\n";
	    print AR "[[SAMPLES]]\n";
	    print AR "SampleName=\"$hla.$pos\"\n";
	    print AR "SampleSize=$rseqcount\n";
	    print AR "SampleData={\n";
	   
	    my $realfreqs = {};
	    foreach my $aa (sort {$aakey{$a} <=> $aakey{$b}} keys %aakey){
		my $freq = 0;
		if (exists ($realpts->{$pos}->{$aa})){
		    $freq = $realpts->{$pos}->{$aa} / $rseqcount;
		    print AR "$aakey{$aa} $freq\n";
		    $realfreqs->{$pos}->{$aa} = $freq;
		}
		else {
		    print AR "$aakey{$aa} $freq\n";
		    $realfreqs->{$pos}->{$aa} = $freq;
		}
	    }
	    print AR "}\n";
	    
	    print AR "SampleName=\"$hla-other.$pos\"\n";
	    print AR "SampleSize=$oseqcount\n";
	    print AR "SampleData={\n";
	    
	    my $randfreqs = {};
	    foreach my $aa (sort {$aakey{$a} <=> $aakey{$b}} keys %aakey){
		my $freq = 0;
		if (exists ($randpts->{$pos}->{$aa})){
		    $freq = $randpts->{$pos}->{$aa} / $oseqcount;
		    print AR "$aakey{$aa} $freq\n";
		    $randfreqs->{$pos}->{$aa} = $freq;
		}
		else {
		    print AR "$aakey{$aa} $freq\n";
		    $randfreqs->{$pos}->{$aa} = $freq;
		}
	    }
	    print AR "}\n\n";
	    
	    print AR "[[STRUCTURE]]\n";
	    print AR "StructureName=\"$hla\"\n";
	    print AR "NbGroups=1\n";
	    print AR "Group={\n";
	    print AR "\"$hla.$pos\"\n";
	    print AR "\"$hla-other.$pos\"\n";
	    print AR "}\n";
	    close (AR);
	    
	    # perform arlequin calculation on this position 
	    # if the frequencies are diff, otherwise bail
	    my $compare = Compare ($realfreqs->{$pos}, $randfreqs->{$pos});
	    
            if ($compare == 1){
                print STDERR "Aborted: populations are the same!\n";
		$datastruc->{$pos}->{'fst'} = "NA";
		$datastruc->{$pos}->{'exact'} = "NA";
                next;
            }
            else {
		`arlecore3513_32bit $outdir/$hla-arl/$hla.$pos.arp $ars`;
	    }
	    
	    # parse the arlequin output
	    my $filename = $hla . "." . $pos . "." . "htm";
	    open (RA, "$outdir/$hla-arl/$hla.$pos.res/$filename");
	    while (my $line = <RA>){
		chomp $line;
		if ($line =~m/P\-value\s*=\s*(.+)\+/i){
		    if ($1 <= 0.05){
			$datastruc->{$pos}->{'fst'} = 1;
		    }
		    else {
			$datastruc->{$pos}->{'fst'} = 0;
		    }
		}
		elsif ($line =~m/P\svalue\s*=\s*(.+)\s*\+/i){
		    if ($1 <= 0.05){
                        $datastruc->{$pos}->{'exact'} = 1;
                    }
                    else {
                        $datastruc->{$pos}->{'exact'} = 0;
                    }
		}
		else {
		    next;
		}
	    }
	}
	
	# print out the final table                                                                           
	open (TB, ">$outdir/$hla.out");
	foreach my $p (sort {$a <=> $b} keys %$datastruc){
	    my @array;
	    for my $type ("nonexp", "exp", "fst", "exact"){
		push (@array, $datastruc->{$p}->{$type});
	    }
	    my $divergentflag;
	    if (($datastruc->{$p}->{'nonexp'} == 1) or ($datastruc->{$p}->{'nonexp'} == 1)){
		if (($datastruc->{$p}->{'fst'} == 1) and ($datastruc->{$p}->{'exact'} == 1)){
		    $divergentflag = 1;
		}
		else {
		    $divergentflag = 0;
		}
	    }
	    else{
		$divergentflag = 0
	    }

	    my $string = join "\t", @array;
	    print TB "$p\t$string\t$divergentflag\n";
	}
	close (TB);
	$pm->finish;
    }
    $pm->wait_all_children;
}


	





#### OLD CODE ####
=head
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
=cut
#    my $randchrom = int(rand($counter)) + 1;
#    my $strand = int(rand(2)) + 1; # 1 is forward; 2 is reverse
#    my $end = int rand($lengths{$ids{$randchrom}}) + 1;
