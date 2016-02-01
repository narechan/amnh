#!/usr/bin/perl

use Getopt::Std;

my %opts = ();
getopts ('g:t:o:h:l:r:d:', \%opts);
my $genelist      = $opts{'l'};
my $genetreesdir  = $opts{'g'};
my $genehmmsdir   = $opts{'h'};
my $tetree        = $opts{'t'};
my $outdir        = $opts{'o'};
my $refptfile     = $opts{'r'};
my $refptdata     = $opts{'d'};

`mkdir -p $outdir/hmmsearch`;

# store the gene tree files in an array
opendir (D, "$genetreesdir");
my @genetreefiles = sort (readdir (D));
shift @genetreefiles;
shift @genetreefiles;
closedir (D);

# store the hmms in an array
opendir (H, "$genehmmsdir");
my @hmms = sort (readdir (D));
shift @hmms;
shift @hmms;
closedir (D);

# parse the reference sequence annotations
$data = {};
open (R, "$refptdata");
while (my $line = <R>){
    chomp $line;
    my ($acc, $start, $end, $strand, $desc) = 
	split (/\t/, $line);
    $data->{$acc}->{'start'} = $start;
    $data->{$acc}->{'end'}   = $end;
}
close (R);

# foreach gene tree calculate its RF distances from the 
# TE tree using hashrf, and calculate each gene's homolog in 
# the reference using hmmsearch
open (L, "$genelist");
while (my $gene = <L>){
    chomp $gene;

    # run the hashrf algorithm
    `cat $tetree $genetreesdir/RAxML_bestTree.$gene > tempcat.tre`;
    `hashrf tempcat.tre 2 -p list > tempcat.rfout`;
    
    # parse the hashrf output
    my $comp;
    my $rfdist;
    open (F, "tempcat.rfout");
    while (my $line = <F>){
	chomp $line;
	
	if ($line =~m/^\<1\,0\>/){
	    ($comp, $rfdist) = split (/\s/, $line);
#	    print "$genetreefile\t$rfdist\n";
	}
	else {
	    next;
	}
    }
    close (F);
 
     # cleanup hashrf                                                                                         
    `rm tempcat.*`;

    # run hmmsearch
    `hmmsearch --notextw --tblout $outdir/hmmsearch/$gene.out $geneshmmdir/$gene.hmm $refptfile`;

    # parse hmmsearch
    my $refseq;
    my $refscore;
    open (A, "$outdir/hmmsearch/$boid.hmm.tblout");
    while (my $line = <A>){
        chomp $line;
        next if ($line =~m/\#/);
        my @line = split (/\s+/, $line);

        if ($line[5] > $score){
            $refseq = $line[0];
            $refscore = $line[5];
        }
        else {
            $refseq = "NONE_SIG";
            $refscore ="NA";
        }
        last; #just want the top hit                                                 
    }
    close (A);
    
    # print final report
    print "$gene\t$rfdist\t$refseq\t$refscore\t$data->{$refseq}->{'start'}\t$data->{$refseq}->{'end'}\n";
}
