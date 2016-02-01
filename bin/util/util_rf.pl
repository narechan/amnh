#!/usr/bin/perl

use Getopt::Std;

my %opts = ();
getopts ('g:t:o:', \%opts);
my $genetreesdir  = $opts{'g'};
my $tetree        = $opts{'t'};
my $outdir        = $opts{'o'};

# store the gene tree files in an array
opendir (D, "$genetreesdir");
my @genetreefiles = sort (readdir (D));
shift @genetreefiles;
shift @genetreefiles;
closedir (D);

# foreach gene tree, calculate its RF distances from the 
# TE tree using hashrf
foreach my $genetreefile (@genetreefiles){

    # run the hashrf algorithm
    `cat $tetree $genetreesdir/$genetreefile > tempcat.tre`;
    `hashrf tempcat.tre 2 -p list > tempcat.rfout`;
    
    # parse the hashrf output
    open (F, "tempcat.rfout");
    while (my $line = <F>){
	chomp $line;
	
	if ($line =~m/^\<1\,0\>/){
	    my ($comp, $rfdist) = split (/\s/, $line);
	    print "$genetreefile\t$rfdist\n";
	}
	else {
	    next;
	}
    }
    close (F);
    
    # cleanup
    `rm tempcat.*`;
}
