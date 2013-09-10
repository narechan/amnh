#!/usr/bin/perl

use Getopt::Std;

my %opts = ();
getopts ('i:', \%opts);
my $indir    = $opts{'i'};

opendir (D, "$indir");
my @exptfiles = sort (readdir(D));
shift @exptfiles;
shift @exptfiles;
closedir (E);

my $data = {};
foreach my $file (@exptfiles){
    my @file = split (/\./, $file);
    my $suffix = pop (@file);
    next unless ($suffix eq "log");
    
    my $svp;
    my $paralog;
    my $ignore;
    my $class;
    my $seq1;
    my $seq2;
    my $icount;
    my $mmrate;
    my $sim;
    my $score;
    open (F, "$indir/$file");
    while (my $line = <F>){
	chomp $line;
	
	# get the designation
	($splice  = $1) if ($line =~m/Number of splice variant pairs\s+\:\s(\d)/i);
	($paralog = $1) if ($line =~m/Number of none splice variants pairs\:\s(\d)/i);
	($ignored = $1) if ($line =~m/Number of ignored sequence pairs\s+\:\s(\d)/i);
	if ($splice == 1){
	    $class = "SpliceVar";
	}
	elsif ($paralog == 1){
	    $class = "Paralog";
	}
	elsif ($ignored == 1){
	    $class = "Ignored";
	}
	else {
	    $class = "Unknown";
	}
	
	# get the data
	($seq1 = $1) if ($line =~m/Sequence\s\#1\s\(\'(.*)\s\'\)/);
	($seq2 = $1) if ($line =~m/\s+\#2\s\(\'(.*)\s\'\)/);
	($icount = $1) if (($line =~m/inverseCBINcount\s([\d.-]+)/) or ($line =~m/Mean\-CBIN\-Len\s([\d.-]+)/));
	($mmrate = $1) if (($line =~m/MatchMismatchRate\s([\d.-]+)/) or ($line =~m/Match\-Rate\s([\d.-]+)/));
	($sim = $1) if (($line =~m/Similarity\s([\d.-]+)/) or ($line =~m/SeqID\s([\d.-]+)/));
	($score = $1) if ($line =~m/Score\:\s([\d.-]+)/);
#	if ($line =~m/\s+\#2\s\(\'(.*)\s\'\)\s+\-\>\s+inverseCBINcount\s(.*)\s\/\sMatchMismatchRate\s(.*)\s\/\sSimilarity\s(.*)\s+\-\>\s+\[Score\:\s(.*)\]/){
#	    $seq2 = $1;
#	    $icount = $2;
#	    $mmrate = $3;
#	    $sim    = $4;
#	    $sim=~s/\s//g;
#	    $score  = $5;
#	}
    }
    close (F);

    print "$seq1\t$seq2\t$class\t$icount\t$mmrate\t$sim\t$score\n";
}
