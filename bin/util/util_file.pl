#!/usr/bin/perl

use Getopt::Std;

my %opts = ();
getopts ('d:', \%opts);
my $dir  = $opts{'d'};

opendir (D, "$dir");
my @resdirs = sort(readdir (D));
shift @resdirs;
shift @resdirs;
closedir (D);

foreach my $resdir (@resdirs){
    my $pileup = file_count ("$dir/$resdir/results/$resdir.pileup");
    my $indels = file_count ("$dir/$resdir/results/$resdir.indels2");
    my $snps = file_count ("$dir/$resdir/results/$resdir.snpsonly2");
    my $snpsf = file_count ("$dir/$resdir/results/$resdir.snpsfiltered2");
    
    my ($strain, $ref) = split (/\-/, $resdir);

    print "$strain\t$ref\t$pileup\t$indels\t$snps\t$snpsf\n";
}


###subs###

sub file_count{
    my $file = shift;

    my $counter = 0;
    open (F, "$file");
    while (my $line = <F>){
	chomp $line;
	$counter++;
    }
    close (F);
    
    return ($counter);
}
