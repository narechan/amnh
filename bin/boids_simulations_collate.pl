#!/usr/bin/perl

use Getopt::Std;
use Statistics::Descriptive;

my %opts = ();
getopts ('d:', \%opts);
my $dir  = $opts{'d'};

opendir (D, "$dir");
my @repdirs = sort(readdir (D));
shift @repdirs;
shift @repdirs;
close (D);

my $data = {};
foreach my $repdir (@repdirs){
    opendir (T, "$dir/$repdir");
    my @treedirs = sort(readdir (T));
    shift @treedirs;
    shift @treedirs;
    closedir (T);
    
    foreach my $treedir (@treedirs){
	open (F, "$dir/$repdir/$treedir/cf.mds.final.stats");
#	open (F, "$dir/$repdir/$treedir/pam.clusters.radmean");
#	open (F, "$dir/$repdir/$treedir/hclust.clusters.radmean");
	while (my $line = <F>){
	    chomp $line;
	    my @line = split (/\t/, $line);
	    my $counter = 0;
	    foreach my $field (@line){
		$counter++;
		push (@{$data->{$treedir}->{$counter}}, $field);
	    }
	}
	close (F);
    }
}

foreach my $t (sort {$a <=> $b} keys %$data){
    print "$t\t";
    my @means;
    foreach my $col (sort {$a <=> $b} keys %{$data->{$t}}){
	my $statobj = Statistics::Descriptive::Full->new();
	$statobj->add_data(@{$data->{$t}->{$col}});
	my $mean = $statobj->mean();
	push (@means, $mean);
    }
    my $string = join "\t", @means;
    print "$string\n";
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
