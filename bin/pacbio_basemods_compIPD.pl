#!/usr/bin/perl -w

#####SETUP#####
use strict;
use warnings;
use Chart::Gnuplot;
use Statistics::Descriptive;
use Statistics::Normality ':all';
use Statistics::Test::WilcoxonRankSum;
use Statistics::TTest;
use Getopt::Long;
use Pod::Usage;

my ($help, $rscript, $outdir, $posfile1, $posfile2, $infile1, $infile2, $ipdrcut);
GetOptions(
    'h|help'          => \$help,
    'r|rscript=s'     => \$rscript, #path to the R script that mines IPDs from pacbio hdf5 files
    'o|outdir=s'      => \$outdir, #output lives here (dist plots)
    'x|posfile1=s'   => \$posfile1, #position file you want to interrogate
    'y|posfile2=s'   => \$posfile2, #position file you want to interrogate 
    'a|infile1=s'     => \$infile1, # *.comp.h5 file path 
    'b|infile2=s'     => \$infile2, # the second h5 file
    'z|ipdrcut=s'    => \$ipdrcut # ipdr cutoff for paying attention to the site
	   ) or pod2usage;
pod2usage if $help;

`mkdir -p $outdir`;

#####MAIN#####
my $undef1 = {};
my $undef2 = {};
open (X1, "$posfile1");
while (my $line1 = <X1>){
    chomp $line1;
    next if ($line1 =~m/^coord/);
    
    # note the for strands, 1 is F/+ and 2 is R/-
    my ($coord1, $position1, $seq1, $strand1, $cov1, $qv1, $ipdr1) = split (/\t/, $line1);
    next unless ($position1 == 4056302);
    # log if ipdr meets cutoff but there is no reference placement
    if ($position1 eq "N"){
	($undef1->{$coord1} = $ipdr1) if ($ipdr1 < $ipdrcut);
	next;
    }
    
    open (Y1, "$posfile2");
    while (my $line2 = <Y1>){
	chomp $line2;
	next if ($line2 =~m/^coord/);
	my ($coord2, $position2, $seq2, $strand2, $cov2, $qv2, $ipdr2) = split (/\t/, $line2);
	next unless ($position2 == 4056302);
	# log if ipdr meets cutoff but there is no reference placement                                 
	if ($position2 eq "N"){
	    ($undef2->{$coord2} = $ipdr2) if ($ipdr2 < $ipdrcut);
	    next;
	}
	
	# bail if we're not looking at the same position
	next if ($position1 != $position2);

	# bail if we are not looking at the same strand
#	next if ($strand1 != $strand2);
	
	# bail and log if both ipdrs are above the threshold we want to analyze
	if (($ipdr1 > $ipdrcut) and ($ipdr2 > $ipdrcut)){
	    print "$coord1\t$position1\t$strand1\t$ipdr1\t\t";
	    print "$coord2\t$position2\t$strand2\t$ipdr2\t\t";
	    print "SKIP_IPDR\n";
	    next;
	}
	
        # mine the ipds using the R script
	# note that the strands here work because, first, the strands in the csv file are
	# 0 for forward and 1 for reverse. we need to look at ipds on the opposite strand
	# (1 for forward and 0 for reverse) which is the input switch
	print STDERR "Mining $infile1:$coord1:$position1:$strand1\n";
	my $ipds1 = ipd_mine($coord1, $rscript, $infile1, $outdir, $strand1);
	print STDERR "Mining $infile2:$coord2:$position2:$strand2\n";
	my $ipds2 = ipd_mine($coord2, $rscript, $infile2, $outdir, $strand2);

	my $ipds1count = @$ipds1;
	my $ipds2count = @$ipds2;
	print STDERR "ipdcount1:$ipds1count\n";
	print STDERR "ipdcount2:$ipds2count\n";
  	
        # get basic stats on the distributions
	print STDERR "Generating basic stats and distribution plot\n";

	open (S1, ">$outdir/$infile1\_$position1\_$strand1.stats");
	my $statobj1 = Statistics::Descriptive::Full->new();
	$statobj1->add_data(@$ipds1);
	my $mean1 =  sprintf("%.3f", $statobj1->mean());
	my $max1  = $statobj1->max();
	my $min1  = $statobj1->min();
	my %freqs1 = $statobj1->frequency_distribution(10);
	my ($pval_w1, $w_statistic1) = shapiro_wilk_test ($ipds1);
	my ($pval_k1, $k_statistic1) = dagostino_k_square_test ($ipds1); 
	print S1 "Mean:$mean1\nMax:$max1\nMin:$min1\n";
	print S1 "SW:$w_statistic1($pval_w1)\n";
	print S1 "K2:$k_statistic1($pval_k1)\n";
	close (S1);

	open (S2, ">$outdir/$infile2\_$position2\_$strand2.stats");
	my $statobj2 = Statistics::Descriptive::Full->new();
	$statobj2->add_data(@$ipds2);
	my $mean2 =  sprintf("%.3f", $statobj2->mean());
	my $max2  = $statobj2->max();
	my $min2  = $statobj2->min();
	my %freqs2 = $statobj2->frequency_distribution(10);
	my ($pval_w2, $w_statistic2) =  shapiro_wilk_test ($ipds2);
	my ($pval_k2, $k_statistic2) =  dagostino_k_square_test ($ipds2);
	print S2 "Mean:$mean2\nMax:$max2\nMin:$min2\n";
	print S2 "SW:$w_statistic2($pval_w2)\n";
	print S2 "K2:$k_statistic2($pval_k2)\n";
	close (S2);
	
	open (SALL, ">$outdir/$infile1\_$position1\_$strand1-v-$infile2\_$position2\_$strand2.stats");
	my $wilcox_test = Statistics::Test::WilcoxonRankSum->new();
	$wilcox_test->load_data($ipds1, $ipds2);
	my $wilcoxprob =  sprintf("%.3f", $wilcox_test->probability());
	print SALL "Wilcoxon:$wilcoxprob\n";
	
	my $ttest = new Statistics::TTest;
	$ttest->load_data($ipds1, $ipds2);
	my $ttestprob =  sprintf("%.3f", $ttest->{t_prob});
	print SALL "TTest:$ttestprob\n";
	close (SALL);
	
	print "$coord1\t$position1\t$strand1\t$ipdr1\t$mean1\t";
	print "$coord2\t$position2\t$strand2\t$ipdr2\t$mean2\t";
	print "ANALYZE_IPDR\t";
	print "$pval_w1\t$pval_k1\t$pval_w2\t$pval_k2\t$wilcoxprob\t$ttestprob\n";

        # plot the frequency distribution
	my $xmax = 0;
	my $xmin = 0;
	my $ymax = 0;

	my @yfreq1;
	my @xfreq1;
	foreach my $bin (sort {$a <=> $b} keys %freqs1){
	    if ($bin > $xmax){
		$xmax = $bin;
	    }
	    if ($bin < $xmin){
		$xmin = $bin;
	    }
	    push (@xfreq1, $bin);
	    if ($freqs1{$bin} > $ymax){
		$ymax = $freqs1{$bin};
	    }
	    push (@yfreq1, $freqs1{$bin});
	}
	
	my @xfreq2;
	my @yfreq2;
	foreach my $bin (sort {$a <=> $b} keys %freqs2){
	    if ($bin > $xmax){
		$xmax = $bin;
	    }
	    if ($bin < $xmin){
		$xmin = $bin;
	    }
	    push (@xfreq2, $bin);
	    if ($freqs2{$bin} > $ymax){
		$ymax = $freqs2{$bin};
	    }
	    push (@yfreq2, $freqs2{$bin});
	}
	
	
	my $chart = Chart::Gnuplot->new(output => "$outdir/$coord1\_$position1\_$strand1-v-$coord2\_$position2\_$strand2.gif",
					xrange=>[$xmin - 0.1, $xmax + 0.1],
					yrange=>[0, $ymax + 1],
					xlabel=>"IPD",
					ylabel=>"Freq",
					title=>"IPD Freq",
					bg=>"white");
	my $data1 = Chart::Gnuplot::DataSet->new(xdata=>\@xfreq1,
						 ydata=>\@yfreq1,
						 style=>"linespoints",
						 pointtype=>"fill-circle",
						 pointsize=>0.5);
	my $data2 = Chart::Gnuplot::DataSet->new(xdata=>\@xfreq2,
						 ydata=>\@yfreq2,
						 style=>"linespoints",
						 pointtype=>"fill-circle",
						 pointsize=>0.5);
	$chart->plot2d($data1, $data2);
    }
}

open (EE, ">$outdir/undefs");
foreach my $c1 (sort {$a <=> $b} keys %$undef1){
    print EE "$c1\t$undef1->{$c1}\n";
} 
foreach my $c2 (sort {$a <=> $b} keys %$undef2){
    print EE "$c2\t$undef2->{$c2}\n";
}
close (EE);

=head
# write the raw ipd file for the position
open (O, ">$outdir/$infile1\_$position1\_$strand1.out");
foreach my $ipd (@$ipds1){
    print O "$ipd\n";
}
close (O);

open (N, ">$outdir/$infile2\_$position2\_$strand2.out");
foreach my $ipd (@$ipds2){
    print N "$ipd\n";
}
close (N);
=cut

sub ipd_mine{
    my $pos = shift;
    my $rscrip = shift;
    my $in = shift;
    my $out = shift;
    my $str = shift;
    
    # run the R script
    `R --slave --args $pos $str $in $outdir/$in\_$pos\_$str.data < $rscrip`;
    
    # parse the output for the positon and strand we want
    my @ret;
    my @ipds;
    open (F, "$outdir/$in\_$pos\_$str.data");
    while (my $line = <F>){
	chomp $line;
	
	my ($refpos, $readbase, $refbase, $idx, $st, $ipd) = split ("\t", $line);
	(push (@ret, log10 ($ipd + 0.01))) unless ($ipd eq "NA");
#	(push (@ret, $ipd)) unless ($ipd eq "NA");
#	(print "$ipd\n") unless ($ipd eq "NA");
    }
    close (F);
#    print "\n";

    # remove the temp file
    `rm $outdir/$in\_$pos\_$str.data`;
    
    return (\@ret)
}

sub log10 {
    my $n = shift;
    return log($n)/log(10);
}
