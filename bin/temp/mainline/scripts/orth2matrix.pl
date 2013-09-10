#!/usr/bin/env perl
#
# orth2matrix.pl
#
# Copyright (C) 2006-2011 Ernest K. Lee
#
#    This file is part of OrthologID.
#
#    OrthologID is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    OrthologID is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with OrthologID.  If not, see <http://www.gnu.org/licenses/>.
#
# Generate matrix from OrthologID orthologs files.
#
# Author: Ernest K Lee <elee@amnh.org>
#
#

use warnings;

my ($OID_HOME, $OID_USER_DIR);

BEGIN {
	$OID_HOME = $ENV{'OID_HOME'};
	die "Environment variable OID_HOME is not defined ... exiting.\n"
	  if !defined($OID_HOME);
	$OID_USER_DIR = $ENV{'OID_USER_DIR'};
	die "Environment variable OID_USER_DIR is not defined ... exiting.\n"
	  if !defined($OID_USER_DIR);

}

use Cwd;
use Getopt::Std;
use lib "$OID_HOME/lib";
use OrthologID;
use OrthologID::Utils;
use strict;

# defaults
my $minNonMissingTaxa = 4;
my %excludeTaxon;   # Exclude texon (key) if defined
my $needOutgroup;   # At least one outgroup must be present in each partition if defined
my $matrixFile = "$OID_USER_DIR/matrix.nex";
my $alignedFile = "FAMILY.aligned";


$Getopt::Std::STANDARD_HELP_VERSION = 1;
sub HELP_MESSAGE() {
	print "Usage: orth2matrix.pl [ -O ] [ -m min_non_missing_taxa ] [ -x comma_delimited_taxon_list_to_exclude ]\n";
}

# Get options
our ($opt_m, $opt_O, $opt_x);
getopts('Om:x:');
$minNonMissingTaxa = $opt_m if $opt_m;
$needOutgroup = 1 if $opt_O;
if ($opt_x) {
	foreach (split(',', $opt_x)) {
		$excludeTaxon{$_} = 1;
		print "Excluding $_ from matrix\n";
	}
}

# Check if matrix.nex already exists
die "\"$matrixFile\" already exists ... will not overwrite\n" if -f $matrixFile;

# Get taxa
my @INGROUP  = getIngroup();
my @OUTGROUP = getOutgroup();
my %isOutgroup;
my @incINGROUP;  # ingroup not excluded
my @incOUTGROUP; # outgroup not excluded
my %taxa;    # taxa to include in matrix if defined
foreach ( @INGROUP ) {
	if (!$excludeTaxon{$_}) {
		$taxa{$_} = 1;
		push @incINGROUP, $_;
	}
}
foreach ( @OUTGROUP ) {
	$isOutgroup{$_} = 1;
	if (!$excludeTaxon{$_}) {
		$taxa{$_} = 1;
		push @incOUTGROUP, $_;
	}
}
my @ogroup  = ();    # list of hashes of alignments
my @oanno = ();      # list of annotations

my $savDir = getcwd;
my $OID_DATADIR = "$OID_USER_DIR/data";
chdir $OID_DATADIR;
my $ogNum = 1;
foreach my $dir (<[1-9]*>) {
	
	my $hasOutgroup = 0;
	next if ! -r "$dir/orthologs";
	print "Processing family $dir\n";
	my $seqH = readFasta("$dir/$alignedFile");
	die "Fail to read/parse alignment.\n" if ! defined $seqH;
		
	# Determine if there is a global set of outgroup orthologs
	open OFH, "$dir/orthologs" or die "Cannot open orthologs file in \"$dir\": $!\n";
	my @globalOutgroup = ();
	my %globalOutgroupSeq = ();
	my $totGene = 0;
	while (<OFH>) {
		chomp;
		my @olist = split;
		my @taxlist = map { /^[^#]*/; $& } @olist;
		$totGene += @olist;
		my $allOutgroup = 1;       # only outgroup in this set
		my $allIngroup = 1;        # only ingroup in this set
		my %seq;
		foreach my $gene (@olist) {

			$gene =~ /#/;
			my $sp = $`;
			$seq{$gene} = $$seqH{$gene};

			if ($isOutgroup{$sp}) {
				$allIngroup = 0;
			}
			else {
				$allOutgroup = 0;
			}
		}
		if (! $allIngroup && ! $allOutgroup) {
			# do not use global outgroup
			@globalOutgroup = ();
			last;
		}
		if ($allOutgroup && @olist > @globalOutgroup) {  # pick the bigger set
			@globalOutgroup = @olist;
			%globalOutgroupSeq = %seq;
		}
	}
	close OFH;
	if (@globalOutgroup && $totGene == keys(%$seqH)) {
		#print " global outgroup = { @globalOutgroup }\n";
	}
	else {
		@globalOutgroup = ();
	}
	
	open OFH, "$dir/orthologs" or die "Cannot open orthologs file in \"$dir\": $!\n";
	while (<OFH>) {
		chomp;
		my @olist = split;
		my @taxlist = map { /^[^#]*/; $& } @olist;
		print " set = { @taxlist }\n";
 		my %numSp = ();
		my %seq;
		my $allOutgroup = 1;       # only outgroup in this set
		foreach my $gene (@olist) {

			$gene =~ /#/;
			my $sp = $`;
			my $geneID = $';
			$seq{$gene} = $$seqH{$gene};

			# $seq{$sp."#".$geneID} = $$seqH{$gene};
			$numSp{$sp}++ if $taxa{$sp};
			if ($isOutgroup{$sp}) {
				$hasOutgroup = 1;
			}
			else {
				$allOutgroup = 0;
			}
		}

		# count no. species in global outgroup
		my $numSpGlobalOutgroup = 0;
		my %numSpGOG;
		foreach my $gene (@globalOutgroup) {
			$gene =~ /#/;
			$numSpGOG{$`} ++;
		}
		
		if ($allOutgroup && @globalOutgroup > 0 ||
			keys(%numSp) + keys(%numSpGOG) < $minNonMissingTaxa ||
			$needOutgroup && ! $hasOutgroup) {
			print " ...skipped\n"
		}
		else {
			# store partition and annotation
			if (@globalOutgroup) {
				print " +global outgroup\n";
				while (my ($g,$s) = each %globalOutgroupSeq) {
					$seq{$g} = $s; 
				}
				push(@olist, @globalOutgroup);
			}
			push( @ogroup, \%seq);
			push( @oanno, "OG$ogNum|Family$dir|@olist");
			print "OG$ogNum: @olist\n";
			$ogNum++;
		}
	}
	
	close OFH;
}

chdir $savDir;

print "Generating matrix from " . scalar(@ogroup) . " ortholog groups ...\n";

# Sort taxa and generate matrix
my @taxaSorted = sort keys %taxa;
my $matrix = genSuperMatrix( \@taxaSorted, \@oanno, @ogroup );

open MFH, ">$matrixFile" or die "Cannot open matrix file for writing: $!\n";
print MFH"$matrix\n";

# Add partition size charsets

# Define outgroup
if ( @incOUTGROUP > 0 ) {
	print MFH "BEGIN PAUP;\n";
	print MFH "outgroup @incOUTGROUP;\n";
	print MFH "END;\n";
}

close MFH;

