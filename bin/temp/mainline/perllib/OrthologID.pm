#
# OrthologID Perl module
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
# NOTE:
# - The OID_HOME and OID_USER_DIR environment variables must be defined.
# - Gene names must be prefixed by species names as defined in config, in the format
#   species#geneID, e.g. Arath#At1g12345 or SpX#NP_123456.1 (DO NOT USE "-" in species/gene names).
# - The command line 'paup' and 'mcl' executables must be in your PATH.
#
# Author: Ernest K Lee <elee@amnh.org>
#

package OrthologID;

use 5.008005;
use strict;
use warnings;

use Cwd;
use DB_File;
use IO::Compress::Gzip qw(gzip $GzipError);
use OrthologID::PhyloTree;
use OrthologID::Utils;

require Exporter;

our @ISA = qw(Exporter);

our %EXPORT_TAGS = (
	'all' => [
		qw(

		  )
	]
);

our @EXPORT_OK = ( @{ $EXPORT_TAGS{'all'} } );

our @EXPORT = qw(
  printConfig
  getIngroup
  getOutgroup
  allBlast
  makeFamily
  alignFamily
  makeTree
  findOrthologs
);

our $VERSION = '0.01';

# Local variables
my $OID_HOME;        # OrthologID main directory
my $OID_USER_DIR;    # OrthologID user run directory
my $OID_CONF;        # OrthologID config file
my $OID_DATADIR;     # OrthologID guide tree data directory
my $OID_BLASTDIR;    # Directory to store OrthologID blast results
my $BLAST_HOME;      # NCBI BLAST installation directory
my @INGROUP;         # Ingroup taxa
my @OUTGROUP;        # Outgroup taxa
my $NCPU;            # Number of CPU to be used
my $VERBOSE;         # Verbosity {0,1,2}
my $DEBUG = 0;       # Debug mode
my $BLAST_RES_DB;    # Database of all-all blast
my ( $DB_USER, $DB_PASSWORD );    # MySQL DB username and password

# BLAST e-value cutoff for potential orthologs
my $FAM_E_CUTOFF     = "1e-10";


# OrthologID "global" variables
die "OID_HOME environment variable undefined!\n" if !$ENV{'OID_HOME'};
die "OID_USER_DIR environment variable undefined!\n" if !$ENV{'OID_USER_DIR'};
$OID_HOME       = $ENV{'OID_HOME'};
$OID_USER_DIR   = $ENV{'OID_USER_DIR'};
$OID_BLASTDIR   = "$OID_USER_DIR/blast";
$BLAST_RES_DB   = "$OID_BLASTDIR/blastres.dbm";
$OID_DATADIR    = "$OID_USER_DIR/data";
$ENV{'BLASTDB'} = "$OID_USER_DIR/blastdb";
$ENV{'PATH'}    = "$OID_HOME/bin:$ENV{'PATH'}";
$OID_CONF       = "$OID_USER_DIR/config";
die "OrthologID config file not found!\n" if !-r $OID_CONF;

# Clustering
my $GENE_LEN_DB = "$OID_BLASTDIR/genelen.dbm";
my $clustersFile = "$OID_BLASTDIR/clusters";

my $unalignedFamily = "FAMILY";
my $alignedFamily = "FAMILY.aligned";

# Initialize
&initOID;

sub initOID {

	# Defaults
	$NCPU = 1;

	# Parse configuration file
	open CONF_FH, "$OID_CONF" or die "Cannot open $OID_CONF: $!\n";
	while (<CONF_FH>) {
		chomp;
		next if (/^#/);
		if (/^\s*BLASTHOME\s*=\s*(.+?)\s*$/)
		{    # Installation directory of NCBI BLAST
			$BLAST_HOME = $1;
		}
		elsif (/^\s*INGROUP\s*=\s*(.+?)\s*$/) {    # Ingroup species prefixes
			$_       = $1;
			@INGROUP = split;
		}
		elsif (/^\s*OUTGROUP\s*=\s*(.+?)\s*$/) {    # Outgroup species prefixes
			$_        = $1;
			@OUTGROUP = split;
		}
		elsif (/^\s*NCPU\s*=\s*(\d+)/) {            # Number of CPUs available
			$NCPU = $1;
		}
		elsif (/^\s*VERBOSE\s*=\s*(\d+)/) {         # Verbosity
			$VERBOSE = $1;
		}
		elsif (/^\s*DEBUG\s*=\s*(\d+)/) {           # Debug mode
			$DEBUG = $1;
		}
		elsif (/^\s*(MYSQL_\w+)\s*=\s*(.+?)\s*$/) {    # MYSQL environment (not used)
			$ENV{$1} = $2;
		}
		elsif (/^\s*DB_USER\s*=\s*(\S+)\s*$/) {        # MYSQL user (not used)
			$DB_USER = $1;
		}
		elsif (/^\s*DB_PASSWORD\s*=\s*(\S+)\s*$/) {    # MYSQL password (not used)
			$DB_PASSWORD = $1;
		}
	}
	close CONF_FH;

	die if !( @INGROUP && @OUTGROUP );

	return 1;
}

# Preloaded methods go here.


#
# Retrieve ingroup list
#
sub getIngroup() {
	return @INGROUP;
}

#
# Retrieve outgroup list
#
sub getOutgroup() {
	return @OUTGROUP;
}

#
# Print configuration info
#
sub printConfig() {
	print "OID_HOME = $OID_HOME\n";
	print "OID_USER_DIR= $OID_USER_DIR\n";
	print "Ingroup taxa = @INGROUP\n";
	print "Outgroup taxa = @OUTGROUP\n";
}

#
# All-all blast
# Optional argument: Ingroup taxon prefix, e.g. allBlast("At")
# No argument means all ingroup taxa.
# Asterisk (*) means combine all previously generated blast results.
#
sub allBlast {
	my $prefix  = shift;
	my $verbose = $VERBOSE;
	#### E-value cutoff
	#my $eCutOff = "1e-16";
	my $eCutOff = $FAM_E_CUTOFF;
	####

	my $outputDir = $OID_BLASTDIR;    # Directory to store blast output
	my @queryGroup;                   # BLAST query taxa
	my $blastResDB;                   # DB file to store blast results
	my $geneLenDB;                    # DB file to store gene length
	my $numProc = $NCPU;
	my %blastRes;
	my %geneLen;

	if ( !-d $outputDir ) {
		mkdir $outputDir, 0755 or die "Cannot create directory: $!\n";
	}

	if ( defined($prefix) && $prefix ne '*' ) {
		$blastResDB = "$OID_BLASTDIR/$prefix" . "_blastres.dbm";
		$geneLenDB = "$OID_BLASTDIR/$prefix" . "_genelen.dbm";
		@queryGroup = ($prefix);
	}
	else {
		$blastResDB = $BLAST_RES_DB;
		$geneLenDB = $GENE_LEN_DB;
		#@queryGroup = (@INGROUP, @OUTGROUP);
		if ( $prefix eq '*' ) {

			# Combine previously generated blast results for all taxa
			warn "$blastResDB exists.  Previous results will be overwritten.\n"
			  if -f $blastResDB;
			tie %blastRes, "DB_File", $blastResDB
			  or die "Cannot open $blastResDB: $!\n";
			tie %geneLen, "DB_File", $geneLenDB
				or die "Cannot open $geneLenDB: $!\n";
			  
			foreach my $sp (@INGROUP, @OUTGROUP) {
				my %spBlastRes;
				my %spGeneLen;
				my $spBlastResDB = "$OID_BLASTDIR/$sp" . "_blastres.dbm";
				if ( !-f $spBlastResDB ) {
					print STDERR "Blast results for $sp not found!\n";
					next;
				}
				my $spGeneLenDB = "$OID_BLASTDIR/$sp" . "_genelen.dbm";
				if ( !-f $spGeneLenDB ) {
					print STDERR "Gene lengths for $sp not found!\n";
					next;
				}

				print "Combining BLAST and length results for $sp\n" if $verbose;
				tie %spBlastRes, "DB_File", $spBlastResDB;
				foreach my $locus ( keys %spBlastRes ) {
					$blastRes{$locus} = $spBlastRes{$locus};
				}
				untie %spBlastRes;
				tie %spGeneLen, "DB_File", $spGeneLenDB;
				foreach my $locus ( keys %spGeneLen ) {
					$geneLen{$locus} = $spGeneLen{$locus};
				}
				untie %spGeneLen;
			}
			untie %blastRes;
			untie %geneLen;
			return 1;
		}
	}

	my ($BLASTALL, $FASTACMD);
	$BLASTALL = defined($BLAST_HOME) ? "$BLAST_HOME/bin/blastall" : "blastall";
	$FASTACMD = defined($BLAST_HOME) ? "$BLAST_HOME/bin/fastacmd" : "fastacmd";

	# Read in BLAST database
	my @qdb = ();
	my @tmpdb;
	foreach ( 0 .. @queryGroup - 1 ) {
		if ($verbose) { print "Reading in ${queryGroup[$_]} sequences\n"; }
		chomp( @tmpdb = `$FASTACMD -d $queryGroup[$_] -D 1` );
		@qdb = ( @qdb, @tmpdb );
	}

	# Store all sequences in hash
	my %qSeq;
	my $locus;
	foreach (@qdb) {
		if (/^>\w*\|(\S+)/) {    # assume 2nd field is id
			$locus = $1;
			$qSeq{$locus} = "";
		}
		else {
			$qSeq{$locus} .= $_;
		}
	}

	# Blast each query sequence against everything
	my $tmpDir = "/tmp/oid-$$";
	if ( !-d $tmpDir ) {
		mkdir $tmpDir, 0755 or die "Cannot create $tmpDir: $!";
	}
	my $tmpFile = "$tmpDir/oid-blast-query";
	my ( $blastCmd, @blastOut );
	my $count;

	warn "$blastResDB exists.  Previous results will be overwritten.\n"
	  if -f $blastResDB;
	tie %blastRes, "DB_File", $blastResDB
	  or die "Cannot open $blastResDB: $!\n";
	tie %geneLen, "DB_File", $geneLenDB
	  or die "Cannot open $geneLenDB: $!\n";

	print "BLASTing ...\n" if $verbose;
	my $nq = keys %qSeq;
	foreach my $locus ( keys %qSeq ) {

		# Record length of sequence
		$geneLen{$locus} = length($qSeq{$locus});
		
		# Write query
		open FH, ">$tmpFile" or die "Cannot open tmp file for writing: $!\n";
		print FH ">$locus\n";
		print FH $qSeq{$locus}, "\n";
		close FH;
		print "Blasting $locus\n" if $verbose > 1;
		$blastRes{$locus} = "" if !defined( $blastRes{$locus} );
		foreach my $sp ( @INGROUP, @OUTGROUP ) {

			# Run blast
			$blastCmd =
			    "$BLASTALL -a $numProc -p blastp -d $sp -I -e $eCutOff "
			  . "-m 8 -i $tmpFile";
			my $blastOutput = `$blastCmd`;
			$blastRes{$locus} .= $blastOutput;
		}
		$count++;
		unlink $tmpFile;
		if ( $verbose == 1 ) {
			$| = 1;    # flush STDOUT
			print "\b" x 16;
			printf "%6.2f%% complete", $count / $nq * 100;
			$| = 0;
		}
	}
	print "\n" if $verbose == 1;

	untie %blastRes;
	untie %geneLen;

}


#
# MCL Clustering
# Edge weight is -log(e-value)*sqrt(alignment_length)
#
sub mclCluster() {
	my $mcl_exe   = "mcl";
	my %lengthOfGene;
	my $log10 = log(10);
	my %blastRes;
	my %geneLen;
	my $numProc = $NCPU;
	my $mciFile = "$OID_BLASTDIR/weights.mci";
	my $tabFile = "$OID_BLASTDIR/weights.tab";
	my $weightFile = "$OID_BLASTDIR/weights.dat.gz"; # for debugging only
	my $gz;  # gzip object
	my $thisFunc = ( caller(0) )[3];
	
	
	print "Generating edge weights for clustering ...\n";
	if ($DEBUG) {
		$gz = new IO::Compress::Gzip $weightFile or die "IO::Compress::Gzip failed :$GzipError\n";
	}
	# Convert directly into undirected graph by using maximum of two edge weights for the same edge
	# and output info mcl matrix format.
	open(MFH, "|mcxload -abc - --stream-mirror -re max -o $mciFile -write-tab $tabFile");
	
	tie %blastRes, "DB_File", $BLAST_RES_DB
	  or die "Cannot open $BLAST_RES_DB: $!\n";
	tie %geneLen, "DB_File", $GENE_LEN_DB
	  or die "Cannot open $GENE_LEN_DB: $!\n";
	foreach my $g ( keys %blastRes ) {
	
		my %seen;
		my $seqLen = $geneLen{$g};
		my @blastOutput = split( /\n/, $blastRes{$g});
		foreach (@blastOutput) {
			my ($target, $alignLen, $eVal) = (split)[1,3,10];

			# only use the first entry of each target
			next if defined $seen{$target};

			$seen{$target} = 1;
			my $targetSeqLen = $geneLen{$target};
			my $alignRatio = $alignLen / (($seqLen > $targetSeqLen) ? $seqLen : $targetSeqLen);
			my $weight;
			if ($eVal != 0) {
				$weight = -log($eVal)/$log10;
			}
			else {
				$weight = 500;
			}
			$weight *= sqrt($alignRatio);
			$gz->print("$g\t$target\t$weight\n") if $DEBUG;
			print MFH "$g\t$target\t$weight\n";
		}
		
	
	}
	untie %blastRes;
	untie %geneLen;
	close MFH;
	$gz->close() if $DEBUG;
	
	# run mcl
	print "Running mcl ...\n";
	#my @mclArgs = ($mclInput, "--abc", "-I", "1.4", "-o", $mclOutput);
	my @mclArgs = ($mciFile, "-scheme", "7", "-I", "1.4", "-use-tab", $tabFile, "-o", $clustersFile);
	push(@mclArgs, "-te", $numProc) if $numProc > 1;
	my $status = system($mcl_exe, @mclArgs);
	die "$thisFunc: $mcl_exe exit status = $?" unless $status == 0;
	print "mcl done.\n";
}

sub makeFamily() {
	my $verbose = $VERBOSE;
	my $familyNum = 1;
	my $singletNum = 0;
	my $familyDir;
	my $singletFile = "$OID_DATADIR/singlets";
	my $FASTACMD = defined($BLAST_HOME) ? "$BLAST_HOME/bin/fastacmd" : "fastacmd";
	
	# Clustering
	mclCluster if ! -f $clustersFile;
	
	# Output gene families
	open MFH, $clustersFile or die "Cannot open $clustersFile: $!\n";

	open SFH, ">>$singletFile"
	  or die "Create open $singletFile for writing: $!\n";
	
	while (<MFH>) {
		chomp;
		my $seq;
		my @out;
		my @genes = split;
		if (@genes == 1) {
			print SFH "$genes[0]\n";
			$singletNum++;
		}
		else {
			print "Writing family $familyNum (" . scalar(@genes) . " sequences)\n" if $verbose == 2;
			$familyDir = "$OID_DATADIR/$familyNum";
			mkdir $familyDir, 0755 or die "Cannot create $familyDir: $!";
			open FFH, ">$familyDir/$unalignedFamily"
			  or die "Cannot open $familyDir/$unalignedFamily for writing: $!";
			foreach my $gene (@genes) {
				my $sp = (parseOIDgeneName($gene))[0];
				if ( grep /^$sp$/, ( @INGROUP, @OUTGROUP ) ) {
					@out = `$FASTACMD -d $sp -s $gene`;
				}
				die "Failed to look up $gene in \"$sp\": $?" if @out == 0;
				shift @out;
				my $seq = join "", @out;
				print FFH ">$gene\n$seq\n";
			}
			close FFH;
			$familyNum++;
		}
	}
	
	close MFH;
	
	print $familyNum - 1 . " families\n";
	print "$singletNum singlets\n";
}


#
# Create alignment for each family (with elision).  Need algn_elide.sh.
# Optional arg: family num pattern
#
sub alignFamily {

	my $dirRE = shift;
	$dirRE = '.*' if !defined($dirRE);
	my $verbose = $VERBOSE;

	# Defaults
	my $thisFunc = ( caller(0) )[3];

	# Change to data dir
	my $savDir = getcwd;
	chdir $OID_DATADIR;

	my $count = 0;
	print "Generating alignments ...\n" if $verbose;

	# Go over each family
	foreach my $dir (<[1-9]*>) {
		next if $dir !~ /$dirRE/;
		$count++;
		chdir $dir
		  or die "$thisFunc: failed to change to family directory $dir: $!\n";
		if ( -s $alignedFamily ) {
			print "Skipping family $dir - alignment file exists\n"
			  if $verbose == 2;
			chdir $OID_DATADIR;
			next;
		}
		else {
			print "Creating alignment for family $dir\n" if $verbose == 2;
		}

		my $status;

		# Aligning
		if ( $verbose == 2 ) {
			$status = system("align_family.sh $unalignedFamily $alignedFamily");
		}
		else {
			if ( $verbose == 1 ) {
				$| = 1;    # flush STDOUT
				print "\b" x 7;
				printf "%6d", $count;
				$| = 0;
			}
			$status = system("align_family.sh $unalignedFamily $alignedFamily >/dev/null 2>&1");
		}
		warn "$thisFunc: alignment error - exit status = $?"
		  unless $status == 0;

		chdir $OID_DATADIR;
	}

	print "\n" if $verbose == 1;
	chdir $savDir;
}

#
# Create guide tree for each family.
# Optional arg: family num pattern
#
sub makeTree {

	my $dirRE = shift;
	$dirRE = '.*' if !defined($dirRE);
	my $verbose = $VERBOSE;

	# Defaults
	my $prefix    = "oid";
	my $paup_exe  = "paup";
	my $nrat      = 20;       # number of ratchets
	my $niter     = 100;      # number of iterations per ratchet
	my $timeLimit = 10800;    # Time limit for a heuristic search
	my $ratTimeLimit = 60;    # Time limit for a search within a ratchet iteration

	my $thisFunc = ( caller(0) )[3];
	my $oldFH;
	my %seq;
	my @outgroup;           # Outgroup taxa for current family

	# Change to data dir
	my $savDir = getcwd;
	chdir $OID_DATADIR;

	# Go over each family
	my $count = 0;
	print "Making trees ...\n" if $verbose;
	foreach my $dir (<[1-9]*>) {
		next if $dir !~ /$dirRE/;
		$count++;
		chdir $dir
		  or die "$thisFunc: failed to change to family directory $dir: $!\n";
		if ( $verbose == 1 ) {
			$| = 1;    # flush STDOUT
			print "\b" x 7;
			printf "%6d", $count;
			$| = 0;
		}
		elsif ( $verbose == 2 ) {
			print "Tree search for family $dir\n";
		}

		%seq      = ();
		@outgroup = ();

		# Check nexus file or tree file does not already exist
		if ( -s $prefix . ".tre" ) {
			print "Tree file already exists ... skipping\n"
			  if $verbose;
			chdir $OID_DATADIR;
			next;
		}

		open( ALIGN, "<", $alignedFamily )
		  or die "$thisFunc: failed to open alignment: $!\n";
		my $nam;
		while (<ALIGN>) {
			chomp;
			if (/^>/) {
				$nam = substr( $_, 1 );
				$seq{$nam} = "";
			}
			elsif ( defined $seq{$nam} ) {
				$seq{$nam} .= $_;
			}
		}
		close ALIGN;
		my $length = length( $seq{$nam} ); # assume they are all the same length
		my $ntaxa  = keys %seq;

		open( NEX, ">", $prefix . ".nex" )
		  or die "$thisFunc: Failed to open file for writing: $!\n";
		$oldFH = select(NEX);
		print "#NEXUS\n";

		# Data block
		print "BEGIN DATA;\n";
		print "dimensions ntax=" . $ntaxa . " nchar=" . $length . ";\n";
		print "format missing=?\n";
		print "symbols=\"ABCDEFGHIJKLMNPQRSTUVWXYZ\"\n";
		print "datatype=PROTEIN gap=-;\n\n";
		print "matrix\n";
		foreach ( keys %seq ) {
			if (length($_) < 32) {
				# Try to line up the characters
				print $_ . " " x ( 32 - length($_) );
			}
			else {
				print "$_ ";
			}
			print $seq{$_} . "\n";

			# Check if taxon is outgroup
			foreach my $og (@OUTGROUP) {
				if (/^$og#/) {
					push @outgroup, $_;
					last;
				}
			}
		}
		print ";\nend;\n";

		# If too few taxa, PAUP won't do searches.  We generate all possible
		# trees and find the shortest.
		if ( $ntaxa < 4 ) {
			print "Begin trees;\n";
			my @tmpSeq = keys %seq;
			if ( $ntaxa == 1 ) {
				print "tree PAUP_1 = [&U] ($tmpSeq[0]);\n";
			}
			elsif ( $ntaxa == 2 ) {
				print "tree PAUP_1 = [&U] ($tmpSeq[0],$tmpSeq[1]);\n";
			}
			elsif ( $ntaxa == 3 ) {
				print
				  "tree PAUP_1 = [&U] ($tmpSeq[0],$tmpSeq[1],$tmpSeq[2]);\n";
				print
				  "tree PAUP_2 = [&U] ($tmpSeq[0],($tmpSeq[1],$tmpSeq[2]));\n";
				print
				  "tree PAUP_3 = [&U] (($tmpSeq[0],$tmpSeq[1]),$tmpSeq[2]);\n";
				print
				  "tree PAUP_4 = [&U] (($tmpSeq[0],$tmpSeq[2]),$tmpSeq[1]);\n";
			}
			print "end;\n";
		}

		# PAUP block
		print "begin paup;\n";
		print "log file=" . $prefix . ".log;\n";
		print "set increase=auto notifybeep=no;\n";
		if (@outgroup) {
			print "outgroup @outgroup;\n";
			print "set root=outgroup;\n";
		}
		else {
			print "set root=midpoint;\n";
		}
		if ( $ntaxa < 4 ) {
			print STDOUT "___no search___\n" if $verbose == 2;
			print "filter best=yes;\n";
			if ( $ntaxa < 3 ) { goto SAVENQUIT }
		}
	#	elsif ( $ntaxa < 13 ) {    # use branch and bound
	#		print STDOUT "___branch and bound___\n" if $verbose == 2;
	#		print "bandb;\n";
	#	}
		elsif ( $ntaxa < 100 ) {    # use tbr swap
			print STDOUT "___tbr swap___\n" if $verbose == 2;
			print "hsearch status=no swap=tbr start=stepwise addseq=random "
			  . "nreps=20 timelimit=$timeLimit allswap=yes;\n";
		}
		else {                      # use ratchet
			print STDOUT "___ratchet___\n" if $verbose == 2;
			my $ratfrac = 0.1;
			my $ratinc  = 0.05 / $nrat;
			for ( my $i = 0 ; $i < $nrat ; $i++ ) {
				$ratfrac += $ratinc;
				my $tree_file = "$prefix.tre.$i";
				my $tmp_file  = $prefix . ".tmp";
				print "[*** ratchet $i ***]\n";
				print "cleartrees nowarn=yes;\n";
				print "hsearch status=no nrep=1 swap=tbr start=stepwise"
				  . " addseq=random timelimit=$timeLimit"
				  . " nchuck=1 chuckscore=1;\n";
				print "savetrees file=$tree_file replace;\n";
				print "savetrees file=$tmp_file replace;\n";

				foreach ( 0 .. $niter ) {
					my ( $nchosen, $weight ) = ( 0, 1 );
					my @chosen;
					$#chosen = $length - 1;
					foreach (@chosen) {
						$_ = 0;
					}
					print "[*** iteration $_ ***]\n";
					while ( $nchosen < scalar(@chosen) ) {
						print "weights $weight :";
						for ( my $n = 0 ; $n < $length ; $n++ ) {
							if ( !$chosen[$n] && rand(1) > $ratfrac ) {
								$chosen[$n] = 1;
								$nchosen++;
								print " " . ( $n + 1 );
							}
						}
						print ";\n";
						$weight++;
					}
					print "pset mstaxa=uncertain;\n";
					print "hsearch status=no start=1 swap=tbr multrees=no "
					  . "timelimit=$ratTimeLimit;\n";
					print "weights 1: 1-$length;\n";
					print "pset mstaxa=uncertain;\n";
					print "hsearch status=no start=1 swap=tbr multrees=no "
					  . "timelimit=$ratTimeLimit;\n";
					print "savetrees file=$tmp_file replace;\n";
					print "gettrees file= $tree_file mode=7;\n";
					print "savetrees file=$tree_file replace;\n";
					print "gettrees file=$tree_file mode=3 warntree=no;\n";
				}
				print "[*** ending ratchet ***]\n";
			}
			print "cleartrees nowarn=yes;\n";

			# find consensus of all independent tatchet searches
			print "gettrees file=$prefix.tre.$_ mode=7;\n"
			  foreach ( 0 .. $nrat - 1 );
			print "hsearch status=no start=current swap=tbr multrees=yes "
			  . "timelimit=$timeLimit allswap=yes nchuck=0 chuckscore=0;\n";
		}
		if (@OUTGROUP) {
			print "contree all/strict=yes treefile=$prefix.tre replace;\n";
		}
		else {
			print "contree all/strict=yes treefile=$prefix.tre replace;\n";
		}

	   # Re-read/write the output consensus tree so that the TRANSLATE command
	   # is used in the tree file -- some other programs assume its existence...
		print "gettrees file=$prefix.tre mode=3 warntree=no;\n";

		#print "roottrees;\n" if (@OUTGROUP);
	  SAVENQUIT:
		print "savetrees file=$prefix.tre replace=yes root=yes;\n";

		# Output tree statistics
		if ( $ntaxa > 3 ) {
			print
"pscores 1 /ci=yes ri=yes rc=yes hi=yes scorefile=$prefix.score replace=yes;\n";
		}

		print "log stop;\n";
		print "quit;\n";
		print "end;\n";
		close(NEX);

		select($oldFH);

		# Execute PAUP

		my $status = system("$paup_exe -n $prefix.nex 1>/dev/null 2>&1");
		die "$thisFunc: PAUP error: $?" unless $status == 0;

		# remove temp files if any

		unlink <$prefix.tre.*>;
		unlink "$prefix.tmp" if -f "$prefix.tmp";
		chdir $OID_DATADIR;
	}

	print "\n" if $verbose == 1;
	chdir $savDir;
}



#
# Find all ortholog sets
# Given a parenthetical tree, return list of ortholog groups (list of list refs).
#
sub findOrthologs($) {
	my $dirRE = shift;
	$dirRE = '.*' if !defined($dirRE);
	my $verbose  = $VERBOSE;
	my $thisFunc = ( caller(0) )[3];
	my $treeFile = "oid.tre";
	my $orthologFile = "orthologs";

	# Change to data dir
	my $savDir = getcwd;
	chdir $OID_DATADIR;

	# Go over each family
	my $count = 0;
	print "Generating ortholog groups ...\n" if $verbose;
	foreach my $dir (<[1-9]*>) {
		next if $dir !~ /$dirRE/;
		$count++;
		if ( $verbose == 1 ) {
			$| = 1;    # flush STDOUT
			print "\b" x 7;
			printf "%6d", $count;
			$| = 0;
		}
		elsif ( $verbose == 2 ) {
			print "Finding orthologs for family $dir\n";
		}

		chdir $dir
		  or die "$thisFunc: failed to change to family directory $dir: $!\n";

		if (! -f $treeFile || -f $orthologFile) {
			chdir $OID_DATADIR;
			next;
		}
		my $parenTree = (untranslateTree($treeFile))[0];
		my $phyTree =
		  new OrthologID::PhyloTree( $parenTree, [ ( @INGROUP, @OUTGROUP ) ] );
		my @orthGroups = $phyTree->orthologGroups( @INGROUP, @OUTGROUP );
		open FH, ">$orthologFile" or die "$thisFunc: failed to open ortholog file for writing: $!\n";
		foreach my $gp (@orthGroups) {
			my @g = @$gp;
			print FH "@g\n";
		}
		close FH;
		chdir $OID_DATADIR;

	}
	print "\n" if $verbose == 1;
	chdir $savDir;
}

1;

__END__


=head1 NAME

OrthologID - Collection of subroutines for the OrthologID framework

=head1 SYNOPSIS

  use OrthologID;


=head1 DESCRIPTION

OrthologID is a parsimony-based ortholog identification framework.
The OrthologID module contains all essential subroutines for
ortholog identification, including all-against-all blast, clustering,
alignment, tree searching, and ortholog extractions.  These subroutines are
intended to be used as part of a pipeline, as in the OrthologID software
package.


=head1 SEE ALSO

Chiu, J. C., Lee, E. K., Egan, M. G., Sarkar, I. N., Coruzzi, G. M., 
and Desalle, R. 2006. OrthologID: automation of genome-scale ortholog 
identification within a parsimony framework. Bioinformatics 22, 6 
(Mar. 2006), 699-707. DOI= http://dx.doi.org/10.1093/bioinformatics/btk040 

=head1 AUTHOR

Ernest K Lee, E<lt>elee@amnh.orgE<gt>

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2007-2010 by American Museum of Natural History

=cut
