#!/usr/bin/perl
#######
# P-Elf
# Phylogenomenclature Classifier 
# This script classifies sequences from a FASTA 
# formatted file, and classifies based on P-Gnome
# rules
#
# Command line Syntax: ./P-Elf -ruleDir=[ruleDir] -input=[fastaFile] -out=txt
#
# I.N. Sarkar
# College of Physicians & Surgeons
# Columbia University
#
# rev June 2004
#######

# get blastall location
$BLASTlocation = `whereis blastall`;
$BLASTlocation =~ s/\/blastall\n$//;

# use colored screen output for interactive mode
use Term::ANSIColor;

# version number
$versNum = "3.2b (20071206)";

# set local current directory
$curDir = `pwd`;
chomp $curDir;

# get current time
($Second, $Minute, $Hour, $Day, $Month, $Year, $WeekDay, $DayOfYear, $IsDST) = localtime (time);
if ($Month < 10) { $Month = "0" . $Month; }
if ($Day < 10) { $Day = "0" . $Day; }
if ($Minute < 10) { $Minute = "0" . $Minute;}
if ($Second < 10) { $Second = "0" . $Second;}
$Month++;
$Year += 1900;

# store time into variable
$todayDate = "$Month.$Day.$Year $Hour:$Minute:$Second";

# set interactive mode to default (1 = yes; 0 = no)
$interactive = 0;

$verboseOutput = 0;

# set default output to macClade (0), 
# text output is (1)
$outputType = 0;

# go into interactive mode when no arguments presented
if ($#ARGV < 1)
{
	# clear screen and print out title lines
	system ("clear");
	print color ("green", "bold", "underline"), "Phylogenomenclature\n", color ("reset");
	print color ("green", "bold"), "$todayDate\nP-Elf $versNum\n", color ("reset");

	# alert use that entering interactive mode
    print "Entering interactive mode....\n\n";
    print "current directory: $curDir\n\n";

	# set interactive flag
    $interactive = 1;
    
	# ask user to type ls...
	print color ("magenta"), "\n-> type 'ls' to see directories\n";
	print color ("magenta", "bold"), "Enter Rule Directory: ", color ("reset");
	
	# get directory input from user
	$ruleDir = <STDIN>;
	chomp $ruleDir;
}
else
{
	# get rule directory and fastafile name from command line
	# go through each input element from the command line and parse out the 
	# ruleDir, fasta file, and outputtype
	# ./P-Elf -ruleDir=[ruleDir] -input=[fastaFile] -out=text
	foreach $commandLineArgument (@ARGV)
	{
		# determine the ruleDir from command line
		if ($commandLineArgument =~ s/^-ruleDir=//i)
		{
			$ruleDir = $commandLineArgument;
		}
		
		# determine input fastA file from command line
		if ($commandLineArgument =~ s/^-input=//i)
		{
			$fastaFile = $commandLineArgument;
		}
		
		# set verbose output flag
		if ($commandLineArgument =~ s/^-verbose=//i)
		{
			$verboseOutput = 1;
		}
		
		# determine output type
		if ($commandLineArgument =~ s/^-out=//i)
		{
			if ($commandLineArgument =~ m/txt/i)
			{
				$outputType = 1;
			}
		}
		
	}	
	
	#$ruleDir   = $ARGV[0];
	#$fastaFile = $ARGV[1];
}

# check and be sure that the group file can be opened
while (!open (groupFile, "$ruleDir/CAOS_groupFile.txt"))
{
	# if user types exit or quit, exit P-Elf
	if ($ruleDir eq "quit" || $ruleDir eq "exit")
	{
		system ("clear");
		exit;
	}
	
	# otherwise if user types 'ls' get directories
	elsif ($ruleDir eq "ls")
	{
		system ("ls -F | grep \/");
	}

	# ask user to type ls...
	print color ("magenta"), "\n-> type 'ls' to see directories\n";
	print color ("magenta", "bold"), "Enter Rule Directory: ", color ("reset");
	
	# get directory input from user
	$ruleDir = <STDIN>;
	chomp $ruleDir;
}

# close the open file
close (groupFile);

# get & open class file
while (!open (classSeqs, $fastaFile))
{
	# if user types quite or exit, exit P-Elf
    if ($fastaFile eq "quit" || $fastaFile eq "exit")
    {
		system ("clear");
		exit;
    }
    
    # otherwise if user types 'ls' get all files
    elsif ($fastaFile eq "ls")
    {
		system ("ls | more");
    }
    
    # otherwise if user types 'lsfa' all files ending with .fa
    elsif ($fastaFile eq "lsfa")
    {
		system ("ls *.fa | more");
    }

	# ask user to type ls...
    print color ("magenta"), "\n-> type 'ls' to see all files, 'lsfa' to see .fa files\n";
    print color ("magenta", "bold"), "Enter FastA Classification File: ", color ("reset");
    
    # get file name input from user
    $fastaFile = <STDIN>;
    chomp $fastaFile;
}

# if in interactive mode, ask user for what they want to 
# name the class file, otherwise default is set to fastaFile.class.txt
if ($interactive == 1)
{
    while (!open (outClassFile, ">$outClassFileName"))
    {
		# get file name
		print color ("magenta"), "\n-> hit enter to use $fastaFile.class.txt as name\n";
		print color ("magenta", "bold"), "Enter output CLASSIFICATION-TEXT file name: ", color ("reset");
		$outClassFileName = <STDIN>;

		chomp $outClassFileName;

		# default (no file name entered)
		if ($outClassFileName eq "")
		{
			$outClassFileName = "$fastaFile.class.txt";
		}
    }
}

# command line; default (no file name entered)
else 
{
	$fastaFile = $ARGV[1];
    open (outClassFile, ">$fastaFile.class.txt");    
}

 

# if in interactive mode, ask user for what they want to 
# name the tree file, otherwise default is set to fastaFile.tree.txt
if ($interactive == 1)
{
	# initialize the allTrees variable
	$allTrees = "";
	$outMacFileName = "$fastaFile.mac.nex";
	
		#     while (!open (outMacFile, ">$outMacFileName"))
		#     {
		# # get file name
		# print color ("magenta"), "\n-> hit enter to use $fastaFile.mac.nex as name\n";
		# print color ("magenta", "bold"), "Enter output MACCLADE file name: ", color ("reset");
		# $outMacFileName = <STDIN>;
		# 
		# chomp $outMacFileName;
		# 
		# # default (no file name entered)
		# if ($outMacFileName eq "")
		# {
		# 	$outMacFileName = "$fastaFile.mac.nex";
		# }
		#     }
}

# command line; default (no file name entered)
else 
{
	$fastaFile = $ARGV[1];
    open (outMacFile, ">$fastaFile.mac.nex");    
}

# set the blast database name, based on the rule directory
$BLASTdbName   = "$ruleDir/blastDB";

# load in the group names 
do "$ruleDir/CAOS_groupFile.txt";

# load in the alignment from the CAOS_input file
open (CAOSinputFile, "$ruleDir/CAOS_input.txt");

# read in the number of taxa and the guide tree
$taxaTree = <CAOSinputFile>;
chomp $taxaTree;
($taxaCount, $guideTree) = split (/\t/, $taxaTree);

# initialize the fullClassTree as the guideTree
$fullClassTree = $guideTree;

# read in the rest of the CAOSinputFile
while (<CAOSinputFile>)
{
	chomp $_;
	
	# seperate out the taxa number, name, and sequence
	($taxaNum, $taxaName, $taxaSeq) = split (/\t/, $_);
	
	# store the sequence with the taxa name into the taxaSeqHash
	# code the first character with a 0, indicating non-CAOS classified
	$taxaSeqHash{$taxaName} = "$taxaSeq";
}

# close the alignment file
close (alignmentFile);

# load in the attribute file
open (attributeFile, "$ruleDir/CAOS_attributesFile.txt");

# pull in header line
$headerLine = <attributeFile>;

# pull in all the attributes
while (<attributeFile>)
{
	chomp $_;
	
	($node, $group, $pos, $state, $conf, $symp) = split (/\t/, $_);
	
	$state = uc($state);
Z	

	push (@{$CAlist{$node}{$group}}, "$pos\t$state");
	
}

# close the attribute file
close (attributeFile);

# Convert class file to UNIX format...
open (classFile, $fastaFile);

# if the file contains carriage returns, 
# go through and fix the file
$_ = <classFile>;
if ($_ =~ /\r/)
{
	# close the class file, and initiate a temp file
    close (classSeqs);
    open (classTemp, ">$ruleDir/ClassTemp");
    
    # convert endlines & send to temp file
    $_ =~ s/\r\n/\n/g; # PC --> UNIX
    $_ =~ s/\r/\n/g;   # MAC -> UNIX
    print (classTemp "$_");

	# destroy old fasta file, and replace with the temp file
	# then reopen the class file
    unlink ($fastaFile);
    system ("mv $ruleDir/ClassTemp $fastaFile");
    close (classFile);
    open (classSeqs, $fastaFile);
}

# read in all FASTA sequences from class file
$classCnt = 0;
while (<classSeqs>)
{
    chomp $_;
    if ($_ =~ /[>]/)
    {
		$fastaComment = $_;
		$fastaComment =~ s/^>//;
    }
    else
    {
		$fastaComment = uc($fastaComment);
    	$classSeq{$fastaComment} .= $_;
    }
}

# classify each sequence 
foreach $classSeqName (sort keys %classSeq)
{
	# create a CAOS name for the classifications
	$CAOSclassSeqName = "CAOS\_$classSeqName";

	# begin classification output to class file
	print  outClassFile "CAOS classification for: $classSeqName\n";  

	if ($interactive == 0 && $outputType == 1)
	{
		print "CAOS classification for: $classSeqName\n";
	}
	elsif ($outputType == 1)
	{
		print "Classifying: $classSeqName\n";
	}
	
	# create temporary file to send output for blast
	open (tempFile, ">$ruleDir/procSeq");
	print tempFile ">$classSeqName\n$classSeq{$classSeqName}\n\n";
	close (tempFile);

	# run BLAST and get the top aligned query sequence
	$BLASTresults = `$BLASTlocation/blastall -p blastp -d $BLASTdbName -i $ruleDir/procSeq -W 0 -F F`;
	#$BLASTresults = `$BLASTlocation/blastall -p blastn -d $BLASTdbName -i $ruleDir/procSeq -W 0 -F F`;
	

	#open (BLASTfileTMP, ">blastFileTmp.txt");
	#print BLASTfileTMP "$BLASTresults";
	#close (BLASTfileTMP);

	unlink ("$ruleDir/procSeq");

	# change the results into an array
	@BLASTresults = split (/\n/, $BLASTresults);	
	
	# go through the blast alignment results, and 
	# just pluck the query sequence and qrystart and subjectstart'
	# values from ONLY the first alignment
	for ($i = 0; $i <= $#BLASTresults; $i++)
	{
		if ($BLASTresults[$i] =~ s/^>//)
		{
			$subjectName = $BLASTresults[$i];
		}
	
		# anchor reading to the top of the alignment with the key word "Length"
		if ($BLASTresults[$i] =~ /\bLength\b/)
		{
			
				# go forward 5 lines
				$i += 5;
				
				# get the qry start position value
				#$qryadjust = substr ($BLASTresults[$i], 6,4);

				# save just the residues from the line
				$BLASTresults[$i] =~ s/\bQuery: \b//;
				$BLASTresults[$i] =~ s/[0-9]//g;
				$BLASTresults[$i] =~ s/ //g;
				
				# store it into the querySequence variable
				$querySequence = "$BLASTresults[$i]";
				
				# go forward two lines, and pull the subject start position value
				$i += 2;
				$subjectStart = substr ($BLASTresults[$i], 6,4);

				# keep reading in lines, until hit the first instance of a double 
				# blank line (end of alignment)
				while (length($BLASTresults[$i]) > 2)
				{
					# go forward two lines & pull the query line
					# adding just the residues to the querySequence variable
					$i += 2;
					if ($BLASTresults[$i] =~ m/\bQuery: \b/)
					{
						$BLASTresults[$i] =~ s/\bQuery: \b//;
						$BLASTresults[$i] =~ s/[0-9]//g;
						$BLASTresults[$i] =~ s/ //g;
						$querySequence .= $BLASTresults[$i];
					}
				}
				
				# force the exit of the for loop
				$i = $#BLASTresults + 1;
		}
	}
	
	#print "subject name: $subjectName\n";
	#print "subject sequence:\n$taxaSeqHash{$subjectName}\n\n";
	#print "querySequence:\n$querySequence\n";
	#print "subjectStart: $subjectStart\n";
	
	# seperate out the subject & query sequences into characters
	@subjectCharArray = split (//, $taxaSeqHash{$subjectName});
	@queryCharArray   = split (//, $querySequence);
	
	
	# initialize the adjustedQuerySequence
	$adjustedQuerySequence = "";
	
	# go through the subject aligned sequence
	# and adjust the querySequence accordingly
	$currPos = 1;
	$querySeqPos = 0;
	foreach $subjectChar (@subjectCharArray)
	{
		
		
		# if the subjectChar is a -, print a gap to the adjustedQuerySequence
		if ($subjectChar eq "-")
		{
			$adjustedQuerySequence .= "-";
		}
		# otherwise, if the subjectChar position is greater or equal to the 
		# current non-gapped position, add the queryCharacter
		elsif ($currPos >= $subjectStart)
		{
			# add the querySequence character to this position
			$adjustedQuerySequence .= "$queryCharArray[$querySeqPos]";
			
			# add the characters to the symbolList, 
			# if they are not -'s or ?'s
			if ($subjectChar ne "-" && $subjectChar ne "?")
			{
				$symbolList{$subjectChar} = 1;
				$symbolList{$queryCharArray[$querySeqPos]} = 1;
			}
			
			# increment the current position
			$currPos++;
			
			# also querySeqPos
			$querySeqPos++;
			
		}
		# otherwise, just add a gap to the adjustedQuerySequence
		# and increment the current position
		else
		{
			$adjustedQuerySequence .= "-";
			$currPos++;
		}
		
	}
	
	# update the query sequence with the adjusted query sequence
	$querySequence = $adjustedQuerySequence;
	
	# add to the taxaCount
	$taxaCount++;
	
	# add the classify taxa to the caosTaxaSeqHash
	$caosTaxaSeqHash{$CAOSclassSeqName} = $querySequence;
	
	# calculate the length of the sequence
	$seqLength = length ($querySequence);
	
	#print "adjustedSequence:\n$adjustedQuerySequence\n";
	#print "querySequence:\n$querySequence\n";
	
	
	
	# pad the query sequence according to the query and subject start values
	#if ($qryadjust >= 1)
    #{
    	# pad for query sequence
    #    for ($x = 1; $x <= $qryadjust; $x++)
    #    {
    #       $querySequence= '-' . $querySequence;
    #	}
    #}
    #else
    #{
    	# pad for sbjct sequence
    #    for ($x = 1; $x <= $subjectadjust; $x++)
    #    {
    #       $querySequence = '-' . $querySequence;
    #    }
    #}

	#print "querySeqName : $classSeqName\n";
	#print "querySequence: $querySequence\n";
	
	# initialize character arrays
	@sequenceChar = ();
	@seqPosCharStates = ();
	
	# split out all the character states from the sequence
	@sequenceChar = split (//, $querySequence);
	
	$seqPos = 0;
	# add all non-gap the character states to list
	foreach $charState (@sequenceChar)
	{
		if ($charState ne "-")
		{
			push (@seqPosCharStates, "$seqPos\t$charState");
		}
		$seqPos++;
	}
	
	
	# starting at base of tree (node = 0), classify sequence
	$curNode = 0;
	while (exists $CAlist{$curNode})
	{
		@grpIsect = ();
	
		#print "processing node $curNode\n";
		# go through each group, and see which group has the greater intersection
		# with the query sequence
		foreach $group (sort {$a <=> $b} keys %{$CAlist{$curNode}})
		{
			# initialize intersection (isect) and union hash and arrays		
			@isectCA = @unionCA = ();			
			%isectCA = %unionCA = ();
			
			#print "\tprocessing group: $group\n";
			
			# determine the intersection and union between the group and the 
			# query sequence character states
			foreach $CA (@{$CAlist{$curNode}{$group}}, @seqPosCharStates)
			{
				$unionCA{$CA}++ && $isectCA{$CA}++
			}
			
			# store into lists the union and intersection
			@unionCA = keys %unionCA;
			@isectCA = keys %isectCA;
			
			# store the number of intersections for each group,
			# ensuring that the minimum number stored is 0
			if ($#isectCA > 0)
			{
				push (@grpIsect, $#isectCA);
			}
			else
			{
				push (@grpIsect, 0);
			}
		}
		
		$grpCnt = 0;
		$classGroupCnt = 0;
		$otherGroupCnt = 0;
		foreach $groupCAcnt (@grpIsect)
		{
			#print "$grpCnt = $groupCAcnt\n";
			
			if ($groupCAcnt > $classGroupCnt)
			{
				$classGroup    = $grpCnt;
				$classGroupCnt = $groupCAcnt;
				$otherGoupCnt += $groupCAcnt;
			}
			else
			{
				$otherGroupCnt += $groupCAcnt;
			}
			
			$grpCnt++; 
		}
		
		$otherGroupCnt = $otherGroupCnt - $classGroupCnt;
		
		$CSIval = $classGroupCnt - $otherGroupCnt;
		
		
		#print "CSI value is: $CSIval\n";
		
		if ($CSIval > 0)
		{
			if ($verboseOutput == 1)
			{
				print outClassFile "$groupName{$curNode}{$classGroup} ($CSIval)\n";
			}
		
			$lastGroup = $groupName{$curNode}{$classGroup};
		
			if ($interactive == 0 && $outputType == 1)
			{
				print "$groupName{$curNode}{$classGroup} ($CSIval)\n";
			}
			
			if (exists $nextNode {$curNode}{$classGroup})
			{
				$curNode = $nextNode {$curNode}{$classGroup};
			}
			elsif ($classGroupCnt = 0)
			{
				last;
			}
		}
		else
		{
			last;
		}
		
	}
	
	# initialize intersection (isect) and union hash and arrays		
	@isectCA = @unionCA = ();			
	%isectCA = %unionCA = ();
	
	# print out the classification stop to the classification file
	if ($verboseOutput == 1)
	{
		print outClassFile "stop.\n\n";
	}
	elsif ($verboseOutput == 0)
	{
		print outClassFile "     ---->     \"$lastGroup\"\n\n";
	}
	elsif ($curNode == 0 && $CSIval == 0)
	{
		print outClassFile "     ---->     (Cannot classify.)\n\n";
	}
		
	# if in interactive mode, print the classification stop to the STDOUT
	if ($interactive == 0 && $outputType == 1)
	{
		print "stop.\n\n";
	}

	# reset CSIval
	$CSIval = 0;

	# insert the classSeqName into the full class tree
	$fullClassTree =~ s/$lastGroup/$lastGroup,$CAOSclassSeqName/;
	
	# initialize the singleClassTree with the guideTree
	$singleClassTree = $guideTree;
	
	# insert the classSeqName into singleClassTree
	$singleClassTree =~ s/$lastGroup/$lastGroup,$CAOSclassSeqName/;
	
	# add the tree to the compilation of allTrees
	$allTrees .= "tree * CLASSIFY_$classSeqName = [&U] $singleClassTree;\n";
	
}

# close output classification file
close (outClassFile);

# compensate for the added guideTree/CAOS character
$seqLength++;

# add 0 = guideTree and 1 = CAOS character states to symbolList
$symbolList{0} = 1;
$symbolList{1} = 1;

# initialize the outMacFileLines
$outMacFileLines = "";

# print the top NEXUS and comment lines
$outMacFileLines .= "#NEXUS\n";
$outMacFileLines .= "[\n";
$outMacFileLines .= " File generated on $todayDate by CAOS P-Gnome $versNum\n\n";
$outMacFileLines .= " Phylogenomenclature (P-Gnome) & the Characteristic Attribute Organization System (CAOS)\n";
$outMacFileLines .= "   is written by I.N. Sarkar\n";
$outMacFileLines .= "     Columbia University College of Physicians & Surgeons Department of Biomedical Informatics\n"; 
$outMacFileLines .= "     in collaboration with the American Museum of Natural History Division of Invertebrate Zoology\n";
$outMacFileLines .= "\n";
$outMacFileLines .= " CAOS is a Columbia University Patent-Pending Technology\n";
$outMacFileLines .= " International Patent Number: WO 02/064813 A2\n";
$outMacFileLines .= " United States Patent Number: PCT/US02/03540\n";
$outMacFileLines .= " United States Priority Number: 60/267,972\n";
$outMacFileLines .= "]\n\n";

# note that the first character is the guideTree versus CAOS characters
$outMacFileLines .= "[ note: First character codes for taxon type -- 0 = guideTree; 1 = CAOS classified ]\n\n";

# begin data matrix specs
$outMacFileLines .= "BEGIN DATA;\n";
$outMacFileLines .= "\tDIMENSIONS  NTAX=$taxaCount NCHAR=$seqLength;\n";
$outMacFileLines .= "\tFORMAT  SYMBOLS = \"";

# print out the symbols, as stored in the symbolList
foreach $symbol (keys %symbolList)
{
	$outMacFileLines .= "$symbol";
}

# denote ?'s and -'s as missing and gaps, respectively
$outMacFileLines .= "\" MISSING=? GAP=- ;\n";

# begin matrix
$outMacFileLines .= "MATRIX\n";

# print out the guideTree taxa
foreach $taxaName (sort keys %taxaSeqHash)
{
	# print the taxa with it's aligned sequence
	$outMacFileLines .= "$taxaName\t0$taxaSeqHash{$taxaName}\n";
}

# print out the caos classified taxa
foreach $caosTaxaName (sort keys %caosTaxaSeqHash)
{
	# print the taxa with it's aligned sequence
	$outMacFileLines .= "$caosTaxaName\t1$caosTaxaSeqHash{$caosTaxaName}\n";
}

# close out the data matrix
$outMacFileLines .= ";\n";
$outMacFileLines .= "END;\n\n";

# begin the tree block
$outMacFileLines .= "BEGIN TREES;\n";

# print out all the individual CAOS classification trees
$outMacFileLines .= "$allTrees";

# print out the full tree of CAOS classification of all the sequences
$outMacFileLines .= "tree * CLASSIFY_ALL = [&U] $fullClassTree;\n";

# close out the tree block
$outMacFileLines .= "END;\n";
	


### create the outputMacClade Nexus file
if ($interactive == 1)
{
	print outMacFile $outMacFileLines;
}

# otherwise, if the output type is set to 0, print the macClade file to stdout
elsif ($outputType == 0)
{
	print $outMacFileLines;
}
