#!/bin/ksh
#
# Script for running OrthologID pipeline
# This script makes use of Sun Grid Engine (SGE)
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
#
# Author: Ernest K Lee <elee@amnh.org>
#

# Maximum number of SGE jobs on queue (global) to limit
# the number of waiting jobs
MAXSGE=200

# Memory size of higher memory node needed for running mcl
HIMEM="8G"

OID_VERSION="1.00"
# Check environment variables and user directory

if ! env | grep -q "^OID_HOME="; then
	print -u2 "OID_HOME not defined ... exiting"
	exit 1
fi
OID_BIN=$OID_HOME/bin

if [[ $# -eq 0 ]] && ! env | grep -q "^OID_USER_DIR="; then
	print -u2 "Must specify run directory as argument or define OID_USER_DIR"
	exit 1
elif [[ $# -eq 1 ]]; then
	export OID_USER_DIR="$1"
elif [[ $# -gt 1 ]]; then
	print -u2 "Usage: $0 [ OID_USER_DIR ]"
	exit 2
fi
if [[ ! -d $OID_USER_DIR ]]; then
	print -u2 "OID_USER_DIR \"$OID_USER_DIR\" is not a directory"
	exit 1
fi

print "== Starting OrthologID v$OID_VERSION pipeline =="
print "Run directory is $OID_USER_DIR"

# Check/parse config file
CONFIG=$OID_USER_DIR/config
if [[ ! -f $CONFIG ]]; then
	print -u2 "Config file does not exist ... exiting"
	exit 1
fi
set -A INGROUP $(grep INGROUP $CONFIG | cut -f2 -d=)
set -A OUTGROUP $(grep OUTGROUP $CONFIG | cut -f2 -d=)
NCPU=$(sed -n 's/^NCPU *= *\([0-9]*\).*/\1/p' $CONFIG)
PEQ=$(sed -n 's/^PEQUEUE *= *\([0-9a-zA-Z_-]*\).*/\1/p' $CONFIG)
if [[ -z $NCPU || -z $PEQ ]]; then
	NCPU=1
fi

# Setup BLAST database
print "Setting up BLAST database ..."
if [[ ! -d $OID_USER_DIR/blastdb ]]; then
	print -u2 "Sequence directory \"$OID_USER_DIR/blastdb\" does not exist ... exiting"
	exit 1
fi
$OID_BIN/setup_blastdb.sh "${INGROUP[@]}" "${OUTGROUP[@]}"
if [[ $? -ne 0 ]]; then
	print -u2 "Unable to setup BLAST db ... exiting"
	exit 1
fi

# Create relevant directories
if [[ ! -d $OID_USER_DIR/blast ]]; then
	mkdir $OID_USER_DIR/blast
fi
if [[ ! -d $OID_USER_DIR/data ]]; then
	mkdir $OID_USER_DIR/data
fi
if [[ ! -d $OID_USER_DIR/log/sge ]]; then
	mkdir -p $OID_USER_DIR/log/sge
fi

# Generate SGE script
SGE_SCRIPT=$OID_USER_DIR/run_oid_sge.sh
cat <<EOF >$SGE_SCRIPT
#!/bin/sh
#
# SGE job script for orthologid
#

#$ -S /bin/bash
#$ -j y
#$ -e $OID_USER_DIR/log/sge
#$ -o $OID_USER_DIR/log/sge

OID_HOME=$OID_HOME
OID_USER_DIR=$OID_USER_DIR
export OID_HOME OID_USER_DIR
PATH=$OID_BIN:$PATH

cd \$OID_USER_DIR
orthologid.pl "\$@"
if [[ "\$1" == "-b" ]]; then
	touch blast/.\$2.done
fi
EOF
# End SGE script
chmod a+x $SGE_SCRIPT

# All-all BLAST
if [[ ! -s $OID_USER_DIR/blast/blastres.dbm ]]; then
	for i in "${INGROUP[@]}" "${OUTGROUP[@]}" ; do
		print "Submitting SGE job for $i-against-all BLAST..."
		if ((NCPU==1)); then
			qsub $SGE_SCRIPT -b $i
		else
			# use PE queue
			qsub -pe $PEQ $NCPU $SGE_SCRIPT -b $i
		fi
	done
	print "Waiting for BLAST jobs to finish ..."

	while sleep 300; do
		for i in "${INGROUP[@]}" "${OUTGROUP[@]}" ; do
			if [[ ! -f $OID_USER_DIR/blast/.$i.done ]]; then
				if [[ $(qstat -u $USER | wc -l) -eq 0 ]]; then
					print -u2 "Error: $i-against-all BLAST not done, but no more SGE jobs running!"
					#exit 1
				fi
				continue 2
			fi
		done
		break
	done

	print "Merging BLAST data ..."
	$OID_BIN/orthologid.pl -B
	if ! [[ -f $OID_USER_DIR/blast/blastres.dbm && -f $OID_USER_DIR/blast/genelen.dbm ]]; then
		print -u2 "Error: failed to generate BLAST results database"
		exit 1
	fi
else
		print "All-aginst-all BLAST results exist ... skipping"
fi

# Create gene families
if ! ls -d $OID_USER_DIR/data/[1-9] >/dev/null 2>&1; then
	print "Creating families ..."
	if [[ -f $OID_USER_DIR/blast/clusters ]]; then
		# Clustering done, just create family directories
		$OID_HOME/bin/orthologid.pl -f
	else
		if ((NCPU==1)); then
			# may need to change memory requirement depending on size of dataset
			qrsh -now n -l virtual_total=$HIMEM $SGE_SCRIPT -f </dev/null
		else
			qrsh -now n -pe $PEQ $NCPU -l virtual_total=$HIMEM $SGE_SCRIPT -f </dev/null
		fi
	fi
	if [[ $? -ne 0 ]]; then
		print -u2 "Family clustering failed!"
		exit 1
	fi
else
	print "Gene families already exist"
fi

# Submit alignment and tree building jobs
print "Submitting alignment and tree search jobs ..."
fnum=1
while [[ -d $OID_USER_DIR/data/$fnum ]]; do
	# No more than MAXSGE jobs running/waiting
	if [[ $(qstat -u $USER | wc -l) -lt $((MAXSGE+2)) ]]; then
		print "Submitting ^$fnum\$"
		qsub $SGE_SCRIPT -at "^$fnum\$"
		((fnum++))
	else
		sleep 60
	fi
done
print "All alignment/tree search jobs submitted.  Waiting for last job to finish."

while sleep 600; do
	for i in $OID_USER_DIR/data/[1-9]*; do
		if [[ ! -f $i/oid.tre ]]; then
			continue 2
		fi
	done
	break
done

# Extract orthologs
print "Extracting orthologs ..."
$OID_BIN/orthologid.pl -O

# Generate big matrix
print 'Generating matrix ...'
qrsh -now n -cwd -V $OID_BIN/orth2matrix.pl

