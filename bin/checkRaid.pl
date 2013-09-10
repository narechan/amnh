#!/usr/bin/perl

# A simple perl program to send alerts if problem found with physical disks
# by using OMSA 5.2 for Dell PowerEdge
# I do not provide any guarantee that this script will work. Use at own risk.

#Written: May 2008
#By: Jerod Hammerstein mods by apurva

my $controller="0";
my $emailaddress="anarechania\@amnh.org";
my $omreport="/opt/dell/srvadmin/sbin/omreport";
my $mail="/usr/local/bin/email";
my $servername="darwin";

my $email=0;
my $subject="";

# run omreport command and put into olist
`$omreport storage pdisk controller=$controller > /tmp/raidstat.out`;
open(LS, "/tmp/raidstat.out");
while (my $line = <LS>){
    chomp $line;

    if ($line =~ /Critical$/i) {
	$email=1;
	$subject="$servername Hard Drive Critical";
    }
    if ($line =~ /Failure$/i) {
	$email=1;
	$subject="$servername Hard Drive Rebuilding";
    }
    if ($line =~ /^Failure\sPredicted\s+:\sYes$/i) {
	$email=1;
	$subject="$servername Hard Drive Predicted to Fail";
    }
}
close (LS);

#If something was found email will = 1 so send email
if ($email==1) {
    system("$omreport storage pdisk controller=$controller | $mail -s '$subject' $emailaddress");

}
`rm /tmp/raidstat.out`;

open (LOG, ">/tmp/raidlastchecked");
my (undef,undef,undef,$mday,$mon,$year) = localtime;
$year = $year+1900;
$mon += 1;
if (length($mon)  == 1) {$mon = "0$mon";}
if (length($mday) == 1) {$mday = "0$mday";}
my $today = "$mon/$mday/$year";
print LOG "$today\n";
close (LOG);

