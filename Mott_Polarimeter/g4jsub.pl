#/usr/bin/env perl
use strict;
use warnings;

use Getopt::Long;
use Data::Dumper;

###############################################################################
##Quick program to submit a series of simulations to 
##Auger. Please see the help menu by running --help.
##
##Author: Josh Magee, magee@jlab.org
##June 7, 2013.
###############################################################################

#declaration of global variables, arrays, and hashes
my $user = "mjmchugh";
my $original_mac;
my $jobname;
my $mac_content;
my $Nevents;
my $Njobs;

#declaration of subroutines
sub helpscreen;
sub print_footer;   #prints macfile header
sub print_xml;      #prints xml file

my ($help);
GetOptions(
 "help|h|?"           =>  \$help,
 "Nevents|events=i"   =>  \$Nevents,
 "Njobs|jobs=i"       =>  \$Njobs,
);

#die helpscreen unless $#ARGV!=0;
die helpscreen if $help;

$Njobs  = 10   unless  $Njobs;
$Nevents=10000 unless $Nevents;
$original_mac = pop @ARGV;

#construct basename
if ($original_mac =~ m/(.*)\.mac/) {
    $jobname = $1;
} else {
    $jobname = $original_mac;
    $original_mac = $original_mac . "\.mac";
}

open ORIG, "<", $original_mac or die "cant open your .mac file: $!\n";
{
  local $/;               #restrict slurp to local block
  $mac_content = <ORIG>;  #SLUUUUUURP
}
close ORIG;

foreach my $number (1..$Njobs) {
  my $basename = "$jobname\_$number";
  my $output = "macros\/$basename\.mac";
  my $xmlout = "xml\/$basename.xml";

  #create individual mapfile
  open my $fh, ">", $output or die "can't open/create $output: $!\n";
  print $fh $mac_content;
  print_footer($fh,$basename,$Nevents);
  close $fh;

  #deal with xml file
  open my $xml, ">", $xmlout or die "can't open/create $xmlout: $!\n";
  print_xml($xml,$basename);
  close $xml;

  my $callAuger = "jsub -xml $xmlout";
  system $callAuger;
} #end foreach over files

print "End of job submissions.\n";

exit;


###############################################################################
### End of main logic. Subroutines defined below.
###############################################################################


sub helpscreen {
my $helpstring = <<EOF;
Program designed to submit multiple jobs to the Auger batchfarm.
Provide this script with a map file, the number of events to
generate, and the number of files to generate, and it will do
the rest.

Calling syntax:
  g4jsub.pl [options]
Example:
  g4sub.pl sample.mac --events 20000 --jobs 5

Options include:
  --help       displays this helpful message
  --events     set number of events in each job
              default is 10k
  --jobs       number of jobs to submit
--              default is 10
NOTE: you MUST make an xml/ and macros/ folder before using.
EOF
die $helpstring if $help;
}

sub print_footer {
  my ($fh,$basename,$Nevents) = @_;

  my $seed1 = int ( rand(1e10) );
  my $seed2 = int ( rand(1e9 ) );

  my $header =
  "
#======================#
# Macro file $basename #
#======================#

# load/execute this macro
/random/setSeeds $seed1 $seed2
/run/beamOn $Nevents

exit
";

print $fh "$header\n";
return;
}

sub print_xml {
  my ($xml,$basename) = @_;

my $xmlfile =
"
<Request>
  <Email email=\"$user\@jlab.org\" request=\"false\" job=\"true\"/>
  <Project name=\"qweak\"/>
  <Track name=\"simulation\"/>
  <Name name=\"$basename\"/>
  <OS name=\"centos62\"/>
  <Command><![CDATA[
source /home/$user/.login
cd /u/home/$user/MottG4/Mott_Polarimeter
build/mott macros/$basename\.mac
  ]]></Command>

  <Memory space=\"1200\" unit=\"MB\"/>
  <TimeLimit unit=\"minutes\" time=\"4320\"/>


  <Job>
    <Stdout dest=\"/u/home/$user/MottG4/output/$basename\.out\"/>
    <Stderr dest=\"/u/home/$user/MottG4/output/$basename\.err\"/>
  </Job>

</Request>
";
print $xml "$xmlfile\n";
return;
} #end print_xml

