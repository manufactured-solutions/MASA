#!/bin/perl -w
# open file, find line, replace with special text
use warnings;

# start prompt
print "\n";
print " * ------------------------------------------------------------ *\n";
print " * Welcome to the MASA solution import system.                  *\n";
print " * This will import a manufactured solution you have generated. *\n";
print " * ------------------------------------------------------------ *\n\n";
print " Warning!! This functionality is experimental.\n"; 
print " It is recommended you create a backup of your MASA local library before proceeding.\n";
print " Continue? (y/n)\n";
$line = <STDIN>;

# error handling
if($line =~ "y")
{
    print "Continuing...\n\n"
}
else
{
    print "Ending\n";
    exit 0;
}
# get solution class name
print " Please input the name of your new solution class:\n";
$line = <STDIN>;
chomp($line);
print "\n";
$soln     = $line;
$new_masa = "  anim.push_back(new $soln<Scalar>());\n\n";

print "Instantiating $soln in masa_core.cpp...\n";
print $new_masa;

# open file(s)
$old = "../masa_core.cpp";
$new = "tmp";
$bak = "masa_core.bak";

open INFILE , "<", $old or die $!;
open OUTFILE, ">", $new  or die $!;

# replace
while($line = <INFILE>)
{
    # error check
    if($line =~ /new $soln<Scalar>()/)
    {
	print "MASA IMPORT CRITICAL ERROR: Solution of that name already exists!\n";
	print "$soln has already been registered!\n";
	exit 1
    }

    # looking for our special moniker
    if($line =~ /-l33t-/)
    {
	print OUTFILE $new_masa;
	print OUTFILE $line;
    }
    else
    {
	print OUTFILE $line;
    }
}

print "Done with masa_core.cpp...\n";
print "Cleaning up...\n\n";

close  INFILE or die $!;
close OUTFILE or die $!;

#rename($old, $bak);
#rename($new, $old);

print "Creating class in masa_internal.h...\n";

print "Editing Makefile.am...\n";


exit 0;


# nick
# 10/11/11
