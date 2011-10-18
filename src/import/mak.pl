#!/bin/perl -w
# open file, find line, replace with special text
use warnings;

print "Editing Makefile.am...\n";

# open file(s)
$old = "../Makefile.am";
$new = "tmp";
$bak = "Makefile.bak";

$name= "nicks_solution_class";

open INFILE , "<", $old or die $!;
open OUTFILE, ">", $new  or die $!;

# edit the Makefile.am to compile our newly generated file

# replace
while($line = <INFILE>)
{
    # error check
    if($line =~ /$name/)
    {
	print "MASA IMPORT CRITICAL ERROR: Solution of that name already exists!\n";
	print "$name has already been registered!\n";
	exit 1;
    }

    # looking for our special moniker
    if($line =~ /-l33t-/)
    {
	print OUTFILE "cc_sources += $name.cpp\n";
	print OUTFILE $line;
    }
    else 
    {
	print OUTFILE $line;
    }

}

# closing and cleaning up
print "Done with Makefile.am...\n";
print "Cleaning up...\n\n";

close  INFILE or die $!;
close OUTFILE or die $!;

#rename($old, $bak);
#rename($new, $old);

# nick
# 10/17/11
