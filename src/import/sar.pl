#!/bin/perl -w
# open file, find line, replace with special text
use warnings;

$soln     = 'blargles_friend';
$new_masa = "  anim.push_back(new $soln<Scalar>());\n\n";

# open file
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

close  INFILE or die $!;
close OUTFILE or die $!;

rename($old, $bak);
rename($new, $old);

# nick
# 10/11/11
