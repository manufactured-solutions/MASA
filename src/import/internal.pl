#!/bin/perl -w
# open file, find line, replace with special text
use warnings;

print "Creating class in masa_internal.h...\n";
print "Please input the new MMS class name.\n";
print "Do not use any special characters or spaces. (underscores are acceptable)\n";

$name= "nicks_solution_class";
$varf = "variables.var";
$srcf = "srcterms.c";
$anaf = "anaterms.c";

# open file(s)
$old   = "../masa_internal.h";
$new   = "tmp";
$bak   = "masa_internal.bak";
$count = 0;

open INFILE , "<", $old or die $!;
open OUTFILE, ">", $new  or die $!;

# down to business...

# replace
while($line = <INFILE>)
{
    # error check
    if($line =~ /class $name/)
    {
	print "MASA IMPORT CRITICAL ERROR: Solution of that name already exists!\n";
	print "$name has already been registered!\n";
	exit 1;
    }

    # looking for our special moniker
    if($line =~ /-l33t-/)
    {
	# print header
	print OUTFILE "// ------------------------------------------------------\n";
	print OUTFILE "// --------------- $name \n";
	print OUTFILE "// ------------------------------------------------------\n";

	# print template and class
	print OUTFILE "template <typename Scalar>\n";
	print OUTFILE "class MASA::$name : public manufactured_solution<Scalar>\n{\n";

	# print pi
	print OUTFILE "  using manufactured_solution<Scalar>::pi;\n";
	print OUTFILE "  using manufactured_solution<Scalar>::PI;\n\n";
	
	# print list of variables
	print OUTFILE "private:\n";

	# open file, get list of variables, populate list here:
	print "Opening $varf...\n";
	open VARFILE , "<", $varf or die $!;
	while($vf = <VARFILE>)
	{	    
	    chomp($vf);
	    @values = split(' ', $vf);
	    print OUTFILE "  Scalar $values[0];\n";
	}
	close  VARFILE or die $!;
	print OUTFILE "\npublic:\n  $name();\n  int init_var();\n";

	# print list of source terms 
	print "Opening $srcf...\n";
	open SRCFILE , "<", $srcf or die $!;
	while($sf = <SRCFILE>)
	{	    
	    if($sf =~ /double/)
	    {

		if($sf =~ /eval_q_/)
		{
		    if($sf =~ /int/)
		    {
			print "Warning: MASA importer only accepts source terms with float double arguments!\n";
		    }
		    
		    if($sf =~ /void/)
		    {
			print "Warning: MASA importer only accepts source terms with float or double arguments!\n";
		    }

		    $sf=~ s/double/Scalar/g;
		    $sf=~ s/float/Scalar/g;
		    @values = split('\(', $sf);
		    print OUTFILE $values[0];

		    # this is indexed at -1 because 
		    # we assume that the function starts with scalar
		    # we are counting the number of arguments in the function call
		    my $size = -1; $size++ while $sf =~ /Scalar/g;

		    # now write the number of Scalars in the source term
		    print OUTFILE "\(";
		    for ($count = 1; $count <= $size; $count++) 
		    {
			if($count eq $size)
			{
			    print OUTFILE "Scalar\);\n";
			}
			else
			{
			    print OUTFILE "Scalar\);\n";
			}
		    }		    
		}   
	    }
	}
	close  SRCFILE or die $!;

	# print list of analytical functions
	print "Opening $anaf...\n";
	open ANAFILE , "<", $anaf or die $!;
	while($af = <ANAFILE>)
	{	    
	    if($af =~ /double/)
	    {

		if($af =~ /eval_exact_/)
		{
		    if($af =~ /int/)
		    {
			print "Warning: MASA importer only accepts source terms with float double arguments!\n";
		    }
		    
		    if($af =~ /void/)
		    {
			print "Warning: MASA importer only accepts source terms with float or double arguments!\n";
		    }

		    $af=~ s/double/Scalar/g;
		    $af=~ s/float/Scalar/g;
		    @values = split('\(', $af);
		    print OUTFILE $values[0];

		    # this is indexed at -1 because 
		    # we assume that the function starts with scalar
		    # we are counting the number of arguments in the function call
		    my $size = -1; $size++ while $af =~ /Scalar/g;

		    # now write the number of Scalars in the source term
		    print OUTFILE "\(";
		    for ($count = 1; $count <= $size; $count++) 
		    {
			if($count eq $size)
			{
			    print OUTFILE "Scalar\);\n";
			}
			else
			{
			    print OUTFILE "Scalar\);\n";
			}
		    }		    
		}   
	    }
	}
	close  ANAFILE or die $!;
	
	# leave moniker in and trip counter
	print OUTFILE "};\n\n\n";
	print OUTFILE $line;
	$count=1;
    }
    else
    {
	print OUTFILE $line;
    }

}

if($count =~ 0)
{
    print "MASA IMPORT CRITICAL ERROR: masa_internal.h corrupted!\n";
    print "Could not find registration area!\n";
    exit 1;
}


# closing and cleaning up
print "Done with masa_core.cpp...\n";
print "Cleaning up...\n\n";

close  INFILE or die $!;
close OUTFILE or die $!;

#rename($old, $bak);
#rename($new, $old);

# nick
# 10/17/11
