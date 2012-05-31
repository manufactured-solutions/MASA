#!/bin/perl -w
# open file, find line, replace with special text
use warnings;
use File::Copy;

$count = 0;

# start prompt
print "\n";
print " * ------------------------------------------------------------ *\n";
print " * Welcome to the MASA solution import system.                  *\n";
print " * This will import a manufactured solution you have generated. *\n";
print " * ------------------------------------------------------------ *\n\n";
print " Warning!! This functionality is experimental.\n"; 
print " It is recommended you create a backup of your local MASA library before proceeding.\n";
print " Continue? (y/n)\n";
$line = <STDIN>;
chomp($line);
if ($line)
{
    # error handling
    if($line =~ "y")
    {
	print "Continuing...\n\n";
    }
    else
    {
	print "Ending\n";
	exit 0;
    }
}
else
{
    print "Continuing...\n\n";
}

# get solution class name
print " Please input the name of your new MMS class (default: mms_import_example):\n";
$line = <STDIN>;
chomp($line);
if ($line)
{
    print "\n";
    $soln     = $line;
}
else
{
    $soln = "mms_import_example";
}
$name     = $soln;
$new_masa = "  anim.push_back(new $soln<Scalar>());\n\n";

# get dimension
print " Please input the MMS dimension (spatial + temporal) (default: 1):\n";
$line = <STDIN>;
chomp($line);
if ($line)
{
    print "\n";
    $dimension  = $line;
}
else
{
    $dimension = 1;
}

print " Please input the source term file (examples/source_terms.cpp):\n";
$line = <STDIN>;
chomp($line);
if ($line)
{
    chomp($line);
    print "\n";
    $srcf     = $line;
}
else
{
    $srcf = "examples/source_terms.cpp";
}

print " Please input the analytical solution file (default: examples/analytical_solution.cpp):\n";
$line = <STDIN>;
chomp($line);
if ($line)
{
    chomp($line);
    print "\n";
    $anaf     = $line;
}
else
{
    $anaf = "examples/analytical_solution.cpp";
}

print " Please input the variable file:\n";
$line = <STDIN>;
chomp($line);
if ($line)
{
    chomp($line);
    print "\n";
    $varf     = $line;
}
else
{
    $varf = "examples/variables.var";
}

print "\n\nInstantiating $soln in masa_core.cpp...\n";
print $new_masa;

# open file(s)
$old = "../masa_core.cpp";
$new = "tmp";
$bak = "masa_core.cpp";

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

rename($old, "backup/$bak");
rename($new, $old);

# --------------------------------------------------------------------------
#
#                 MASA INTERNAL GENERATOR
#
# --------------------------------------------------------------------------
print "Creating class in masa_internal.h...\n";
# open file(s)
$old   = "../masa_internal.h";
$new   = "tmp.internal";
$bak   = "masa_internal.h";
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
	print OUTFILE "namespace MASA{\ntemplate <typename Scalar>\n";
	print OUTFILE "class $name : public manufactured_solution<Scalar>\n{\n";

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

	    if($values[0]=~/^\s*$/)
	    {
		# whitespace, which we ignore
	    }
	    else
	    {  
		print OUTFILE "  Scalar $values[0];\n";
	    }
	    
	} # done with while
	close  VARFILE or die $!;
	print OUTFILE "\npublic:\n  $name();\n  int init_var();\n";

	# print list of source terms 
	print "Opening $srcf...\n";
	open SRCFILE , "<", $srcf or die $!;
	while($sf = <SRCFILE>)
	{	    

	    # look for source terms or helper functions
	    if($sf =~ /eval_q_/ || $sf =~ /helper_/)
	    {

		# ignore if they are being called as functions
		if($sf =~ /=/ || $sf =~ /\/\//)
		{
		    
		}
		else
		{

		    if($sf =~ /int/)
		    {
			print "Warning: MASA importer only accepts source terms with float double arguments!\n";
			exit 1;
		    }
		    
		    if($sf =~ /void/)
		    {
			print "Warning: MASA importer only accepts source terms with float or double arguments!\n";
			exit 1;
		    }

		    $sf=~ s/double/Scalar/g;
		    $sf=~ s/float/Scalar/g;
		    @values = split('\(', $sf);

		    # CHECK if values[0] has a scalar in it:
		    if($values[0] =~ /Scalar/)
		    {
			print OUTFILE "  $values[0]";
		    }
		    else
		    {
			print OUTFILE "  Scalar $values[0]";
		    }

		    # this is indexed at -1 because 
		    # we assume that the function starts with scalar
		    # counting the number of arguments in the function call
		    my $size = -1; $size++ while $sf =~ /Scalar/g;

		    # now write the number of Scalars in the source term
		    print OUTFILE "\(";
		    for ($count = 1; $count <= $size; $count++) 
		    {
			if($count eq $size)
			{
			    # see if the function is const
			    @bracket = split('\)', $sf);
			    $bsize = scalar (@bracket);

			    if($bsize > 1)
			    {
				if($bracket[1] =~ /const/)
				{
				    print OUTFILE "Scalar\) const;\n";
				}
				else
				{
				    print OUTFILE "Scalar\);\n";
				}

			    }
			    else
			    {
				print OUTFILE "Scalar\);\n";
			    }

			}
			else
			{
			    print OUTFILE "Scalar,";
			}
		    }		    

		    if($size eq 0 )
		    {
			# see if the function is const
			@bracket = split('\)', $sf);
			$bsize = scalar (@bracket);

			if($bsize > 1)
			{
			    if($bracket[1] =~ /const/)
			    {
				print OUTFILE "\) const;\n";
			    }
			    else
			    {
				print OUTFILE "\);\n";
			    }

			}
			else
			{
			    print OUTFILE "\);\n";
			}
		    }
		    
		    
		} # done with else{
	    } # done with line matching source term
	} # done with file
	close  SRCFILE or die $!;

	# print list of analytical functions
	print "Opening $anaf...\n";
	open ANAFILE , "<", $anaf or die $!;
	while($af = <ANAFILE>)
	{	    
	    if($af =~ /eval_exact_/ || $af =~ /helper_/)
	    {

		# ignore if they are being called as functions
		# or if they are mentioned in a comment
		if($af =~ /=/ || $af =~ /\/\//)
		{
		    
		}
		else
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
		    
		    # CHECK if values[0] has a scalar in it:
		    if($values[0] =~ /Scalar/)
		    {
			print OUTFILE "  $values[0]";
		    }
		    else
		    {
			print OUTFILE "  Scalar $values[0]";
		    }
		    
		    # this is indexed at -1 because 
		    # we assume that the function starts with scalar
		    # we are counting the number of arguments in the function call
		    my $size = -1; $size++ while $af =~ /Scalar/g;

		    # now write the number of Scalars in the analytical solution
		    print OUTFILE "\(";
		    for ($count = 1; $count <= $size; $count++) 
		    {
			# this is the last argument to the function
			if($count eq $size)
			{
			    # see if the function is const
			    @bracket = split('\)', $af);
			    $bsize = scalar (@bracket);

			    # something is passed after the function ends:
			    # check if it is const
			    if($bsize > 1)
			    {
				# it is a const!
				if($bracket[1] =~ /const/)
				{
				    print OUTFILE "Scalar\) const;\n";
				}
				else # it isnt, ignore it
				{
				    print OUTFILE "Scalar\);\n";
				}
			    }
			    else # it isnt, ignore it
			    {
				print OUTFILE "Scalar\);\n";
			    }
			}
			else # print scalar and loop
			{
			    print OUTFILE "Scalar,";
			}
		    }	

		    # this is when the user has no Scalars passed as an argument
		    if($size eq 0 )
		    {
			# see if the function is const
			@bracket = split('\)', $af);
			$bsize = scalar (@bracket);

			if($bsize > 1)
			{

			    if($bracket[1] =~ /const/)
			    {
				print OUTFILE "Scalar\) const;\n";
			    }
			    else
			    {
				print OUTFILE "Scalar\);\n";
			    }
			    
			}
			else
			{
			    print OUTFILE "\);\n";
			}
			
		    }
		} # done with else {
	    } # done with analytical term
	} # done with file
	close  ANAFILE or die $!;
	
	# leave moniker in and trip counter
	print OUTFILE "};}\n\n\n";
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

close  INFILE or die $!;
close OUTFILE or die $!;

rename($old, "backup/$bak");
rename($new, $old);

# closing and cleaning up
print "Done with masa_core.cpp...\n";
print "Cleaning up...\n\n";

# --------------------------------------------------------------------------
#
#         MAKEFILE.AM EDITOR
#
# --------------------------------------------------------------------------
print "Editing Makefile.am...\n";

# open file(s)
$old = "../Makefile.am";
$new = "tmp.mak";
$bak = "Makefile.am";

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

copy($old, "backup/$bak");
rename($new, $old);

# ---------------------------------------------------------------------------
#
#         $name.cpp GENERATOR
#
# ---------------------------------------------------------------------------

print "Generating MMS file...\n";

# open file(s)
$new = "../$name.cpp";
$bak = "$name.bak";

open OUTFILE, ">", $new  or die $!;

# start with headers:
print OUTFILE "// -*-c++-*-\n";
print OUTFILE "//-----------------------------------------------------------------------bl-\n";
print OUTFILE "//-----------------------------------------------------------------------el-\n";
print OUTFILE "\n#include <masa_internal.h>\n";
print OUTFILE "\n#include <ad_masa.h>\n";

# typedefs
print OUTFILE "\ntypedef ShadowNumber<double, long double> RawScalar;\n";
print OUTFILE "\nconst unsigned int NDIM = $dimension;\n";
print OUTFILE "\ntypedef DualNumber<RawScalar, NumberArray<NDIM, RawScalar> > FirstDerivType;\n";
print OUTFILE "\ntypedef DualNumber<FirstDerivType, NumberArray<NDIM, FirstDerivType> > SecondDerivType;\n";
print OUTFILE "\ntypedef SecondDerivType ADType;\n";

#print OUTFILE "using std::cos;\n";
#print OUTFILE "using std::sin;\n";
#print OUTFILE "using std::pow;\n";
print OUTFILE "using namespace MASA;\n\n";

# build constructor
print OUTFILE "template <typename Scalar>\n";
print OUTFILE "MASA::$name<Scalar>::$name()\n\{\n";
print OUTFILE "  this->mmsname = \"$name\";\n";
print OUTFILE "  this->dimension = $dimension;\n\n";

open VARFILE , "<", $varf or die $!;
while($vf = <VARFILE>)
{	    
    chomp($vf);
    @values = split(' ', $vf);

    if($values[0]=~/^\s*$/)
    {
	# whitespace, which we ignore
    }
    else
    {  
	print OUTFILE "  this->register_var(\"$values[0]\",&$values[0]);\n";
    }
    
}
close  VARFILE or die $!;

print OUTFILE "\n  this->init_var();\n\n";
print OUTFILE "\} // done with constructor\n\n";

# build init_var
print OUTFILE "template <typename Scalar>\n";
print OUTFILE "int MASA::$name<Scalar>::init_var()\n\{\n";
print OUTFILE "  int err = 0;\n\n";

open VARFILE , "<", $varf or die $!;
while($vf = <VARFILE>)
{	    
    chomp($vf);
    @values = split(' ', $vf);


    if($values[0]=~/^\s*$/)
    {
	# whitespace, which we ignore
    }
    else
    {  	
	$arraySize = scalar (@values);

	# no default given -- assume 12
	if($arraySize < 2)
	{
	    print OUTFILE "  err += this->set_var(\"$values[0]\",12);\n";
	}
	else # use the value we were given
	{
	    print OUTFILE "  err += this->set_var(\"$values[0]\",$values[1]);\n";
	}
    }# done with else
} # done with while
close  VARFILE or die $!;

print OUTFILE "\n  return err;\n\n\} // done with init_var\n\n";

# dump list of source terms 
print OUTFILE "// ----------------------------------------\n";
print OUTFILE "// Source Terms\n";
print OUTFILE "// ----------------------------------------\n";

# open file and find source terms
print "Opening $srcf...\n";
open SRCFILE , "<", $srcf or die $!;
while($sf = <SRCFILE>)
{	    

    # find location of each source term start
    if($sf =~ /eval_q_/ || $sf =~ /helper_/)
    {

	# ignore if they are being called as functions
	# or as comment
	if($sf =~ /=/ || $sf =~ /\/\//)
	{
	    
	}
	else
	{
	    
	    if($sf =~ /int/)
	    {
		print "Warning: MASA importer only accepts source terms with float double arguments!\n";
	    }
	    
	    if($sf =~ /void/)
	    {
		print "Warning: MASA importer only accepts source terms with float or double arguments!\n";
	    }
	    
	    # print template information
	    print OUTFILE "template <typename Scalar>\n";
	    $sf=~ s/eval_q_/MASA::$name<Scalar>::eval_q_/;
	    $sf=~ s/helper_/MASA::$name<Scalar>::helper_/;
	    
	} # done with else...
    } # done with source term
    $sf=~ s/double/Scalar/g;
    $sf=~ s/float/Scalar/g;
    
    # print line after replacements
    print OUTFILE $sf;
    
}
close  SRCFILE or die $!;

# dump list of analytical terms 
print OUTFILE "\n\n// ----------------------------------------\n";
print OUTFILE "// Analytical Terms\n";
print OUTFILE "// ----------------------------------------\n";

# build analytical terms
print "Opening $anaf...\n";
open ANAFILE , "<", $anaf or die $!;
while($af = <ANAFILE>)
{	    

    # find all analytical terms
    if($af =~ /eval_exact_/ || $af =~ /helper_/)
    {

	# ignore if they are being called as functions
	if($af =~ /=/ || $af =~ /\/\//)
	{
	    
	}
	else
	{
	    if($af =~ /int/)
	    {
		print "Warning: MASA importer only accepts source terms with float or double arguments!\n";
	    }
	    
	    if($af =~ /void/)
	    {
		print "Warning: MASA importer only accepts source terms with float or double arguments!\n";
	    }
	    
	    # print template information
	    print OUTFILE "template <typename Scalar>\n";
	    $af=~ s/eval_exact_/MASA::$name<Scalar>::eval_exact_/;	
	    $af=~ s/helper_/MASA::$name<Scalar>::helper_/;	
	    
	} # done with else
    }   
    
    $af=~ s/double/Scalar/g;
    $af=~ s/float/Scalar/g;
    
    # print line after replacements
    print OUTFILE $af;

}

# dump list of analytical terms 
print OUTFILE "\n\n// ----------------------------------------\n";
print OUTFILE "// Template Instantiation(s)\n";
print OUTFILE "// ----------------------------------------\n";
print OUTFILE "\nMASA_INSTANTIATE_ALL(MASA::$name);\n\n";

print OUTFILE "\n\n";
print OUTFILE "//---------------------------------------------------------\n";
print OUTFILE "// AUTOMASA\n"; 

# print date
($sec,$min,$hour,$mday,$mon,$year,$wday,
 $yday,$isdst)=localtime(time);
printf OUTFILE "// Generated on: %4d-%02d-%02d %02d:%02d:%02d\n",
    $year+1900,$mon+1,$mday,$hour,$min,$sec;
print OUTFILE "//---------------------------------------------------------\n";

# clean up 
print "Done with $name.cpp...\n";
print "Cleaning up...\n\n";
close OUTFILE or die $!;

print "Finished with no errors.\n\n";
print "Exiting: Have a Well Verified Day.\n\n";
exit 0;

# nick
# 10/11/11
