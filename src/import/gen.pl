#!/bin/perl -w
# generate .cpp file for src and analytical terms
use warnings;

print "Generating MMS file...\n";

$name= "nicks_solution_class";

# open file(s)
$new = "$name.cpp";
$bak = "Makefile.bak";
$varf = "variables.var";
$srcf = "source_terms.cpp";
$anaf = "analytical_solution.cpp";

open OUTFILE, ">", $new  or die $!;

# start with headers:
print OUTFILE "// -*-c++-*-\n";
print OUTFILE "//-----------------------------------------------------------------------bl-\n";
print OUTFILE "//-----------------------------------------------------------------------el-\n";
print OUTFILE "\n#include <masa_internal.h>\n";
print OUTFILE "using namespace MASA;\n\n";

# build constructor
print OUTFILE "template <typename Scalar>\n";
print OUTFILE "MASA::$name<Scalar>::$name()\n\{\n";
print OUTFILE "  this->mmsname = \"$name\";\n";
print OUTFILE "  this->dimension = 1;\n\n";

open VARFILE , "<", $varf or die $!;
while($vf = <VARFILE>)
{	    
    chomp($vf);
    @values = split(' ', $vf);
    print OUTFILE "  this->register_var(\"$values[0]\",&$values[0]);\n";
}
close  VARFILE or die $!;

print OUTFILE "\n  this->init_var();\n\n";
print OUTFILE "\} // done with constructor\n\n";

# build init_var
print OUTFILE "template <typename Scalar>\n";
print OUTFILE "MASA::$name<Scalar>::init_var()\n\{\n";
print OUTFILE "  int err = 0;\n\n";

open VARFILE , "<", $varf or die $!;
while($vf = <VARFILE>)
{	    
    chomp($vf);
    @values = split(' ', $vf);
    print OUTFILE "  err += this->set_var(\"$values[0]\",$values[1]);\n";
}
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

	# print template information
	print OUTFILE "template <typename Scalar>\n";
	$sf=~ s/eval_q_/MASA::$name<Scalar>::eval_q_/;
	
    }
    
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
open ANAFILE , "<", $anaf or die $!;
while($af = <ANAFILE>)
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
	
	# print template information
	print OUTFILE "template <typename Scalar>\n";
	$af=~ s/eval_exact_/MASA::$name<Scalar>::eval_exact_/;	
	
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
print OUTFILE "// generated using AUTOMASA\n"; 
print OUTFILE "//---------------------------------------------------------\n";

# clean up 
print "Done with $name.cpp...\n";
print "Cleaning up...\n\n";
# #rename($old, $bak);
# #rename($new, $old);

close OUTFILE or die $!;


# nick
# 10/17/11
