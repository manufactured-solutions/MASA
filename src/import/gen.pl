#!/bin/perl -w
# generate .cpp file for src and analytical terms
use warnings;

print "Generating MMS file...\n";

$name= "nicks_solution_class";

# open file(s)
$new = "$name.cpp";
$bak = "Makefile.bak";
$varf = "variables.var";

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

# print "Opening $srcf...\n";
# open SRCFILE , "<", $srcf or die $!;
# while($sf = <SRCFILE>)
# {	    
#     if($sf =~ /double/)
#     {

# 	if($sf =~ /eval_q_/)
# 	{
# 	    if($sf =~ /int/)
# 	    {
# 		print "Warning: MASA importer only accepts source terms with float double arguments!\n";
# 	    }
	    
# 	    if($sf =~ /void/)
# 	    {
# 		print "Warning: MASA importer only accepts source terms with float or double arguments!\n";
# 	    }

# 	    $sf=~ s/double/Scalar/g;
# 	    $sf=~ s/float/Scalar/g;
# 	    @values = split('\(', $sf);
# 	    print OUTFILE "  $values[0]";

# 	    # this is indexed at -1 because 
# 	    # we assume that the function starts with scalar
# 	    # we are counting the number of arguments in the function call
# 	    my $size = -1; $size++ while $sf =~ /Scalar/g;

# 	    # now write the number of Scalars in the source term
# 	    print OUTFILE "\(";
# 	    for ($count = 1; $count <= $size; $count++) 
# 	    {
# 		if($count eq $size)
# 		{
# 		    print OUTFILE "Scalar\);\n";
# 		}
# 		else
# 		{
# 		    print OUTFILE "Scalar,";
# 		}
# 	    }		    

# 	    if($size eq 0 )
# 	    {
# 		print OUTFILE "\);\n";
# 	    }
# 	}   
#     }
# }
# close  SRCFILE or die $!;

# # build analytical terms





# # clean up 
# print "Done with $name.cpp...\n";
# print "Cleaning up...\n\n";
# #rename($old, $bak);
# #rename($new, $old);

# close OUTFILE or die $!;

print OUTFILE "\n\n";
print OUTFILE "//---------------------------------------------------------\n";
print OUTFILE "// generated using AUTOMASA\n"; 
print OUTFILE "//---------------------------------------------------------\n";

# nick
# 10/17/11
