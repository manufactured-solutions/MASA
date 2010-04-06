#
# Spec file for GRVY Toolkit
#
Summary:      GRVY is a toolkit library for HPC Application Development
Name:         grvy
Version:      0.10.0
Release:      1
License:      GPL
Group:        applications
Source:       grvy-%{version}.tar.gz
Distribution: Koomie Linux
Vendor:       Koomie
Packager:     karl@ices.utexas.edu

BuildRoot: /var/tmp/%{name}-%{version}-buildroot

%define _topdir /h1/karl/build/rpms/
%define APPS /h1/karl/public
%define MODULES modulefiles
%define BOOST_DIR /org/centers/pecos/LIBRARIES/BOOST/boost-1.37.0-gcc-4.3.2-ubuntu-amd64

%define INSTALL_DIR %{APPS}/%{name}/%{version}
%define MODULE_DIR  %{APPS}/%{MODULES}/%{name}

# PECOS Library settings

%define APPS /org/centers/pecos/LIBRARIES
%define MODULES /org/centers/pecos/modulefiles

%define INSTALL_DIR %{APPS}/%{name}/%{version}
%define MODULE_DIR  %{MODULES}/%{name}


%description

The HPCT Toolkit is a library to support common functionality desired
across HPC application development done in house.


%prep

rm -rf $RPM_BUILD_ROOT/%{INSTALL_DIR}
mkdir -p $RPM_BUILD_ROOT/%{INSTALL_DIR}

%setup

%build

./configure FC=ifort F77=ifort CC=gcc CXX=g++ --prefix=%{INSTALL_DIR} --with-boost=%{BOOST_DIR}
make
make DESTDIR=$RPM_BUILD_ROOT install

rm -rf $RPM_BUILD_ROOT/%{MODULE_DIR}
mkdir -p $RPM_BUILD_ROOT/%{MODULE_DIR}
cat > $RPM_BUILD_ROOT/%{MODULE_DIR}/%{version} << 'EOF'
#%Module1.0###################################################################
#
# This modulefile setups the environment for the HPCT Toolkit.
#
##############################################################################


proc ModulesHelp { } {
puts stderr "The %{name} module file defines the following environment variables:"
puts stderr "HPCT_DIR, HPCT_LIB, and HPCT_INC for the location of the "
puts stderr "HPCT distribution.\n"

puts stderr "To use the HPCT library, compile the source code with the option:"
puts stderr ""
puts stderr "\t-I\$HPCT_INC "
puts stderr "\nand add the following options to the link step: "
puts stderr ""
puts stderr "\t-L\$HPCT_LIB -lhpct"
puts stderr ""
puts stderr "\npkg-config may also be used to find compile and link options."

puts stderr "\nVersion %{version}"

}

prepend-path    LD_LIBRARY_PATH   %{INSTALL_DIR}/lib
prepend-path    PKG_CONFIG_PATH   %{INSTALL_DIR}/lib/pkgconfig
prepend-path    INCLUDE           %{INSTALL_DIR}/include

setenv GRVY_DIR %{INSTALL_DIR}
setenv GRVY_INC %{INSTALL_DIR}/include
setenv GRVY_LIB %{INSTALL_DIR}/lib

EOF

cat > $RPM_BUILD_ROOT/%{MODULE_DIR}/.version.%{version} << 'EOF'
#%Module1.0#################################################
##
## version file for GRVY
##
 
set     ModulesVersion      "%{version}"
EOF

%files
%defattr(-,karl,pecos)
%{INSTALL_DIR}
%{MODULE_DIR}

%post

%clean
rm -rf $RPM_BUILD_ROOT
