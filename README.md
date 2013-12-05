MASA
====

Method of Manufactured Solutions Repository


MASA (Manufactured Analytical Solution Abstraction) is a library written in C++ (with C, python and Fortran90 interfaces) which provides a suite of manufactured solutions for the software verification of partial differential equation solvers in multiple dimensions. Example formulations include:

    Heat Equation
    Laplace's Equation
    Euler Equations (with and without thermal equilibrium chemistry)
    Navier-Stokes Equations
    Reynolds Averaged Navier Stokes with Various Turbulence Models

====
Citing MASA
====

Should you use MASA in your research, please provide a citation of the library through the MASA paper, 
MASA: a library for verification using manufactured and analytical solutions published in Engineering with Computers, 
doi:10.1007/s00366-012-0267-9

Bibtex:

    @article{
    year={2012},
    issn={0177-0667},
    journal={Engineering with Computers},
    doi={10.1007/s00366-012-0267-9},
    title={MASA: a library for verification using manufactured and analytical solutions},
    url={http://dx.doi.org/10.1007/s00366-012-0267-9},
    publisher={Springer-Verlag},
    keywords={Verification; Manufactured solutions; Partial differential equations; Finite elements},
    author={Malaya, Nicholas and Estacio-Hiroms, Kemelli C. and Stogner, Roy H. and Schulz, Karl W. and Bauman, Paul T. and Carey, Graham F.},
    pages={1-10},
    language={English}
    }

For more information on the use of MASA and a brief introduction to verification, consult the following presentation:

https://red.ices.utexas.edu/attachments/download/1687/masa_intro.pdf

====
Adding Manufactured Solutions to MASA
====
MASA provides two methods to import manufactured solutions into the library. Users can either generate their own 
source terms, or they can use the automatic differentiation capabilities provided in MASA. The method by which 
solutions can be added to is provided by the "MASA-import" script. Please consult the MASA documentation for more 
details. The MASA development team will gladly add manufactured solutions to the library, for use in the 
verification community. If you would like to incorporate an MMS in the library, please contact 
masa-dev \@ ices.utexas.edu with the following information:

    Model document (written in LaTeX) detailing the equations, manufactured analytical terms
    C-code used to import the manufactured solutions into masa, using the "masa-import" feature.

Please recall that MASA is an open-source library, and all MMS are publicly available.

====
Questions/Problems?
====

We welcome any feedback regarding MASA and bugs in the code and errors or omissions in the documentation can be 
reported to masa-dev \@ ices.utexas.edu. Requests and contributions are welcome at the same e-mail address. 
Ideally, please include the following information:

    the version number of the MASA library (versioning information can be obtained by running the masa_version binary located in the bin/ directory of a local MASA installation)
    the hardware and operating system
    the local compiler version number
    a description of the bug behavior
    a short program which reproduces the bug.
    

