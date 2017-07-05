MASA
====

[![Build Status](https://travis-ci.org/manufactured-solutions/MASA.png?branch=master)](https://travis-ci.org/manufactured-solutions/MASA) [![Join the chat at https://gitter.im/manufactured-solutions/](https://badges.gitter.im/Join%20Chat.svg)](https://gitter.im/manufactured-solutions?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)


Method of Manufactured Solutions Repository

MASA (Manufactured Analytical Solution Abstraction) is a library written in C++ (with C, python and Fortran90 interfaces) which provides a suite of manufactured solutions for the software verification of partial differential equation solvers in multiple dimensions. Example formulations include:

    Heat Equation
    Laplace's Equation
    Euler Equations (with and without thermal equilibrium chemistry)
    Navier-Stokes Equations
    Reynolds Averaged Navier Stokes with Various Turbulence Models

Summary
====
Code verification focuses on identifying failures of the code to correctly implement a desired numerical algorithm. When performing code verification, analytical solutions to mathematical equations are used to calculate error in a corresponding approximate solution.

There are many techniques common in the software engineering community that can assist code verification; for example, in design and construction of unit and module tests that exercise specific subsets of the software, regression tests that exercise fixes for previously discovered software errors, code coverage analysis, etc. In the realm of mathematical modeling, one should also exercise the software to ensure that exact solutions to the desired equations can be recovered when the code is given the corresponding inputs. This process was termed the “method of exact solutions”. However, one is often limited by the availability of analytical solutions; there will almost certainly be no non-trivial analytical solutions available when the physics and/or geometry become complex.

An important tool that has emerged over the past decade to assist in the code verification process is the Method of Manufactured Solutions (MMS). MMS, instead of relying upon the availability of an exact solution to the governing equations, specifies a solution. This artificial solution is then substituted into the equations. Naturally, there will be a residual term since the chosen function is unlikely to be an exact solution to the equations. This residual can then be added to the code as a source term; the MMS test then uses the code to solve the modified equations and checks that the chosen function is recovered. Although previous work has focused mainly on partial differential equations (PDE’s), this idea applies to a broad range of systems in mathematical physics including nonlinear equations, systems of algebraic equations, and ordinary differential equations.

MASA began as a centralized repository for the MMS generated across the [PECOS](http://pecos.utexas.edu/) Center at [ICES](ices.utexas.edu) at the [University of Texas](utexas.edu)for use with verification. Given that there appears to be no openly available, application-independent software package that provides generated MMS source terms, solutions, etc., it was decided to centralize the Center’s MMS efforts into one library to enhance reusability and consistency across the various software packages. The library is written in C++ (with C and Fortran90 interfaces) and provides a suite of manufactured solutions for the software verification of partial differential equation solvers in multiple dimensions.


Citing MASA
====

Should you use MASA in your research, please provide a citation of the library through the MASA paper, 
[MASA: a library for verification using manufactured and analytical solutions published in Engineering with Computers](http://link.springer.com/article/10.1007%2Fs00366-012-0267-9#page-1)

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

For more information on the use of MASA and a brief introduction to verification, consult the following [repository of recent MASA presentations](https://github.com/manufactured-solutions/presentations) 


Online Documentation
====

Please check out github.io pages at: http://manufactured-solutions.github.io/MASA/


Adding Manufactured Solutions to MASA
====
MASA provides two methods to import manufactured solutions into the library. Users can either generate their own 
source terms, or they can use the automatic differentiation capabilities provided in MASA. The method by which 
solutions can be added to is provided by the "MASA-import" script. Please consult the MASA documentation for more 
details. The MASA development team will gladly add manufactured solutions to the library, for use in the 
verification community. If you would like to incorporate an MMS in the library, start a github ticket with the 
following information:

    Model document (written in LaTeX) detailing the equations, manufactured analytical terms
    C-code used to import the manufactured solutions into masa, using the "masa-import" feature.

Please recall that MASA is an open-source library, and all MMS are publicly available.

Questions/Problems?
====

We welcome any feedback regarding MASA and bugs in the code and errors or omissions in the documentation can be 
reported as a ticket on github. Requests and contributions are welcome at the location. 
Ideally, please include the following information:

    the version number of the MASA library (versioning information can be obtained by running the masa_version binary located in the bin/ directory of a local MASA installation)
    the hardware and operating system
    the local compiler version number
    a description of the bug behavior
    a short program which reproduces the bug.
    

