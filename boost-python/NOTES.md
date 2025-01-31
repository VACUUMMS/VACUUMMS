<H1> BOOST::PYTHON Interface to VACUUMMS </H1>

<H2> How it works </H2>

In boost-python directory, there is a file bd_def.cc which defines the interface from the python classes to the C++ classes. This interface is compiled to build the vacuumms.so target, which in turn depend on VACUUMMS runtime and C++ libraries.

To use this interface, the vacuumms package must be imported into python. This ostensibly means vacuumms.so being in the PYTHONPATH or potentially in the current working directory. 

Since it is C++, all the interfaces bind to C++ implementations of the underlying functionality. E.g., The vacuumms.Parameters type is declared in bp_def.cpp with the backend declared as the VACUUMMS Parameters type.

<H2> Example </H2>

<H4> Import the interface </H4>

    import vacuumms

<H4> Create an empty configuration object, which reads data in from the specified file </H4>

    c = vacuumms.Configuration('filename')

<H4> Set any parameters </H4>

    p = vacuumms.Parameters(["-N", "100", "-box", "3", "4", "5"])
    
<H4> Create a DDX (CESA algorithm) object to operate on the configuration object </H4>

    o = vacuumms.DDX(c, p)

<H4> Call the execute method to execute the CESA operation </H4>

    o.execute()

<H4> Call the getResult method to get the result, in this case a CavityConfiguration object. </H4>

    cavs = o.getResult()

<H4> Call another operation (csd in this case) to get the distribution based on result </H4>

    dist = vacuumms.csd(cavs)

<H4> Print the output </H4>

    dist.print()

