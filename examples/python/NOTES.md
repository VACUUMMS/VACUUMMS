<H1> PYBIND11 Interface to VACUUMMS </H1>

<H2> How it works </H2>

In the pybind directory, there is a file vacuumms_pybind.cc which defines the interface from the python classes to the C++ classes. This interface is compiled to build the vacuumms.cpython-313-x86_64-linux-gnu.so target, which in turn depend on the VACUUMMS runtime and C++ libraries.

To use this interface, the vacuumms package must be imported into python. This ostensibly means vacuumms.cpython-313-x86_64-linux-gnu.so being in the PYTHONPATH or the current working directory. LD_LIBRARY_PATH must also include paths to vacuumms_cpp.so and vacuumms_rt.so. 

Since it is C++, all the interfaces bind to C++ implementations of the underlying functionality. E.g., The vacuumms.Parameters type is declared in vacuumms_pybind.cc with the backend declared as the VACUUMMS Parameters type.

<H2> Example </H2>

<H4> Import the interface: </H4>

    import vacuumms as v

<H4> Create an empty configuration object, which reads data in from the specified file: </H4>

    c = v.Configuration('fcc.gfg')

<H4> Set any parameters: </H4>

    p = v.Parameters(["-n", "10", "-box", "4.24264", "4.24264", "4.24264"])
    
<H4> Create a DDX (CESA algorithm) operation object to operate on the configuration object: </H4>

    o = v.DDX(c, p)

<H4> Call the execute method to execute the CESA operation </H4>

    o.execute()

<H4> Call the getResult method to get the result, in this case a CavityConfiguration object. </H4>

    cavs = o.getResult()

<H4> Call another operation (csd in this case) to get the distribution based on result </H4>

fix this

    params v.Parameters(["-n", "10", "-box", "4.24264", "4.24264", "4.24264"])
    d = v.CavitySizeDistribution(r, params)
    d.execute()
    h = d.getResult() #will be a histogram type, not yet implemented


<H4> Print the output </H4>

    dist.print()

