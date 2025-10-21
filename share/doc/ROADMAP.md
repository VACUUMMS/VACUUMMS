# ROADMAP

Future development will focus exclusively on the C++ interface. 

- HDF5 support:
  - The Configuration type is an obvious candidate for implementation as an HDF5 datatype. 
  - FVI is another candidate, as an alternative to writing ASCII or TIFF data. 
  - implementation can be conditional compilation (#ifdef VACUUMMS_HDF5 ...) of IO methods.

- Removal of all C includes from C++ code

- Help functions on python interface, to clarify API while developing.

- New CLI utilities (if any) should use C++ library.

- Expand on functionality of variational modeling, forming cavity networks, etc.
