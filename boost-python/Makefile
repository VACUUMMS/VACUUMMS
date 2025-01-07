

# The prototype is working, where I'm bringing in the includes from the original source, and linking against the libvacuumms_cpp library for the implementation. Some environment stuff needs to be set up to make all this work, but this is what it needs to use, above. 
vacuumms: vacuumms.cc
	g++ -fPIC vacuumms.cc -c -o vacuumms.o 
	g++ vacuumms.o -shared -o vacuumms.so -L/home1/fwillmor/VACUUMMS/build/libraries/ -lvacuumms_cpp -lboost_python312 -lpython3.12

clean:
	rm -rf *.o *.so

