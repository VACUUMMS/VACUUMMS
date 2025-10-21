# Notes on building the conda package for VACUUMMS from source directory. 

VACUUMMS will build from either the modules on an HPC cluster or 
from a suitable conda environment, which means using the compiler, etc. 
in that environment. The latter is used here.

## Building VACUUMMS:

From modules:

     module load gcc cuda libtiff python/3.13 ...

From conda environment:

     mamba create -n vacuumms gcc=12 python=3.13 cuda=12 pybind11 xorg-libx11 xorg-xproto cmake numpy libtiff 

Activate the environment:

     mamba activate vacuumms

Enable plots and jupyter kernel (optional):

     mamba install matplotlib ipykernel

Then:

     mkdir VACUUMMS/build
     cd VACUUMMS/build
     export PREFIX=/usr/local/vacuumms # or wherever you want it installed.

     cmake .. -DPython3_EXECUTABLE=$(which python3) \
              -DCMAKE_INSTALL_PREFIX=${PREFIX}      \
              -DBUILD_CUDA_COMPONENTS=Y             \
              -DBUILD_TIFF_UTILS=Y                  \
              -DBUILD_PYTHON_BINDINGS=Y
     make -j4 install


## Building the conda package:

    mamba create -n build-env
    mamba activate build-env
    mamba install -c conda-forge conda-build
    cd conda-recipe
    conda build .

The build environment must include conda-build. Note that when testing install from a local channel, 
the channel must first be indexed, in order to make the package discoverable, e.g.:

    conda index ~/miniforge3/conda-bld/linux-64/
    mamba create -n install-env
    mamba activate install-env
    mamba install -c local vacuumms
 
The package itself is either a tarball vacuumms-<version>-<build_string>.tar.bz2
or a conda archive vacuumms-<version>-<build_string>.conda depending on which version
of conda-build is used.

The built conda package can be added to the appropriate channel directory (channel/linux-64/) in the VACUUMMS source.

### Some other things to consider if build doesn't go smoothly: 

    conda clean --all
   
If the build environment or install environment are corrupted, it's best to remove them and start fresh.
On some HPC systems, building script runs best on a compute node. YMMV. 

### Alternatively, installation of a pre-built standalone package can be accomplished from the github site:
 
    mamba install https://raw.githubusercontent.com/frankwillmore/VACUUMMS/scene/channel/linux-64/vacuumms-1.2.1-py313h2bc3f7f_0.conda 
    mamba install https://raw.githubusercontent.com/frankwillmore/VACUUMMS/scene/channel/linux-aarch64/vacuumms-1.3.0-py312h025b047_0.conda

## To run VACUUMMS from the python command line:

Use an conda environment in which all of the needed dependencies are installed, and the VACUUMMS python binding is present in the PYTHONPATH. Then test with:

    import vacuumms as v

## To run VACUUMMS from a Jupyter kernel:

Make sure vacuumms installed to the environment you will be using (e.g. vacuumms here):

     conda activate vacuumms

Install a kernel so that jupyter can use this environment:

     python -m ipykernel install --user --name vacuumms --display-name "VACUUMMS Jupyter kernel"
 
Then start a Jupyter server using this environment and create a notebook using this kernel. Note that it may be a necessary workaround to put 
the directory with shared object vacuumms.cpython-313-aarch64-linux-gnu.so in the PYTHONPATH for vacuumms 
to be visible/importable. The vacuumms classes need to be visible to the server.

    export PYTHONPATH=/home/frankwillmore/vacuumms/lib
    jupyter server

### Make sure the jupyter server is configured: 

    cat ~/.jupyter/jupyter_notebook_config.py 
    c.ServerApp.ip = '0.0.0.0'  # Allow connections from outside the VM
    c.ServerApp.allow_remote_access = True
    c.ServerApp.token = ''

### Connect to the server:

    Point browser to http://192.168.254.124:8888/tree (IP will differ)

