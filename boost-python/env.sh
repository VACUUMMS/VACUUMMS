export CPATH=/home1/fwillmor/miniconda3/include/python3.12

# Get the python side of boost
# conda install boost
conda activate

# Get the C++ side from libraries in the module
ml boost

# Get the python and boost_python libs in the path
export LIBRARY_PATH=/home1/fwillmor/miniconda3/lib:${LIBRARY_PATH} 
export LD_LIBRARY_PATH=/home1/fwillmor/miniconda3/lib:${LD_LIBRARY_PATH} 
