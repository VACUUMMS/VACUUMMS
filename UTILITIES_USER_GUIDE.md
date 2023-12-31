# Utilities User Guide

The utilities included here are typically quick text-processing utilities that perform simple manipulations such as smoothing data, sorting/histogramming, translating one format to another, etc.  These utilities are UNIX-based, are designed to be used in conjunction with the standard UNIX utilities sed, awk, and grep to build elegant pipelines for the processing of data streams.  A basic understanding of these standard UNIX utilities will enable the user to perform analyses of data from xyz, car, psf, and other common formats.

## VERSION and BUILD OPTIONS

vacuumms        Provides the release version of VACUUMMS and the options with which it was built. 

## REDUCING OPERATORS:  These operate on a list of values and return a single value

avg             calculates average value of input set
sum             calculates sum of input set
max             returns maximum value from input set
min             returns minimum value from input set
std             returns standard deviation (provide mean on command line) from input set
miss            Finds the missing number in a series of values from input set

## SCALAR OPERATORS:  These utils operate on a stream of values (one per line) to generate:

loge            the natural logs
log10           the logs base 10
exp             the exponentials
sqrt            the square roots
sq              the squares
pow             the input values raised to the power specified on the command line
recip           output the reciprocal of the input 
expr_add        the input values plus the quantity specified on the command line
expr_multiply   the input values multiplied by the quantity specified on the command line
stream_multiply like expr_multiply, but accepts multiple quantities on the command line

## TRANSFORMS:  Function or distribution in, function out

dst2hst         Transforms a distribution to a histogram
wdst2hst        Transforms a weighted distribution to a histogram

## TABULATED FUNCTION OPERATORS:

add             Adds a list of tabulated functions.
                In:  list of files (two fields, separated by tabs) on the command line
                Out: A tabulated function the sum of the input functions.
normalize       In:  Tabulated function
                Out: Normalized tabulated function
sew             Accepts a list of tabulated functions, sews them together columnwise to generate a single output

## MISCELLANEOUS:

uniq            In:  List of cavities in .cav format
                Out: List of unique cavities in .cav format
a2b             In:  2 files of .cfg format, same number of lines
                Out: a list of distances from pt A (1st file) to pt B (2nd file)
dwf             For estimating Debye-Waller factor... generates a list of squares of distances from first to second
smooth          In:  A tabulated function
                Out: Smoothed version of function
cram            In:  A .gfg format list of atoms, plus box dimensions (at command line)
                Out: A .gfg format list of atoms, wrapped for the periodic boundary condition.
povheader       Writes a povray header, according to params specified on command line
stack_tiffs     Takes a stack of 2D tiff images and converts them to a single 3D tiff
replicate_gfg   Takes a 3D configuration and copies it into the 26 surrounding periodic boxes.
replicate_cav   Takes a 3D configuration of cavities and copies it into the 26 surrounding periodic boxes.
ftw_cuda        Checks to see if a CUDA-capable device is present
truncate        Capture a subset of gfg data by specifying two corners of a bounding box.
ck              Count number of lines .cav format read
clustercat      In:  Individual clusters of cavities
                Out: Multiple clusters of cavities

## FORMAT CONVERTERS:  Operate on lists, typically of (x, y, z) positions

stream2slice    slices a data stream into files with the specified number of lines
fvi2tiff        generates a 3D TIFF from fvi data
cfg2gfg         adds columns (default values of 1.0f for sigma and epsilon)
gfg2pov         generates povray format from gfg
gfgc2pov        generates povray format from gfgc
fvi2pov         generates povray format from fvi data
cav2pov         generates povray format from cavity list
vis2pov         generates povray format from vis
cav2vis         generates vis format from cav
cfg2vis         generates vis format from cfg
gfg2fvi         A utility for generating free volume index/FVI data from configuration/gfg; requires GPU and CUDA utils
sgfg2fvi        Serial version of utility for generating free volume index/FVI data from configuration/gfg

