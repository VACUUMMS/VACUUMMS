/* libraries/vacuumms_cpp/fvi.cc */

#include <vacuumms/ddx.hh>
#include <vacuumms/configuration.hh>
#include <vacuumms/cavity.hh>
#include <vacuumms/fvi.hh>
#include <vacuumms/exports.hh>

#include <vacuumms/limits.h>
#include <vacuumms/rng.h>

#ifdef BUILD_TIFF_UTILS
#include <tiffio.h>
#endif

#include <math.h>


FVIX::FVIX(Configuration c, Parameters p) : 
    c{c}, p{p} 
{
    setParameters(p);
}

FVIX::FVIX()
{
}

FVIX::FVIX(Configuration _c)
{
    c = _c;
}


void FVIX::setParameters(Parameters _p)
{
    p = _p;

    //p.getIntParam((char*)"-resolution", &resolution);
    resolution = p.getIntParam("-resolution");

    p.getFloatParam((char*)"-attenuator", &attenuator);
    p.getFloatParam((char*)"-preexponential", &preexponential);
    p.getFloatParam((char*)"-sigma", &sigma);
    p.getFloatParam((char*)"-epsilon", &epsilon);
    p.getFloatParam((char*)"-temperature", &temperature);

    if (p.getFlagParam((char*)"-usage")) printUsage();

    p.getVectorParam((char*)"-box", &c.box_x, &c.box_y, &c.box_z);
}


Parameters FVIX::getParameters()
{
    return p;
}

void FVIX::setConfiguration(Configuration _c)
{
    c = _c;
}

Configuration FVIX::getConfiguration()
{
    return c;
}

void FVIX::setDimensions(std::vector<size_t> _dimensions)
{
    dimensions = _dimensions;
}

std::vector<size_t> FVIX::getDimensions()
{
    return dimensions;
}

#ifdef BUILD_PYBIND_BINDINGS

pybind11::array_t<vacuumms_float>FVIX::getAttraction()
{
    return pybind11::array_t<vacuumms_float>(dimensions, attraction.data());
}

pybind11::array_t<vacuumms_float>FVIX::getRepulsion()
{
    return pybind11::array_t<vacuumms_float>(dimensions, repulsion.data());
}

pybind11::array_t<vacuumms_float>FVIX::getEnergy()
{
    return pybind11::array_t<vacuumms_float>(dimensions, energy.data());
}

pybind11::array_t<vacuumms_float>FVIX::getFVI()
{
    return pybind11::array_t<vacuumms_float>(dimensions, FVI.data());
}

#endif


/*
void FVIX::printResult()
{
    for (int i=0; i<resolution; i++) 
    for (int j=0; j<resolution; j++) 
    for (int k=0; k<resolution; k++)
        printf("%f\t%f\t%f\t%f\n", 
               i*c.box_x / resolution, 
               j*c.box_y / resolution, 
               k*c.box_z / resolution, 
               0); // FTW this is from the hard coded 16 version, 
//               preexponential * exp(ea->energy[i][j][k]/(-temperature * attenuator))); 
}
*/

void* FVIX::getResult()
{
    void* result;
    return result;   
}


void FVIX::printUsage()
{
    printf("usage:     gfg2fvi        -box [10.0 10.0 10.0]\n");
    printf("                          -potential [69] \n");
    printf("                          -resolution [256] \n");
    printf("                          -temperature [298.0] \n");
    printf("                          -attenuator [1.0] \n");
    printf("                          -preexponential [1.0] \n");
    printf("                          -sigma [0.0] \n");
    printf("                          -epsilon [1.0] \n");
    printf("                          -check_device \n");
}

#ifdef BUILD_PYBIND_BINDINGS
pybind11::str FVIX::__repr__()
{
    pybind11::str retval;
    retval += c.__repr__();
    retval += p.__repr__();
    return retval;
}
#endif

FVIX::~FVIX()
{
}

#ifdef BUILD_TIFF_UTILS

/********************************************************************************/
/*                                                                              */
/*  Reads an fvi format (%f\t%f\t%f\t%f\n") and generates a tif file            */
/*  note: height is z, width is y, depth is x                                   */
/*                                                                              */
/********************************************************************************/
void FVIX::generateTIFF(char* filename)
{ 
    //double _dim_x=256, _dim_y=256, _dim_z=256; // to capture command line args as double, to then convert to int 
    // to capture command line args as double, to then convert to int 
    int depth=dimensions[0], width=dimensions[1], height=dimensions[2]; 
    size_t pixel_volume = depth * width * height;
    if (pixel_volume > 4294967296) // (1024*1024*1024) 
    {
        fprintf(stderr, "provided diemsions of %d x %d x %d ");
        fprintf(stderr, "are larger than %ld\n", depth, width, height, 4294967296);
        fprintf(stderr, "and not supported by standard TIFF library.\n");
        fprintf(stderr, "Declining to generate.\n");
        return;
    }
    int alpha = 255;
    int green = 1, red = 1, blue = 1;
    int sampleperpixel = 4;
    char *image;


    // This was a fun bug. For dims of 1024, this works out to 2^32, and therefore zero as an int.
    // So malloc(0) returns a valid pointer, but can't write to the memory. 
    // Added the cast to long to fix 32-bit arithmetic issue.
    long blocksize = (long)depth * (long)width * (long)height * (long)sampleperpixel;
    image = (char*)malloc(blocksize);

    if (image==NULL) { fprintf(stderr, "Couldn't allocate memory."); exit(137); }
  
    for (int i=0; i<depth; i++)
    for (int j=0; j<width; j++)
    for (int k=0; k<height; k++)
    {
        //size_t which = idx * dim_x * dim_y + idy * dim_y + idz;
        size_t which = i * depth * width + j * height + k;
        vacuumms_float fvi = FVI[which];

        long voxel = (long)sampleperpixel * ((long)(i*width*height) + (long)(j*height) + (long)k);
//printf("voxel #%ld: %f\n", fvi);
        unsigned int fvid = floor(fvi*256);
//printf("voxel #%ld: %d\n", voxel, fvid);
//printf("writing (%d, %d, %d) = %f\n", i,j,k, fvi);

        if (red) image[0 + voxel] = fvid;
        else image[0 + voxel] = 0;
        if (green) image[1 + voxel] = fvid;
        else image[1 + voxel] = 0;
        if (blue) image[2 + voxel] = fvid;
        else image[2 + voxel] = 0;

        image[3 + voxel] = alpha;
    }

    TIFF *out = TIFFOpen(filename, "w");

    tsize_t linebytes = sampleperpixel * width;     // length in memory of one row of pixel in the image.
    unsigned char *buf = NULL;        // buffer used to store the row of pixel information for writing to file

    // We set the strip size of the file to be size of one row of pixels
    TIFFSetField(out, TIFFTAG_ROWSPERSTRIP, TIFFDefaultStripSize(out, linebytes));

        int page;
        for (page = 0; page < depth; page++)
        {
            TIFFSetField(out, TIFFTAG_IMAGEWIDTH, width);
            TIFFSetField(out, TIFFTAG_IMAGELENGTH, height);
            TIFFSetField(out, TIFFTAG_BITSPERSAMPLE, 8);
            TIFFSetField(out, TIFFTAG_SAMPLESPERPIXEL, sampleperpixel);
            TIFFSetField(out, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
            TIFFSetField(out, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_RGB);
            TIFFSetField(out, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT);

            /* We are writing single page of the multipage file */
            TIFFSetField(out, TIFFTAG_SUBFILETYPE, FILETYPE_PAGE);
            /* Set the page number */
            TIFFSetField(out, TIFFTAG_PAGENUMBER, page, depth);

            char *pbuffer;

            int row;
            for (row = 0; row < height; row++) 
            {
              pbuffer = image + linebytes * (page*height+ row);
              TIFFWriteScanline(out, pbuffer, row, 0);
            }

            TIFFWriteDirectory(out);

        } // next page

    TIFFClose(out);

    free(image);
}

#endif
