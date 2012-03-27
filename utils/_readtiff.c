/* readtiff.c */

#include <stdio.h>
#include <stdlib.h>

#include <tiffio.h>

main(int argc, char* argv[])
{
    TIFF* tif = TIFFOpen(argv[1], "r");
    if (tif) 
    {
        uint32 imagelength;
        tdata_t buf;
        uint32 row;
        uint16 s, nsamples;
	int pixel;

        TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &imagelength);
        TIFFGetField(tif, TIFFTAG_SAMPLESPERPIXEL, &nsamples);
printf("nsamples = %d\n", nsamples);

        buf = _TIFFmalloc(TIFFScanlineSize(tif));
        for (row = 0; row < imagelength; row++)
        {
            TIFFReadScanline(tif, buf, row, 0);
// need to get width here...
            for(pixel =0;pixel<256; pixel++) 
            {
              printf("%d", ((char)(*(char*)(pixel + buf)))/26 );
            }
printf("\n");
        }
        _TIFFfree(buf);
        TIFFClose(tif);
    }
    else printf("Could not open %s\n", argv[1]);
    exit(0);
}

#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include <ftw_makeTIFF.h>

/****************************************************************************/
/*                                                                          */
/* Takes a 1D array, and initializes it with a stack of color gradients.    */
/*                                                                          */
/* The image is sent to makeTIFF, which creates the output file.            */
/*                                                                          */
/****************************************************************************/

main()
{
  int width = 128;
  int height = 128;
  int depth = 128;

  int i,j,k;

  int sampleperpixel = 4;    // or 3 if there is no alpha channel, you should get a understanding of alpha in class soon.
  char image[width*height*depth*sampleperpixel];
  for (i=0; i<width; i++)
    for (j=0; j<height; j++)
      for (k=0; k<depth; k++)
      {
        int base= sampleperpixel * (i + j*width + k*width*height);
        image[0 + base] = i;
        image[1 + base] = j;
        image[2 + base] = k;
        image[3 + base] = 0xff;
      }
//  memset(image, 0x7F, sizeof(image));

  makeTIFF("new.tif", width, height, depth, image, sampleperpixel);
}
