# FCC CLI example

The example configuration contains 108 atoms. There are 4 atoms in a unit cell, so 108/4 = 27
Since 3^d = 27, this is unit cell repeated three deep in each direction.
The atoms here are unit diameter, so unit cell dimension for FCC is sqrt(2) and box dim is 3*sqrt(2) = 4.24264

## Generate samples

    ddx                                 \
        -n 10000                        \
        -box 4.24264 4.24264 4.24264    \
        -n_steps 100000                 \
        -precision_parameter 0.00001    \
        -verlet_cutoff 5                \
        < fcc.gfg                       \
        > fcc.cav

## Scrub for uniqueness

    uniq                                \
        -box  4.24264 4.24264 4.24264   \
        < fcc.cav                       \
        > fcc.unq

## Generate number distribution

    awk '{print $4}'                    \
        < fcc.unq                       \
        | dst2hst                       \
            -width 0.0125               \
            -n_bins 80                  \
            > fcc.hst 

## Generate volume-weighted distribution

    awk '{print $4"\t"$4 * $4 * $4}'    \
        < fcc.unq                       \
        | wdst2hst                      \
        -width 0.0125                   \
        -n_bins 80                      \
        | normalize                     \
        > fcc-vol.hst

## Generate FVI information

    gfg2fvi                             \
        -box 4.24264 4.24264 4.24264    \
        -resolution 256                 \
        -potential 612                  \
        < fcc.gfg                       \
        > fcc.fvi

## Use FVI output to generate TIFF 

     fvi2tiff                           \
        -box 4.24264 4.24264 4.24264    \
        -dims 256 256 256               \
        -alpha 255                      \
        -o fcc.tif                      \
        < fcc.fvi

