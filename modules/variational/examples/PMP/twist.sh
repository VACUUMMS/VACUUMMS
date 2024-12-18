#!/bin/bash

# take the template and use it to generate and render all of the frames
for i in {0..359}
do
    number=$(printf "%03d" ${i})
    filename=$(printf "pov/PMP-${number}.pov")
    sed -e "s/camera_light/pov\/camera_light-${number}.pov/" < PMP.pov > ${filename}
    povray -W1920 -H1080 ${filename}
done
