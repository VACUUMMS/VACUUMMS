for i in {0..359}
#for i in {0..1}
do
    number=$(printf "%03d" ${i})
    filename=$(printf "pov/lj-${number}.pov")
    sed -e "s/camera_light/pov\/camera_light-${number}.pov/" < lj.pov > ${filename}
    povray -W1920 -H1080 ${filename}
done
