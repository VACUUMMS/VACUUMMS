
# 4.568648259267288 25.431351740730385 xlo xhi
# 4.568648259267288 25.431351740730385 ylo yhi
# 4.568648259267288 25.431351740730385 zlo zhi
# --> -box 20.862703481463092 20.862703481463092 20.862703481463092

cat PMP.lmps | head -1207 | tail -990 | awk '{print $5"\t"$6"\t"$7"\t"$3}' > xyz_atoms

while read -r -a array
do 
    atom_type=${array[3]}

    case ${atom_type} in

    '1'|'2'|'6'|'10'|'12'|'15'|'21'|'22')
        printf "%f\t%f\t%f\t3.500000\t0.066000\n" ${array[0]} ${array[1]} ${array[2]}
        ;;
    '3'|'4'|'5'|'7'|'8'|'9'|'11'|'13'|'14'|'16'|'17'|'18'|'19'|'20')
        printf "%f\t%f\t%f\t2.500000\t0.030000\n" ${array[0]} ${array[1]} ${array[2]}
        ;;

    esac

done < xyz_atoms |cram -box 20.862703481463092 20.862703481463092 20.862703481463092 > PMP.gfg

# Find cavities via random insertion
# pddx -box 20.862703481463092 20.862703481463092 20.862703481463092 -n_samples 1000000 -n_threads 16 < PMP.gfg > PMP.cav 

