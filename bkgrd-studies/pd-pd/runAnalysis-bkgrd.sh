#!/bin/bash
for X in `seq 1 1 5`; do
    input="/data10/users/khreptak/SIMULATION/WMC/pd-pd-$X.ems.bz2"
    output="${OUTPUT_MC}/LUMIN/MC-pd-momcut-$X"
    if [ -e ${input} ]; then
        echo "MC number $X..."
        if [ -e ${output}.root ];then
            echo "...was already done."
        else
            scriptname="analysis-mc-$X.sh"
            if [ -e ${scriptname} ]; then
                echo "...already in process."
            else
                echo "... PROCESSING MC $X ..."
                echo "#!/bin/bash" >> ${scriptname}
                echo "#PBS -l walltime=12:00:00" >> ${scriptname}
                echo "cd $PWD" >> ${scriptname}
                echo "./main -mode mc -fin file:${input} -n ${output} -abort" >> ${scriptname}
                echo >> ${scriptname}
                echo "rm -f $PWD/${scriptname}" >> ${scriptname}
                chmod u+x ${scriptname}
                qsub -q batch ${scriptname}
                echo "... run"
                sleep 2
            fi
        fi
    fi
done
