#!/bin/bash
for X in `seq 1 1 10`; do
    input="${WMC_DATA}/pd-ppn_qf-CDBONN-$X.ems.bz2"
    output="${OUTPUT_MC}/LUMIN/MC-ppn_qf-CDBONN-x6-$X"
    if [ -e ${input} ]; then
        echo "MC number $X..."
        if [ -e ${output}.root ];then
            echo "...was already done."
        else
            scriptname="analysis-mc-$X.sh"
            if [ -e ${scriptname} ]; then
                echo "...already in process."
            else
                echo "... PROCESSING MC QUASI-FREE $X ..."
                echo "#!/bin/bash" >> ${scriptname}
		echo "#PBS -l walltime=12:00:00" >> ${scriptname}
                echo "cd $PWD/p_$X" >> ${scriptname}
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
