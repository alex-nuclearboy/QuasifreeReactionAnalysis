#!/bin/bash
for X in `seq 45934 1 46884`; do
    #input="/home/wasasoft/AllRuns"
    input="${PRESEL}/run-${X}-Tr21.presel.bz2"
    output="/data7/users/khreptak/OUTPUT/DATA/LUMIN/DATA_ppn_qf_offset-run-${X}"
    if [ -e ${input} ]; then
        echo "run number $X..."
        if [ -e ${output}.root ];then
            echo "...was already done."
        else
            scriptname="analysis-run-$X.sh"
            if [ -e ${scriptname} ]; then
                echo "...already in process."
            else
                echo "... PROCESSING DATA $X ..."
                echo "#!/bin/bash" >> ${scriptname}
		echo "request_run ${X}" >> ${scriptname}
		echo >> ${scriptname}
                echo "cd $PWD" >> ${scriptname}
		echo "./main -mode raw -fin file:${input} -r ${X} -n ${output} -abort " >> ${scriptname}
                #echo "./main -mode raw -fin cluster:${input}/run_${X} -r ${X} -n ${output} -abort" >> ${scriptname}
                echo >> ${scriptname}
		#echo "release_run ${X}" >> ${scriptname}
                echo "rm -f $PWD/${scriptname}" >> ${scriptname}
                chmod u+x ${scriptname}
                qsub -q batch ${scriptname}
                sleep 2
                echo "...done."
            fi
        fi
    fi
done
