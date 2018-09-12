#!/bin/bash
# step_list=(499 487 475 464 453 442 432 421 411 401 392 382 373 365 355 347 338 331 323 315 307 300 293 286 279 272 266 259 253 247 241 235) #230 224 219 213 208 203 198 194 189 184 180 176 171 167 163 159 155 151 148 144 141 137 134 131 127 124 121)
# gltcs_internal=(32   31     30     29     28     27     26     25     24     23     22     21     20     19     18     17     16     15     14     13     12     11     10     9      8      7      6      5      4      3     2    1)
# gltcs_str_z=(0.0000 0.0245 0.0502 0.0749 0.1008 0.1279 0.1538 0.1837 0.2123 0.2423 0.2705 0.3035 0.3347 0.3636 0.4017 0.4337 0.4714 0.5022 0.5391 0.5777 0.6184 0.6557 0.6948 0.7358 0.7788 0.8240 0.8646 0.9143 0.9591 1.0060 1.x 1.y)
step_list=(247       241        235      230      224      219      213      208      203      198      194      189      184      180      176      171      167      163      159      155      151      148      144      141      137      134      131      127   124   121)
gltcs_internal=(30   29         28       27       26       25       24       23       22       21       20       19       18       17       16       15       14       13       12       11       10       9        8        7        6        5        4        3        2        1)
gltcs_str_z=(1.0060  1.0552     1.1069  1.1520     1.2088   1.2584   1.3210   1.3759   1.4334   1.4938   1.5443   1.6104   1.6800   1.7384   1.7994   1.8797   1.9472   2.0180   2.0923   2.1703   2.2524   2.3168   2.4068   2.4775   2.5765   2.6545   2.7361   2.8506   2.9412   3.0361)
echo "#!/bin/bash" >submit_all.sh
chmod +x submit_all.sh
echo "#!/bin/bash" >login_run_all.sh
chmod +x login_run_all.sh

echo "${!step_list[@]}"
for i in "${!step_list[@]}"
do
    echo "$i"
    sed "s/#step_list#/${step_list[$i]}/g; s/#gltcs_str_z#/${gltcs_str_z[$i]}/g; s/#gltcs_internal#/${gltcs_internal[$i]}/g" <param.template >${step_list[$i]}.param
    #    sed "s/#step#/${step_list[$i]}/g" <run.template >run_${step_list[$i]}.sh
    #chmod +x run_${step_list[$i]}.sh
    echo "qsub -n 1 -t 60 -o params/gltcs_library_high_z/logs/${step_list[$i]}.out -e params/gltcs_library_high_z/logs/${step_list[$i]}.err   --debuglog=params/gltcs_library_high_z/logs/${step_list[$i]}.colbalt  ./match_up3.py params/gltcs_library_high_z/${step_list[$i]}.param ">>submit_all.sh
    echo "./match_up3.py params/gltcs_library_high_z/${step_list[$i]}.param ">>login_run_all.sh
done
