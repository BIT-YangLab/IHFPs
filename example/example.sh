#!/bin/bash

bashrc_dir=/home/jinlong/.bashrc
init_config=/home/jinlong/VDisk1/Jinlong/2024_homologous_parcellation/1example/init_config
MATLAB=/usr/local/bin/matlab
work_dir=/home/jinlong/VDisk1/Jinlong/2024_IHFPs/example/


cat  $init_config >> $bashrc_dir
$bashrc_dir

# start procedure (matlab)
$MATLAB -nodesktop -nodisplay -nosplash -r "cd $work_dir; IHFPs_procedure;"









