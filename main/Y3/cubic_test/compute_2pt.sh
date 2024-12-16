# !/bin/bash
source /global/common/software/desi/users/adematti/cosmodesi_environment.sh main

srun -n 1 -c 128 -C gpu -t 04:00:00 --gpus 4 --qos interactive --account desi python compute_cubic_2pt.py