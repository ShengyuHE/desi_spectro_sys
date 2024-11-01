#!/bin/bash

# load the desi environment
source /global/common/software/desi/users/adematti/cosmodesi_environment.sh main

No=1
TRACER=ELG_LOPnotqso #
RAN_MOCK_NUM=10 # ELG: 10; LRG: 8 ;QSO: 4 
REGION='NGC SGC'  # NGC or SGC or NGC SGC
MOCKS_FN=${SCRATCH}/altml${No}/mock${No}/
PK_FN=${SCRATCH}/altml${No}/data${No}/
PK_RUN=${HOME}/project_rc/main/desihub/pkrun.py

# use pkrun.py to calculate the power spectrum
srun -N 4 -n 32 -C cpu -t 04:00:00 --qos interactive --account desi python ${PK_RUN} --tracer ${TRACER} --region ${REGION} --weight_type default_FKP --rebinning y --nran ${RAN_MOCK_NUM} --basedir ${MOCKS_FN} --outdir ${PK_FN} --thetacut 0.05  --calc_win y