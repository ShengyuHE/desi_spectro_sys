#!/bin/bash
#SBATCH --nodes=4
#SBATCH --ntasks=16
#SBATCH --time=24:00:00
#SBATCH -C cpu
#SBATCH --exclusive
#SBATCH --qos=regular
#SBATCH --account desi

# load the environment
source /global/common/software/desi/users/adematti/cosmodesi_environment.sh main

# load the LSS standard scripts
export LSSCODE=${HOME}/project_rc/jiaxi # LSS package from jiaxi, set the catastrophics branch
cd ${LSSCODE}/LSS
PYTHONPATH=$PYTHONPATH:${LSSCODE}/LSS/py
PATH=$PATH:${LSSCODE}/LSS/bin

survey=Y1
mockver=v4_2
specver=iron
tracers=(ELG_LOPnotqso LRG QSO)
names=(ELG_LOP LRG QSO)
nran=(10 8 4)
weight=default_FKP
#_angular_pip maybe for later

MOCK_DIR=mocks/${survey}/Abacus_${mockver}
RUN_DIR=${LSSCODE}/LSS/scripts/

mock_fn=${SCRATCH}/${MOCK_DIR}/altmtl${MOCKNUM}/${specver}/mock${MOCKNUM}/LSScats
pk_fn=${SCRATCH}/ctm_data/

"""
Pk calculation \n
--option elgzcatas : to calculate the pk without 1.31<z<1.33 for ELGS \n
--catas_type ${catas[$j]} : to calculate the pk of mocks with redshift catastrophics \n
catas=('realistic' 'failures' 'slitless') \n
--calc_win y : to calculate the pk with window function or not \n
"""

# implement the spectroscopic systematics
catas=('realistic' 'failures' 'slitless')

for tp in `seq 0 0`; do
    for j in `seq 0 1`; do
        echo "calculate pk for catastrophics ${tracers[$tp]} ${catas[$j]} mocks with window clustering in ${mock_fn}"
        srun python ${RUN_DIR}/pkrun.py --tracer ${tracers[${tp}]} --region NGC SGC --weight_type ${weight} --rebinning y --nran ${nran[${tp}]} --basedir ${mock_fn} --outdir ${pk_fn} --catas_type ${catas[$j]} --option elgzcatas --calc_win y
    done
done 









