#!/bin/bash
#SBATCH --nodes=8
#SBATCH --ntasks=32
#SBATCH --time=24:00:00
#SBATCH -C cpu
#SBATCH --qos=regular
#SBATCH --account desi

# load the desi environment
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

MOCKNUM=0
MOCK_DIR=mocks/${survey}/Abacus_${mockver}
RUN_DIR=${LSSCODE}/LSS/scripts/

"""
Pk calculation \n
--catas_type ${catas[$j]} : to calculate the pk of mocks with redshift catastrophics
--calc_win y : to calculate the window function
--thetacut 0.05: to avoid the fiber assignement influence
"""

# implement the spectroscopic systematics
mock_fn=${SCRATCH}/${MOCK_DIR}/altmtl${MOCKNUM}/${specver}/mock${MOCKNUM}/LSScats
catas=('realistic' 'failures' 'slitless')

# calculating the 2-point correlation function
for tp in `seq 0 0`; do
    for j in `seq 0 1`; do
        2pcf_fn=${SCRATCH}/statistics/${tracers[$tp]}
        echo -e "\ncalculate 2pcf for catastrophics ${tracers[$tp]} ${catas[$j]} mocks in ${mock_fn}"
        srun python ${RUN_DIR}/xirunpc.py --tracer ${tracers[${tp}]} --region NGC SGC --corr_type smu --weight_type ${weight} --njack 0 --nran ${nran[${tp}]} --basedir ${mock_fn}  --outdir ${2pcf_fn} --catas_type ${catas[$j]}
        echo "${tracers[$tp]} ${catas[$j]} mocks 2pcf statistics are complete"
    done
done

# for calculating the power spectrum
for tp in `seq 0 0`; do
    for j in `seq 0 1`; do
        pk_fn=${SCRATCH}/statistics/${tracers[$tp]}
        echo -e "\ncalculate pk for catastrophics ${tracers[$tp]} ${catas[$j]} mocks with window function in ${mock_fn}"
        srun python ${RUN_DIR}/pkrun.py --tracer ${tracers[${tp}]} --region NGC SGC --weight_type ${weight} --rebinning y --nran ${nran[${tp}]} --basedir ${mock_fn} --outdir ${pk_fn} --thetacut 0.05 --catas_type ${catas[$j]} --calc_win y
        echo "${tracers[$tp]} ${catas[$j]} mocks pk statistics are complete"
    done
done

# for tp in `seq 0 0`; do
#     # thetacut: avoid the fiber assignement influence
#     mock_fn=${SCRATCH}/${MOCK_DIR}/altmtl${MOCKNUM}/${specver}/mock${MOCKNUM}/LSScats
#     pk_fn=${SCRATCH}/statistics/${tracers[$tp]}
#     echo "calculate pk for catastrophics ${tracers[$tp]} ${catas[$j]} mocks with window clustering in ${mock_fn}"
#     srun python ${RUN_DIR}/pkrun.py --tracer ${tracers[${tp}]} --region NGC SGC --weight_type ${weight} --rebinning y --nran ${nran[${tp}]} --basedir ${mock_fn} --outdir ${pk_fn}  --thetacut 0.05 --calc_win y
# done
