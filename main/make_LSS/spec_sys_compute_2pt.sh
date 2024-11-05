#!/bin/bash

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
runDir=${LSSCODE}/LSS/scripts/

# implement the spectroscopic systematics
catas=('realistic' 'failures' 'slitless')

for tp in `seq 0 0`; do
    for MOCKNUM in `seq 0 0`; do
        mock_fn=${SCRATCH}/${MOCK_DIR}/altmtl${MOCKNUM}/${specver}/mock${MOCKNUM}/LSScats
        pk_fn=${SCRATCH}/ctm_data/

        # standard 2pcf
        # outputs=${fn}/smu/xipoles_${tracers[${tp}]}_GCcomb_1.1_1.6_${weight}_lin10_njack0_nran18_split20.txt
        # if [ ! -f "${outputs}" ]; then
        #     echo "calculate clustering ${outputs}"
        #     echo srun -N 1 -C gpu -t 04:00:00 --gpus 4 --qos interactive --account desi_g python xirunpc.py --gpu --tracer ${tracers[${tp}]} --region NGC SGC --corr_type smu --weight_type ${weight} --njack 0 --nran ${nran[${tp}]} --basedir ${fn}  --outdir ${fn}
        # fi

        # 2pcf without 1.31<z<1.33
        # outputs=${fn}/smu/xipoles_${tracers[${tp}]}_GCcomb_1.1_1.6elgzcatas_${weight}_lin10_njack0_nran18_split20_thetacut0.05.txt
        # if [ ! -f "${outputs}" ]; then
        #     echo "calculate clustering ${outputs}"
        #     echo srun -N 1 -C gpu -t 04:00:00 --gpus 4 --qos interactive --account desi_g python xirunpc.py --gpu --tracer ${tracers[${tp}]} --region NGC SGC --corr_type smu --weight_type ${weight} --njack 0 --nran ${nran[${tp}]} --basedir ${fn}  --outdir ${fn} --option elgzcatas
        # fi

        # standard pk
        # outputs=${fn}/pk/pkpoles_${tracers[${tp}]}_GCcomb_1.1_1.6_${weight}_lin_thetacut0.05.npy
        # if [ ! -f "${outputs}" ]; then
        #     if [ ${MOCKNUM} == "0" ]; then
        #         echo "calculate with window clustering ${outputs}"
        #         srun -N 4 -n 16 -C cpu -t 04:00:00 --qos interactive --account desi python ${computedir}/pkrun.py --tracer ${tracers[${tp}]} --region NGC SGC --weight_type ${weight} --rebinning y --nran ${nran[${tp}]} --basedir ${mock_fn} --outdir ${pk_fn}  --calc_win y
        #     else
        #         echo "calculate without window clustering ${outputs}"
        #         srun -N 4 -n 16 -C cpu -t 04:00:00 --qos interactive --account desi python ${computedir}/pkrun.py --tracer ${tracers[${tp}]} --region NGC SGC --weight_type ${weight} --rebinning y --nran ${nran[${tp}]} --basedir ${mock_fn} --outdir ${pk_fn}        
        #     fi
        # fi
        # pk without 1.31<z<1.33
        # outputs=${fn}/pk/pkpoles_${tracers[${tp}]}_GCcomb_1.1_1.6elgzcatas_${weight}_lin_thetacut0.05.npy
        # if [ ! -f "${outputs}" ]; then
        #     if [ ${MOCKNUM} == "0" ]; then
        #         echo "calculate with window clustering ${outputs}"
        #         echo srun -N 1 -n 16 -C cpu -t 04:00:00 --qos interactive --account desi python pkrun.py --tracer ${tracers[${tp}]} --region NGC SGC --weight_type ${weight} --rebinning y --nran ${nran[${tp}]} --basedir ${fn} --outdir ${fn} --option elgzcatas  --calc_win y
        #     else
        #         echo "calculate without window clustering ${outputs}"
        #         echo srun -N 1 -n 16 -C cpu -t 04:00:00 --qos interactive --account desi python pkrun.py --tracer ${tracers[${tp}]} --region NGC SGC --weight_type ${weight} --rebinning y --nran ${nran[${tp}]} --basedir ${fn} --outdir ${fn} --option elgzcatas        
        #     fi
        # fi
        # catastrophics 2pt
        for j in `seq 0 0`; do
            # outputs=${fn}/smu/xipoles_${tracers[${tp}]}_GCcomb_1.1_1.6_${weight}_lin10_njack0_nran18_split20_${catas[$j]}_thetacut0.05.txt
            # if [ ! -f "${outputs}" ]; then 
            #     echo "calculate clustering ${outputs}"
            #     echo srun -N 1 -C gpu -t 04:00:00 --gpus 4 --qos interactive --account desi_g python xirunpc.py --gpu --tracer ${tracers[${tp}]} --region NGC SGC --corr_type smu --weight_type ${weight} --njack 0 --nran ${nran[${tp}]} --basedir ${fn}  --outdir ${fn} --catas_type ${catas[$j]}
            # fi
            outputs=${fn}/pk/pkpoles_${tracers[${tp}]}_GCcomb_1.1_1.6_${weight}_lin_${catas[$j]}_thetacut0.05.npy
            if [ ! -f "${outputs}" ]; then 
                if [ ${MOCKNUM} == "0" ]; then
                    echo "calculate with window clustering ${outputs}"
                    srun -N 1 -n 16 --exclusive -C cpu -t 04:00:00 --qos interactive --account desi python ${runDir}/pkrun.py --tracer ${tracers[${tp}]} --region NGC SGC --weight_type ${weight} --rebinning y --nran ${nran[${tp}]} --basedir ${mock_fn} --outdir ${pk_fn} --catas_type ${catas[$j]} --option elgzcatas --calc_win y
                else
                    echo "calculate without window clustering ${outputs}"
                    srun -N 1 -n 16 --exclusive -C cpu -t 04:00:00 --qos interactive --account desi python ${runDir}/pkrun.py --tracer ${tracers[${tp}]} --region NGC SGC --weight_type ${weight} --rebinning y --nran ${nran[${tp}]} --basedir ${mock_fn}  --outdir ${pk_fn} --catas_type ${catas[$j]} --option elgzcatas
                fi
            fi
        done 
    done
done
# remove the effect of the redshift uncertainties







