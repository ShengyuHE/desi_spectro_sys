#!/bin/bash
source /global/common/software/desi/users/adematti/cosmodesi_environment.sh main
source activate spec_sys
export LSSCODE=${HOME}/project_rc/jiaxi # LSS package from jiaxi, set the catastrophics branch

survey=Y1
specver=iron
mockver='v4_2'
mockdir=/dvs_ro/cfs/cdirs/desi/survey/catalogs/${survey}/mocks/SecondGenMocks/AbacusSummit_${mockver}

tracers=(ELG_LOPnotqso LRG QSO)
notqso=(y n n)
names=(ELG_LOP LRG QSO)
nran=(10 8 2)

# the directory of mocks
MOCK_DIR=mocks/${survey}/Abacus_${mockver}

# the suffix of redshift column without the redsihft error
remove_zerror=None

for tp in `seq 2 2`; do
    echo -e "\nGenerate DESI ${tracers[${tp}]} mocks with spectroscopic systematics and update the FKP weight"
    if [ "$tracers[${tp}]" == "ELG_LOPnotqso" ]; then
        catas="realistic failures slitless"
    else
        catas="realistic failures"
    fi
    for MOCKNUM in `seq 0 0`; do
        fn=${SCRATCH}/${MOCK_DIR}/altmtl${MOCKNUM}/${specver}/mock${MOCKNUM}/LSScats
        # get the SecondGen AbacusSummit mocks (25 realisations in total)
        if [ ! -e "${fn}" ]; then
            echo "home directory ${fn} does not exists, creating it..."
            mkdir -p ${fn}
        fi
        if [ ! -e "${fn}/${tracers[${tp}]}_full_HPmapcut.dat.fits" ]; then
            echo "${tracers[${tp}]} Abacus mocks do not exists, copying them from the official directory ..."
            cp ${mockdir}/altmtl${MOCKNUM}/mock${MOCKNUM}/LSScats/${tracers[${tp}]}_*full_HPmapcut.*.fits ${fn}
            cp ${mockdir}/altmtl${MOCKNUM}/mock${MOCKNUM}/LSScats/${tracers[${tp}]}_frac_tlobs.fits ${fn}
        fi
        echo -e "\nAll ${tracers[${tp}]} files are ready."

        # From the AbacusSummit mocks (as the real Universe without any selection) to LSS clustering mocks with spectroscopic systematics
        outputs=${fn}/${tracers[${tp}]}_clustering.dat.fits
        if [ ! -e "${outputs}" ]; then
            echo -e "Create the ${tracers[${tp}]} clustering mocks with spectroscopic systematics"
            export OMP_NUM_THREADS=16
            # can not use MPI (multiple tasks or nodes, error in writting files)
            srun -N 1 -n 1 -c 16 -C cpu -t 04:00:00 --qos interactive --account desi python ${LSSCODE}/LSS/scripts/mock_tools/mkCat_SecondGen_amtl.py --base_output ${SCRATCH}/${MOCK_DIR}/altmtl${MOCKNUM} --mockver ab_secondgen --mocknum ${MOCKNUM}  --survey ${survey} --add_gtl y --specdata Y3 --tracer ${names[${tp}]} --notqso ${notqso[${tp}]} --minr 0 --maxr ${nran[${tp}]} --fulld n --fullr n --apply_veto n --use_map_veto _HPmapcut --mkclusran y  --nz y --mkclusdat y --splitGC y --targDir ${mockdir} --outmd 'notscratch' --addcatas ${catas} --remove_zerror ${remove_zerror}
        fi
        if [ -e "${outputs}" ]; then
            echo -e "${tracers[${tp}]} clustering mocks are complete"
        fi
    
        # update the WEIGHT_FKP column in the LSS catalogues 
        srun -N 1 -n 1 -c 16 -C cpu -t 04:00:00 --qos interactive --account desi python catas_FKP.py --tracer ${tracers[${tp}]} --targDir ${fn} 
    done
done
