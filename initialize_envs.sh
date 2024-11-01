# DESI enviroment
source /global/common/software/desi/users/adematti/cosmodesi_environment.sh main

# rc_env enviroment (used for the introduction)
module load PrgEnv-gnu cray-mpich cudatoolkit craype-accel-nvidia80 python
conda activate rc_env
# export MPICH_GPU_SUPPORT_ENABLED=1 

# gpu test enviroment
module load PrgEnv-gnu cray-mpich cudatoolkit craype-accel-nvidia80 python
conda activate gpu-aware-mpi
export MPICH_GPU_SUPPORT_ENABLED=1 

#spec_sys enviroment from jiaxi
source /global/common/software/desi/users/adematti/cosmodesi_environment.sh main
conda create -n spec_sys 
source activate spec_sys

export LSSCODE=${HOME}/project_rc/jiaxi
if [ ! -e "${LSSCODE}" ]; then
    mkdir -p ${LSSCODE}
fi

cd ${LSSCODE}
git clone git@github.com:Jiaxi-Yu/LSS.git
cd ${LSSCODE}/LSS
git pull
git checkout catastrophics

source LSS_path.sh ${LSSCODE}
python setup.py develop --user