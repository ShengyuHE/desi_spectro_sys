# DESI enviroment
source /global/common/software/desi/users/adematti/cosmodesi_environment.sh main

# rc_env enviroment
module load PrgEnv-gnu cray-mpich cudatoolkit craype-accel-nvidia80 python
conda activate rc_env
export MPICH_GPU_SUPPORT_ENABLED=1 

# gpu test enviroment
module load PrgEnv-gnu cray-mpich cudatoolkit craype-accel-nvidia80 python
conda activate gpu-aware-mpi
export MPICH_GPU_SUPPORT_ENABLED=1 

# different path
DESI_mocks_path = '/global/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/AbacusSummit_v4_2/altmtl0/mock0/LSScats/'
scrath_path = '/pscratch/sd/s/shengyu' or $SCRATCH