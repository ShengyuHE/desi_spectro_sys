cp /pscratch/sd/s/shengyu/mocks/Y1/Abacus_v4_2/altmtl0/iron/mock0/LSScats/ELG_LOPnotqso_NGC_clustering.dat.fits ~/project_rc/main/make_LSS/example_files/
cp /pscratch/sd/s/shengyu/mocks/Y1/Abacus_v4_2/altmtl0/iron/mock0/LSScats/ELG_LOPnotqso_NGC_0_clustering.ran.fits ~/project_rc/main/make_LSS/example_files/
cp /pscratch/sd/s/shengyu/mocks/Y1/Abacus_v4_2/altmtl0/iron/mock0/LSScats/ELG_LOPnotqso_NGC_1_clustering.ran.fits ~/project_rc/main/make_LSS/example_files/
cp /pscratch/sd/s/shengyu/mocks/Y1/Abacus_v4_2/altmtl0/iron/mock0/LSScats/ELG_LOPnotqso_NGC_2_clustering.ran.fits ~/project_rc/main/make_LSS/example_files/
cp /pscratch/sd/s/shengyu/mocks/Y1/Abacus_v4_2/altmtl0/iron/mock0/LSScats/ELG_LOPnotqso_NGC_3_clustering.ran.fits ~/project_rc/main/make_LSS/example_files/



source /global/common/software/desi/users/adematti/cosmodesi_environment.sh main


srun -N 1 -n 1 -c 16 -C cpu -t 04:00:00 --qos interactive --account desi python catas_FKP.py --region 'NGC' --tracer 'ELG_LOPnotqso' --targDir /global/homes/s/shengyu/project_rc/main/make_LSS/example_files --addcatas 'realistic' 'failures' --parallel 'n'
srun -N 2 -n 4 -c 16 -C cpu -t 04:00:00 --qos interactive --account desi python MPI_v1_catas_FKP.py --region 'NGC' --tracer 'ELG_LOPnotqso' --targDir /global/homes/s/shengyu/project_rc/main/make_LSS/example_files --addcatas 'realistic' 'failures' --parallel 'y'
srun -N 2 -n 4 -c 16 -C cpu -t 04:00:00 --qos interactive --account desi python MPI_v2_catas_FKP.py --region 'NGC' --tracer 'ELG_LOPnotqso' --targDir /global/homes/s/shengyu/project_rc/main/make_LSS/example_files --addcatas 'realistic' 'failures' --parallel 'y'