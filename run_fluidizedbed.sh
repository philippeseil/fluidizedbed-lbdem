#! /bin/bash

#PBS -N out_fluidizedbed
#PBS -l walltime=110:00:00
#PBS -m n
#PBS -l nodes=8:ppn=8

module unload mvapich2
module load openmpi

echo $LD_LIBRARY_PATH

cd $PBS_O_WORKDIR

REAL_OUTDIR=$OUTDIR/uin-${U_IN}_ced-${COHESION}/

rm -rf $REAL_OUTDIR
mkdir -p $REAL_OUTDIR/post
mkdir -p $REAL_OUTDIR/tmp

mpirun -hostfile ${PBS_NODEFILE} \
    ./fluidizedbed \
    -N 8 --u_in ${U_IN} --u_max_lb 0.05 --dt_dem_max 1e-6 --outdir $REAL_OUTDIR \
    --rho_f 1.205 --nu_f 1.511e-5 --rho_s 2500 --radius 0.001 --vol_frac 0.4 \
    --lx 0.02 --ly 0.02 --lz 0.16 --inlet_distance 0.01 \
    --max_t 5 --write_t 0.01 --write_flowfield true --ramp_t 5 \
    --cohesion_energy_density ${COHESION} --use_smagorinsky true --c_smago 0.15

