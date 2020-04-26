#!/bin/bash

PARTITION=development
while [ $# -gt 0 ] ; do
  if [ $1 = "-p" ] ; then
    shift
    PARTITION=$1
    shift
  fi
done

# first investigate just the system
HOSTNAME=`hostname -f | cut -d. -f2`
IBRUN=ibrun
if [ "${HOSTNAME}" = "longhorn" ] ; then
    PARTITION=v100
    IBRUN="mpiexec -mca common_pami_use_umr 0 -n 2"
fi

# from now on use system + MPI
HOSTNAME=`hostname -f | cut -d. -f2`-${TACC_FAMILY_MPI}
PROGRAM=packperf.${HOSTNAME}

BATCHFILE=sbatch.${HOSTNAME}
rm -f ${BATCHFILE}
touch ${BATCHFILE}

cat >>${BATCHFILE} <<EOF
#!/bin/bash

#SBATCH -J packingtest-${HOSTNAME}
#SBATCH -A A-ccsc
#SBATCH -o packingtest-${HOSTNAME}.o%j
#SBATCH -e packingtest-${HOSTNAME}.o%j
#SBATCH -n 2
#SBATCH -N 2
#SBATCH -p ${PARTITION}
#SBATCH -t 00:20:00
#SBATCH --mail-user=eijkhout@tacc.utexas.edu
#SBATCH --mail-type=end

if [ ! -f "${PROGRAM}" ] ; then 
   echo "Please build your program <<${PROGRAM}>> first"
   exit 1
fi
EOF

##
## local patches
##

if [ "${HOSTNAME}" = "frontera" ] ; then
cat >>${BATCHFILE} <<EOF
# module reset
# module use ~cazes/modulefiles/intel19
EOF
fi

if [ "${HOSTNAME}" = "ls5" ] ; then
cat >>${BATCHFILE} <<EOF
export MPICH_GNI_MAX_EAGER_MSG_SIZE=256
EOF
fi

if [ "${TACC_FAMILY_MPI}" = "spectrum_mpi" ] ; then
cat >>${BATCHFILE} <<EOF
export OMPI_MCA_common_pami_remote_umr_limit=0
ompi_info --param all all
EOF
fi

##
## execution loop
##

cat >>${BATCHFILE} <<EOF
date
echo "Job: \${SLURM_JOBID}, running on \$SLURM_NODELIST"
echo "Measuring packing performance on `hostname -f`"
echo "with loaded modules"
module list

n=8
while [ \$n -lt 150000000 ] ; do

  # exact power of 2
  echo "Sending \$n words"
  ${IBRUN} ./${PROGRAM} \$n

  # in between points
  m=\$(( 3*n/2 ))
  echo "Sending \$m words"
  ${IBRUN} ./${PROGRAM} \$m

  n=\$(( 2*n ))
done
echo "All tests on `hostname -a` completed"
EOF

echo 
echo "generated: ${BATCHFILE}"
echo
