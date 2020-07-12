# ARGS1: number of cores to use per job
# ARGS2: tissue
# ARGS3: config middle name

ncore=$1
tissue=$2
myconfig=$3

mkdir -p logs

qsub -N mixfine-$tissue -l nodes=1:ppn=$ncore -v TISSUE=$tissue,NCORE=$ncore,MYCONF=$myconfig run.qsub


