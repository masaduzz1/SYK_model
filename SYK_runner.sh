#!/bin/bash

N=6
instance_id=3
mybeta=(0.1 0.2 0.3 0.5 0.8 1.0)
delt=1.5
nstep=6

#for N in {2..11}
for test in ${mybeta[@]}
do
	
    echo "$test"
    papermill betaOTOC_ED_original_def.ipynb temp_beta${test}.ipynb -y "
    N: ${N}
    instance_id: ${instance_id}
    mybeta: ${test}
    delt: ${delt}
    nstep: ${nstep}
    "
done


dt=1.5
#nshots=20000
RR=2000
myid=(3 4 7)

#for N in {2..11}
for N in ${myid[@]}
do
	
    echo "$N"
    papermill global_protocol_ckt_efficient.ipynb temps${N}.ipynb -y "
    instance_id: ${N}
    myseed: $((${N}+233423))
    random_run: ${RR}
    delta_t: ${dt}
    # myshots: ${nshots}
    "
done


