#!/bin/bash


mydt=0.3
N=6

for i in {2..11}
do
	
    echo "$i"
    cp N${N}_configs/QC_N${N}_${i}.qasm N${N}_configs/QC_N${N}_${i}_dt${mydt}.qasm
    sed -i "s/dt/${mydt}/g" N${N}_configs/QC_N${N}_${i}_dt${mydt}.qasm
    sed -i '/barrier/d' N${N}_configs/QC_N${N}_${i}_dt${mydt}.qasm   
done

