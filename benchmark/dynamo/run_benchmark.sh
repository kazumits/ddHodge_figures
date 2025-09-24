#!/bin/bash

PYTHON="/home/kazumitsu.maehara/work/ddhodge2025/envs/dynamo/bin/python"

for ncell in 100 500 1000 5000 10000 50000 100000
    do for nbasis in 100 500 1000 1500 2000
        do for rep in 1 2 3
            do echo "Performing dynamo with M=$nbasis N=$ncell, trial $rep"
            $PYTHON ./bench_dynamo.py $nbasis $ncell | tee results/bench_dynamo_M${nbasis}_N${ncell}_R${rep}.json
        done
    done
done
