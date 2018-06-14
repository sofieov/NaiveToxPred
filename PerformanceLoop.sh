#!/bin/bash                                                                                                                                                            
module load perl/5.20.2 ncbi-blast/2.6.0+ signalp/4.1c anaconda3/4.4.0

if [[ "$@" == *"--BenchmarkFile"* ]]; then
    echo > PerformanceUnited.txt #A new file is created for the performance data for all runs
    for i in 1e-10 1e-11 1e-12 1e-13 1e-14 1e-15 1e-16 1e-17 1e-20 1e-30 1e-45 1e-60 
    do
	Evalue="--Evalue $i"
	./NaiveToxPred.v1.0.py $@ $Evalue | awk '/##/{sub(/.*##/, ""); print}'	
    done >> PerformanceUnited.txt
else
    ./NaiveToxPred.v1.0.py "$@";
fi
