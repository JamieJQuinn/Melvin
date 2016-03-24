#/bin/bash
nZ=101
nN=51
dt=5e-6
Pr=0.5
a=3
T=0.5
S=0.01
Ra=1e5
ICs=Stable_Ra1e5_ICn1_nZ101_nN51_PERTURBED
for amp in 0.025 0.05 0.075 0.1; do
	python addPerturbation.py 101 51 2 $amp Stable_Ra1e5_ICn1_nZ101_nN51 Stable_Ra1e5_ICn1_nZ101_nN51_PERTURBED
	folder=StabilityAnalysis/n1$amp
	mkdir -p $folder
	rm -f $folder/*
	./box "-nZ" $nZ "-nN" $nN "-dt" $dt "-Ra" $Ra "-Pr" $Pr "-a" $a "-T" $T "-S" $S "-i" $ICs "-o" $folder/ | tee $folder/log
done
