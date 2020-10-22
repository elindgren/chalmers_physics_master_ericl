#!/bin/bash 
# First make and clean
make clean
make -B

if [[ ! -d "./bench" ]]; then
	mkdir ./bench
fi

BL=true # Run baselines - toggle this 
b_naive="bench.csv"
WARMUP=1

clear
# 100x100
cp ./test_data/diffusion_100_100 diffusion
if [[ $BL ]]; then
	echo "Running baseline calculation for 100x100 case"
	hyperfine --export-csv "bench/naive_100_100.csv" --time-unit millisecond --warmup $WARMUP "./naive -n100000 -d0.01"
	hyperfine --export-csv "bench/handin_100_100.csv" --time-unit millisecond --warmup $WARMUP "./handin -n100000 -d0.01"
fi
hyperfine --export-csv "bench/optimized_100_100.csv" --time-unit millisecond --warmup $WARMUP "./optimized -n100000 -d0.01"

clear
# 100000x100
cp ./test_data/diffusion_100000_100 diffusion
if [[ $BL ]]; then
	echo "Running baseline calculation for 100000x100 case"
	hyperfine --export-csv "bench/naive_100000_100.csv" --time-unit millisecond --warmup $WARMUP "./naive -n200 -d0.6"
	hyperfine --export-csv "bench/handin_100000_100.csv" --time-unit millisecond --warmup $WARMUP "./handin -n200 -d0.6"
fi
hyperfine --export-csv "bench/optimized_100000_100.csv" --time-unit millisecond --warmup $WARMUP "./optimized -n200 -d0.6"

clear
# 10000x10000
cp ./test_data/diffusion_10000_10000 diffusion
if [[ $BL ]]; then
	echo "Running baseline calculation for 10000x10000 case"
	hyperfine --export-csv "bench/naive_10000_10000.csv" --time-unit millisecond --warmup $WARMUP --max-runs 3 "./naive -n1000 -d0.02"
	hyperfine --export-csv "bench/handin_10000_10000.csv" --time-unit millisecond --warmup $WARMUP --max-runs 3 "./handin -n1000 -d0.02"
fi
hyperfine --export-csv "bench/optimized_100_100.csv" --time-unit millisecond --warmup $WARMUP --max-runs 3 "./optimized -n1000 -d0.02"
