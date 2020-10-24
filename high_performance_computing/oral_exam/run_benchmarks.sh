#!/bin/bash 
# First make and clean
make clean
make -B

if [[ ! -d "./bench" ]]; then
	mkdir ./bench
fi

BL=false # Run baselines - toggle this 
b_naive="bench.csv"
WARMUP=1

Sizes={1, 5, 8, 10, 15, 16, 20, 25, 32}

clear
# 100x100 
cp ./test_data/diffusion_100_100 diffusion
if [[ $BL = true ]]; then
	echo "Running baseline calculation for 100x100 case"
	hyperfine --export-csv "bench/overhead_100_100.csv" --time-unit millisecond --warmup $WARMUP "./overhead -n100000 -d0.01"
	hyperfine --export-csv "bench/naive_100_100.csv" --time-unit millisecond --warmup $WARMUP "./naive -n100000 -d0.01"
	hyperfine --export-csv "bench/handin_100_100.csv" --time-unit millisecond --warmup $WARMUP "./handin -n100000 -d0.01"
fi
echo "Running calculations for 100x100"
hyperfine --export-csv "bench/optimized_ls_100_100.csv" --time-unit millisecond --warmup $WARMUP --parameter-scan l 1 32 "./optimized_ls -n100000 -d0.01 -l{l}"
clear

# 128x128 
cp ./test_data/diffusion_128_128 diffusion
if [[ $BL = true ]]; then
	echo "Running baseline calculation for 128x128 case"
	hyperfine --export-csv "bench/overhead_128_128.csv" --time-unit millisecond --warmup $WARMUP "./overhead -n100000 -d0.01"
	hyperfine --export-csv "bench/naive_128_128.csv" --time-unit millisecond --warmup $WARMUP "./naive -n100000 -d0.01"
	hyperfine --export-csv "bench/handin_128_128.csv" --time-unit millisecond --warmup $WARMUP "./handin -n100000 -d0.01"
fi
echo "Running calculations for 128x128"
hyperfine --export-csv "bench/optimized_ls_128_128.csv" --time-unit millisecond --warmup $WARMUP --parameter-scan l 1 32 "./optimized_ls -n100000 -d0.01 -l{l}"
clear

# 100000x100
cp ./test_data/diffusion_100000_100 diffusion
if [[ $BL = true ]]; then
	echo "Running baseline calculation for 100000x100 case"
	hyperfine --export-csv "bench/overhead_100000_100.csv" --time-unit millisecond --warmup $WARMUP "./overhead -n200 -d0.6"
	hyperfine --export-csv "bench/naive_100000_100.csv" --time-unit millisecond --warmup $WARMUP "./naive -n200 -d0.6"
	hyperfine --export-csv "bench/handin_100000_100.csv" --time-unit millisecond --warmup $WARMUP "./handin -n200 -d0.6"
fi
# Local size
echo "Running calculations for 100000x100"
hyperfine --export-csv "bench/optimized_ls_100000_100.csv" --time-unit millisecond --warmup $WARMUP --parameter-scan l 1 32 "./optimized_ls -n200 -d0.6 -l{l}"
clear

# 100000x128
cp ./test_data/diffusion_100000_128 diffusion
if [[ $BL = true ]]; then
	echo "Running baseline calculation for 100000x128 case"
	hyperfine --export-csv "bench/overhead_100000_128.csv" --time-unit millisecond --warmup $WARMUP "./overhead -n200 -d0.6"
	hyperfine --export-csv "bench/naive_100000_128.csv" --time-unit millisecond --warmup $WARMUP "./naive -n200 -d0.6"
	hyperfine --export-csv "bench/handin_100000_128.csv" --time-unit millisecond --warmup $WARMUP "./handin -n200 -d0.6"
fi
# Local size
echo 
echo "Running calculations for 100000x128"
hyperfine --export-csv "bench/optimized_ls_100000_128.csv" --time-unit millisecond --warmup $WARMUP --parameter-scan l 1 32 "./optimized_ls -n200 -d0.6 -l{l}"
clear

# 10000x10000
cp ./test_data/diffusion_10000_10000 diffusion
if [[ $BL = true ]]; then
	echo "Running baseline calculation for 10000x10000 case"
	hyperfine --export-csv "bench/overhead_10000_10000.csv" --time-unit millisecond --warmup $WARMUP --max-runs 3 "./overhead -n1000 -d0.02"
	hyperfine --export-csv "bench/naive_10000_10000.csv" --time-unit millisecond --warmup $WARMUP --max-runs 3 "./naive -n1000 -d0.02"
	hyperfine --export-csv "bench/handin_10000_10000.csv" --time-unit millisecond --warmup $WARMUP --max-runs 3 "./handin -n1000 -d0.02"
fi
# Local size
echo "Running calculations for 10000x10000"
hyperfine --export-csv "bench/optimized_ls_10000_10000.csv" --time-unit millisecond --warmup $WARMUP --max-runs 3 --parameter-scan l 1 32  "./optimized_ls -n1000 -d0.02 -l{l}"
clear

cp ./test_data/diffusion_10016_10016 diffusion
if [[ $BL = true ]]; then
	echo "Running baseline calculation for 10016x10016 case"
	hyperfine --export-csv "bench/overhead_10016_10016.csv" --time-unit millisecond --warmup $WARMUP --max-runs 3 "./overhead -n1000 -d0.02"
	hyperfine --export-csv "bench/naive_10016_10016.csv" --time-unit millisecond --warmup $WARMUP --max-runs 3 "./naive -n1000 -d0.02"
	hyperfine --export-csv "bench/handin_10016_10016.csv" --time-unit millisecond --warmup $WARMUP --max-runs 3 "./handin -n1000 -d0.02"
fi
# Local size
echo "Running calculations for 10016x10016"
hyperfine --export-csv "bench/optimized_ls_10016_10016.csv" --time-unit millisecond --warmup $WARMUP --max-runs 3 --parameter-scan l 1 32  "./optimized_ls -n1000 -d0.02 -l{l}"
clear

echo "Calculations completed"
