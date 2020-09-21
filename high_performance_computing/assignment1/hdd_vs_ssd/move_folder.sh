setup_include() {
	if [[ ! -d "./include" ]]; then
		echo "'include' folder not in current directory. Copying..."
		cp -r /usr/include .
		echo "'include' copied to current location."
	else
		echo "'include' folder is already in current directory."
	fi
}

copy_folder(){
	start=`date +%s`
	for((i=0; i<10; i++)); do
		cp -r include include_copy
	done
	end=`date +%s`
	time=$((end-start))
	avgtime=$(bc <<< "scale=2;$time/10")
	echo "Average copy time for $1: $avgtime s" 
}

### HDD

# First copy the folder to home directory if it doesn't already exist
setup_include
# Copy include to another folder, include_copy, 10 times. Measure the time. 
copy_folder "HDD"
# Remove include folder:
rm -r ./include
rm -r ./include_copy
echo "'include' and 'include_copy' removed from current location."

### SSD

# Move to new location
cd /run/mount/scratch/hpcuser272
# First copy the folder to home directory if it doesn't already exist
setup_include
# Copy include to another folder, include_copy, 10 times. Measure the time. 
copy_folder "SSD"
# Remove include folder:
rm -r ./include
rm -r ./include_copy
echo "'include' and 'include_copy' removed from current location."

