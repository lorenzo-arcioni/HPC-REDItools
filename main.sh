#!/bin/bash
#Main script

#Use it to automate the running of the start.sh
usage="$(basename "$0") -i inputpath (bam/bams DIR path) [-h] [-r reference genome] [--slurm] [--htcondor] [-t threads]
#-- parallel REDItools executor"

man="This is a parallelization tool for reditools python module. Make sure you've installed the right version (3) of reditools (on the launch machine). 
Note that, for an error free execution, the command 'python -m reditools' must run without errors on the launching shell.

For more information: https://github.com/RediTools/reditools
"

wlm='none'

while getopts 'i:hr:t:-:' option;
do
    case "${option}" in
        h) echo $usage
		   echo ""
		   printf "$man"
           exit 0;;
        i) inputpath=$OPTARG;;
		t) threads=$OPTARG;;
        r) reference=$OPTARG;;
        
		-)
			case "${OPTARG}" in
				slurm)
				wlm='slurm'
				;;
				htcondor)
				wlm='htcondor'
				;;
				* )
				echo printf "illegal option: -%s\n" "$OPTARG" >&2
					echo "$usage" >&2
					exit 1
				;;
			esac
			;;

        :) printf "missing argument for -%s\n" "$OPTARG" >&2
			echo "$usage" >&2
	   		exit 1;;
		?) printf "illegal option: -%s\n" "$OPTARG" >&2
            echo "$usage" >&2
            exit 1;;
		
    esac
done

if [ -z "$inputpath" ]
then
	echo "Insert an input file directory path (BAM/BAMS)" >&2
	echo $usage
	exit 2
fi

if [ -z "$threads" ]
then
	echo "Insert the number of threads" >&2
	echo $usage
	exit 2
fi

if [ -z "$reference" ]
then
	echo "Insert the reference genome" >&2
	echo $usage
	exit 2
fi

echo "Current settings"
echo
echo Input path: $inputpath
echo Threads: $threads
echo Reference: $reference
echo Workload manager: $wlm

echo
echo "Continue? (y/n)"
read choice

if [ "$choice" = "y" ]
then
	python3  creator.py -i "$inputpath" -t $threads -r "$reference" -w "$wlm"

	echo "All codes are correctly generated!"
	echo "Do you want to: (enter the number choice and press ENTER)"
	echo "1) Exec your code rigth now (on this machine)"
	echo "2) Generate tar file to upload on a remote machine"
	echo "3) Exit"
	read choice

	if [ "$choice" = "1" ]; then
			case "$wlm" in
				slurm)
				sbatch --export=ALL,threads=$threads,inputpath=$inputpath,reference=$reference start.sh
				;;
				htcondor)
				echo "Work in progress!! Please select an other workload manager."
				;;
				none)
				export threads=$threads inputpath=$inputpath reference=$reference && ./start.sh
				;;
			esac
	else if [ "$choice" = "2" ]; then
		case "$wlm" in
			slurm)
			tar -cf hpc_reditools.tar read.py start.sh control_script.sh
			;;
			htcondor)
			echo "Work in progress!! Please select an other workload manager."
			;;
			none)
			tar -cf hpc_reditools.tar read.py start.sh control_script.sh
			;;
		esac 
		rm  read.py start.sh control_script.sh

		fi
		
	fi

fi
