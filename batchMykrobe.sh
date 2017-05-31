#!/bin/bash

usage="usage: $0 -s (summarize) -i from_dir -o to_dir -f out_file(.tsv)"
summarize=false
from_dir=.
to_dir=./MykrobeOut
out_file=summary.tsv

export PATH=$HOME/code/Mykrobe-predictor/mccortex/bin/:$PATH

if [ $? != 0 ]
then
	   echo $usage
	   exit 2
fi

# A POSIX variable
OPTIND=1         # Reset in case getopts has been used previously in the shell.

while getopts "h?sp:i:o:f:" opt; do
    case "$opt" in
    h|\?)
        echo $usage
        exit 0
        ;;
    # v|--verbose)  verbose=1
    #     ;;
	s)	summarize=true
		;;
	p)	panel_in=$OPTARG
		;;
    i)  from_dir_in=$OPTARG
        ;;
    o)  to_dir_in=$OPTARG
        ;;
    f)  out_file_in=$OPTARG
        ;;
    esac
done

shift $((OPTIND-1))

[ "$1" = "--" ] && shift

echo "
Batch processing with Mykrobe predictor"
echo "
Input arguments"
echo "summarize=$summarize, panel_in='$panel_in', from_dir_in='$from_dir_in', 
to_dir_in='$to_dir_in', out_file_in='$out_file_in', Leftovers: $@
"
echo "Info:"

# check inputs
# check panel ('bradley-2015', 'walker-2015')
if [ ! -z $panel_in ] ; then
	if [ ! $panel_in == bradley-2015 -a ! $panel_in == walker-2015 ] ; then
		echo "Panel (-p): $panel_in invalid choice. (choose from 'bradley-2015', 'walker-2015')"
		exit 1
	else 
		panel=$panel_in
	fi
else 
	panel="bradley-2015"
	echo "Panel is $panel"
fi

# check from_dir
if [ ! -d $from_dir_in ]; then
	echo "$from_dir_in is not a valid directory"
	echo $usage
	exit 1
else
	if [ ! -z $from_dir_in ]; then
		from_dir="$( cd "${from_dir_in}" && pwd )"
		echo "Search in: $from_dir"
	fi
fi

# check to_dir
if [ -z $to_dir_in ]; then
	to_dir_in=$to_dir
	echo "Save to default: $to_dir_in"
fi
if [ -d $to_dir_in ]; then
	echo "
**WARNING**: Output directory \"$to_dir_in\" exists. 
Existing files will be overwritten and summary will include all .json files.
"
else 
	echo "Output directory \"$to_dir_in\" will be created."
	create_to_dir=true
fi

# check summary file name
if $summarize; then
	if [ ! -z $out_file_in ]; then
		if [ ! ${out_file_in: -4} == ".tsv" ]; then
			echo "$out_file_in must end with .tsv"
			echo $usage
			exit 1
		else
			out_file=$out_file_in
			echo "Summary in $to_dir_in/$out_file"
		fi
	else
		echo "No summary file given. Summary in default: $to_dir_in/$out_file"
	fi
fi

# check for redundant arguments
if [ ! -z $4 ]; then
	echo "Unused argument $4"
	echo $usage
	exit 1
fi

# Argument summary
echo "
Job summary:"
echo "Summarize: $summarize"
echo "Panel is $panel"
echo "From: $from_dir"
echo "To: $to_dir"
echo "Summary in $to_dir_in/$out_file"

read -r -p "Do you want to proceed? [y/N] " response
case $response in
	[yY][eE][sS]|[yY])
		if [ $create_to_dir ]; then
			mkdir -p $to_dir_in
		fi
		to_dir="$( cd "${to_dir_in}" && pwd )"
		echo "Output directory: $to_dir"
		;;
	*)
		exit 1
		;;
esac

## MAIN: Run Mykrobe genotype
ls $from_dir/*.fastq.gz | awk -F '_R[12]' '{print $1}' | xargs basename | uniq | parallel "mykrobe predict {} staph --panel $panel -1 $from_dir/{}*.fastq.gz > $to_dir/{.}_MpOut.json"

# Summarize all .json files into one .tsv
if $summarize; then
	python /Users/rsieber/code/my/py/json_summary.py -o $to_dir/$out_file -i $to_dir
fi
exit 0;
