#!/bin/sh

usage="usage: $0 -s (summarize) -p panel(.fasta) -i from_dir -o to_dir -f out_file(.tsv)"
summarize=false
recoursive=false
panel=
from_dir=.
to_dir=./MykrobeOut
out_file=summary.tsv

export PATH=$HOME/code/Mykrobe-predictor/mccortex/bin/:$PATH

if [ $? != 0 ]; then
	echo $usage
	exit 2
fi

# A POSIX variable
OPTIND=1         # Reset in case getopts has been used previously in the shell.

while getopts "h?srRp:i:o:f:" opt; do
    case "$opt" in
    h|\?)
        echo $usage
        exit 0
        ;;
    # v|--verbose)  verbose=1
    #     ;;
	s)	summarize=true
		;;
    p)  panel_in=$OPTARG
        ;;
    i)  from_dir_in=$OPTARG
        ;;
    o)  to_dir_in=$OPTARG
        ;;
    f)  out_file_in=$OPTARG
        ;;
    R|r)  recoursive=true
        ;;
    esac
done

shift $((OPTIND-1))

[ "$1" = "--" ] && shift

echo "
Batch processing with Mykrobe genotype"
echo "
Input arguments"
echo "summarize=$summarize, panel_in='$panel_in', from_dir_in='$from_dir_in', 
to_dir_in='$to_dir_in', out_file_in='$out_file_in', Leftovers: $@
"
echo "Info:"


# check inputs
# check panel (.fasta)
if [ -z $panel_in -o ! -f $panel_in -o ! ${panel_in: -6} == ".fasta" ] ; then
	echo "Panel (-p) is required and must be a .fasta file. \"$panel_in\" is not valid"
	echo $usage
	exit 1
else
	# panel=$panel_in
	panel="$( cd "$(dirname $panel_in)" && pwd )"/"$(basename $panel_in)"
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
			out_file=$(basename $out_file_in)
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
if $recoursive; then 
		find $from_dir -name '*.fastq.gz' | awk -F '_[Rr][12]' '{print $1}' | sort | uniq | parallel "mykrobe genotype {.} $panel -1 {}*.fastq.gz > $to_dir/{/.}_MpOut.json"
	else
		find $from_dir -d 1 -name '*.fastq.gz' | awk -F '_[Rr][12]' '{print $1}' | sort | uniq  | parallel "mykrobe genotype {.} $panel -1 {}*.fastq.gz > $to_dir/{/.}_MpOut.json"
fi
# ls $from_dir/*.fastq.gz | awk -F '_[Rr][12]' '{print $1}' | xargs basename | uniq  | parallel "mykrobe genotype {} $panel -1 $from_dir/{}*.fastq.gz > $to_dir/{.}_MpOut.json"
# find $from_dir/ -name *.fastq.gz | awk -F '_R[12]' '{print $1}' | xargs basename | uniq  | parallel "mykrobe genotype {} $panel -1 $from_dir/{}*.fastq.gz > $to_dir/{.}_MpOut.json"
# ls $from_dir/*.fastq.gz | awk -F '_R[12]' '{print $1}' | xargs basename | uniq  | parallel "mykrobe genotype {} $panel -1 $from_dir/{}*.fastq.gz > $to_dir/{.}_MpOut.json"

# Summarize all .json files into one .tsv
if $summarize; then
	python $HOME/code/my/py/json_summary_genotype.py -o $to_dir/$out_file -i $to_dir
	echo "Make sure to consider coverage when interpreting results."
fi
exit 0;
