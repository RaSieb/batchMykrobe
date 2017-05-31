# batchMykrobe
These scripts allow to run the commandline version of [Mykrobe Predictor](https://github.com/iqbal-lab/Mykrobe-predictor) in parallel for batches of sequence reads on macOSX. 

## batchMykrobe.sh
batchMykrobe.sh will allow to run the command `mykrobe predict staph` on a set of samples with paired read files (fastq.gz). Using the -s flag, the results will be summarized in a table (requires adjusting the path to json_summary.py). 

## json_summary.py
Summarizes the outputs from batchMykrobe.sh.

## batchMykrobe_genotype.sh
batchMykrobe.sh will allow to run the command `mykrobe genotype` on a set of samples with paired raw read files (fastq.gz). Using the -s flag, the results will be summarized in a table (requires adjusting the path to json_summary_genotype.py). 

## json_summary_genotype.py
Summarizes the outputs from batchMykrobe_genotype.sh.

## Dependencies
All scripts are tested on macOSX. Adjustments might be necessary for adapting them to linux. They depend on a working installation of the commandline version of [Mykrobe Predictor](https://github.com/iqbal-lab/Mykrobe-predictor) with its dependencies. Bash scripts also depend on a working installation of [GNU parallel](https://www.gnu.org/software/parallel/)
