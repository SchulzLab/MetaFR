#! /bin/bash

#### NF variables ####
path="$1" #input directory with binned activity data
model_type="$2" #for now only "RF" is supported
seed_value="$3"  #random seed for reproducibility
expr_threshold="$4" #threshold for variance of log2 transformed gene expression
min_nonzero_expr="$5" #minimum percentage of non-zero expression counts
min_mc_total="$6" #minimum number metacells required
min_nonzero_feature_mat="$7" #minimum percentage non-zero feature values
test_cells_path="$8" #path to file with test metacell - nf pipeline train-test split in activity binning step
plot_dir="$9" #path to output directory for plots
log_file_path="${10}" #path to log file
models_out_dir="${11}" #path to output directory for model files
performance_out_dir="${12}" #path to output directory for model performance files
input_extension="${13}" #input file extension, e.g. "_lsi_activity_top_100.txt.gz
n_jobs="${14}" #maximum number of parallel jobs to run
timeout="${15}" #timeout duration for each job, e.g. "10m"

########################
#### Set environment variables to limit the number of threads used by libraries ####
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

########################
#### Set fixed paths ####
python_script=/projects/single_cell_stitchit/work/binned_RF/nf_pipeline/bin/setup_RF_individualGene_lsi_pbmc_meta_cell_all_features.py
########################
#path="/projects/single_cell_stitchit/work/binned_RF/input/paired_PBMC10k/mc100_1MB500bp_binned_lsi_activity_q0"
START_TIME=$(date +%s)

#TODO: optimize processing of input files
echo "Processing files from: $path"
declare -a x=()
for file in $path"/"*"${input_extension}"; do
    fn=$(basename "$file")         
    fn_no_ext="${fn%.*}"     
    gene="${fn_no_ext%%_*}"   
    echo $gene    

    if [ -n "$gene" ]; then
        x+=("$gene")
    fi
done

run_job() {
	local filename="$1" #file name = gene name
	local input_file="$path/${filename}${input_extension}"
	local logfile="logs/${filename}.log"
	mkdir -p "$(dirname "$logfile")"
	timeout "$timeout" python3 "$python_script" \
        "$input_file" \
        "$model_type" \
        "$seed_value" \
        "$expr_threshold" \
        "$min_nonzero_expr" \
        "$min_mc_total" \
        "$min_nonzero_feature_mat" \
        "$test_cells_path" \
        "$plot_dir" \
        "$log_file_path" \
        "$models_out_dir" \
        "$performance_out_dir" \
        2>&1 | tee -a "$logfile"
}

mkdir -p logs
for filename in "${x[@]}"; do
	while (( $(jobs -rp | wc -l) >= n_jobs )); do
		sleep 1  
	done
	echo "$filename"
	run_job "$filename"  &
done
wait

### Calculate total runtime ###
END_TIME=$(date +%s)
TOTAL_TIME=$((END_TIME - START_TIME))

HOURS=$((TOTAL_TIME / 3600))
MINUTES=$(( (TOTAL_TIME % 3600) / 60 ))
SECONDS=$((TOTAL_TIME % 60))

echo "----------------------------------------"
echo "Total runtime: ${HOURS}h ${MINUTES}m ${SECONDS}s"
echo "----------------------------------------"
