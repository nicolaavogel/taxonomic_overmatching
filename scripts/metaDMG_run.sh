name=$1
THREADS=30
THREADSP=8  # Adjust parallel jobs accordingly
OUTPUT_PATH="/projects/caeg/data/pp_analysis/Iceland/species_id/${name}"  

SAMPLE_LIST=$2

log_step() {
    echo "$(date) - $1"
}



log_step "Sorting merged BAM file..."
cat "$SAMPLE_LIST" | parallel -j "$THREADSP" "samtools sort -n -@ $THREADS -m 10G -o $OUTPUT_PATH/{}.sort.merged.bam" "$OUTPUT_PATH/{}.merged.bam"


log_step "Running taxonomic classification with metaDMG..."
cat "$SAMPLE_LIST" | parallel -j "$THREADSP" "/projects/wintherpedersen/apps/metaDMG_28Nov24/metaDMG-cpp lca \
  --names /datasets/caeg_dataset/taxonomy/20250210/names.dmp \
  --nodes /datasets/caeg_dataset/taxonomy/20250210/nodes.dmp \
  --acc2tax /projects/caeg/people/bfj994/hashilan_marsh/newDBall.acc2taxid.gz \
  --sim_score_low 0.95 --sim_score_high 1.0 --how_many 30 --weight_type 1 \
  --fix_ncbi 0 --threads 10 --filtered_acc2tax $OUTPUT_PATH/{}.acc2tax \
  --bam $OUTPUT_PATH/{}.sort.merged.bam --out_prefix $OUTPUT_PATH/{}.sort.merged"

log_step "Running damage estimation with metaDMG..."
cat "$SAMPLE_LIST" | parallel -j "$THREADSP" "/projects/wintherpedersen/apps/metaDMG_28Nov24/metaDMG-cpp dfit \
	  $OUTPUT_PATH/{}.sort.merged.bdamage.gz --threads 6 \
  	  --names /datasets/caeg_dataset/taxonomy/20250210/names.dmp \
  	  --nodes /datasets/caeg_dataset/taxonomy/20250210/nodes.dmp \
      --showfits 2 --nopt 10 \
      --nbootstrap 20 --doboot 1 --seed 1234 --lib ds \
      --out_prefix $OUTPUT_PATH/{}.sort.merged"


log_step "Aggregating lca and dfit metaDMG..."
cat "$SAMPLE_LIST" | parallel -j "$THREADSP" "/projects/wintherpedersen/apps/metaDMG_28Nov24/metaDMG-cpp aggregate \
	  $OUTPUT_PATH/{}.sort.merged.bdamage.gz \
  	  --names /datasets/caeg_dataset/taxonomy/20250210/names.dmp \
  	  --nodes /datasets/caeg_dataset/taxonomy/20250210/nodes.dmp \
      --lcastat $OUTPUT_PATH/{}.sort.merged.stat.gz --dfit $OUTPUT_PATH/{}.sort.merged.dfit.gz --out_prefix $OUTPUT_PATH/{}.sort.merged.agg"
