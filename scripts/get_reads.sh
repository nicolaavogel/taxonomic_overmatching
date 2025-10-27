module load seqtk

# First loop: Extract read IDs
for basedir in $(cut -f1-10 -d/ library_list.txt); do
    while read id; do
        lca=$(find "$basedir/results/metadmg/lca/" -name "*.lca.gz")
        genus=$(echo "$id" | cut -f2 -d":" | tr -d '"')
        ids="../readids/$(basename "$lca" _collapsed.lca.gz).$genus.txt"

        # Debugging: Check values
        echo "Processing: basedir=$basedir, lca=$lca, genus=$genus, id=$id"
        
        # Fix: Use double quotes around variables and avoid single quotes inside grep
        zcat "$lca" | grep -F "$id" | cut -f1 > "$ids"
        zcat "$lca" | grep -F "$id" > "../readids/$(basename "$lca" _collapsed.lca.gz).$genus.lca"

        # Check if IDs were written
        if [[ ! -s "$ids" ]]; then
            echo "Warning: No matching reads for $id in $lca"
        fi
    done < taxa
done


# Second loop: Extract sequences
# for basedir in $(cut -f1-10 -d/ library_list.txt); do
#     fastq=$(find "$basedir/results/reads/low_complexity/" -name "*_collapsed.fastq.gz")
#
#     while read id; do
#         lca=$(find "$basedir/results/metadmg/lca/" -name "*.lca.gz")
#         genus=$(echo "$id" | cut -f2 -d":" | tr -d '"')
#         ids="../readids/$(basename "$lca" _collapsed.lca.gz).$genus.txt"
#
#         # Debugging: Check values
#         echo "Processing FASTQ extraction: fastq=$fastq, ids=$ids"
#
#         # Fix: Ensure ID file is not empty before running seqtk
#         if [[ -s "$ids" ]]; then
#             seqtk subseq "$fastq" "$ids" > "../fastqs/$(basename "$lca" _collapsed.lca.gz).$genus.fastq"
#         else
#             echo "Warning: Skipping seqtk for $genus - ID file is empty."
#         fi
#     done < taxa
# done
