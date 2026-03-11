#!/bin/bash

# ================================================
# RNA-seq Quantification Pipeline
# Dataset: Airway (dexamethasone treatment)
# Tool: Salmon
# ================================================

INDEX=reference/salmon_index
DATA=data
RESULTS=results

mkdir -p $RESULTS

# ── Sample metadata ───────────────────────────
declare -A SAMPLES
SAMPLES["SRR1039508"]="untreated"
SAMPLES["SRR1039509"]="treated"
SAMPLES["SRR1039512"]="untreated"
SAMPLES["SRR1039513"]="treated"

# ── Quantification loop ───────────────────────
for SAMPLE in "${!SAMPLES[@]}"; do

    CONDITION=${SAMPLES[$SAMPLE]}

    echo "========================================="
    echo "Quantifying: $SAMPLE ($CONDITION)"
    echo "========================================="

    salmon quant \
        -i $INDEX \
        -l A \
        -1 $DATA/${SAMPLE}_1.fastq.gz \
        -2 $DATA/${SAMPLE}_2.fastq.gz \
        -p 4 \
        --validateMappings \
        --gcBias \
        -o $RESULTS/$SAMPLE

    echo "Done: $SAMPLE"
    echo ""

done

echo "========================================="
echo "All samples quantified!"
echo "Results in: $RESULTS/"
echo "========================================="
