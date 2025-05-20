#!/bin/sh
echo "Getting FastQC report for $1"
fastqc $1

echo "Indexing ref $2..."
minimap2 -d ref.mmi $2

echo "Aligning $1..."
minimap2 -a ref.mmi $1 > aligned.sam

echo "Converting to .bam..."
samtools view -b aligned.sam > aligned.bam

echo "Getting stats..."
samtools flagstat aligned.bam > alignment_report.txt

mapped_percent=$(grep -oP '\d+\.\d+(?=% mapped)' alignment_report.txt)
if (( $(echo "$mapped_percent > 90" | bc -l) )); then
    echo "OK"
else
    echo "Not OK"
fi