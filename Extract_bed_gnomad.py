import pysam
import csv

"""
Script goes through the Annovar file of gnomAD variants 
and selects for variants that fall within the regions 
of the FMI bed file - targeted sequencing panel of 
~300 genes and several thousand SNP sites. 
"""

#hg19 data from gnomAD 
#(bgzip compressed and indexed)
exome_file = "hg19_gnomad211_genome_sorted.txt.gz"
index_file = "hg19_gnomad211_genome_sorted.txt.gz.tbi"
bed_file = "FMI_CDx.bed"
output_file = "filtered_genome_data.txt"

#read the intervals
bed_intervals = []
with open(bed_file, 'r') as bed:
    for line in bed:
        parts = line.strip().split()
        chrom = parts[0].replace('chr', '')  # Remove 'chr' prefix for compatibility
        start = int(parts[1])
        end = int(parts[2])
        bed_intervals.append((chrom, start, end))

# Open the exome data file using pysam 
exome_tabix = pysam.TabixFile(exome_file, index=index_file)

# Open output file for writing filtered data
with open(output_file, 'w') as output:
    writer = csv.writer(output, delimiter='\t')

    # Go through each interval in the BED file and fetch matching records from exome data
    for chrom, start, end in bed_intervals:
        try:
            # Fetch records within the specified region
            for record in exome_tabix.fetch(chrom, start, end):
                # Write to output file
                writer.writerow(record.split('\t'))
        except ValueError:
            # Handle cases where the region is not found in the exome data
            print(f"Region {chrom}:{start}-{end} not found in exome data.")

