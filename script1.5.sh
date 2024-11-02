#!/bin/bash
# Enable error handling and strict mode
set -eou pipefail

# Load required modules for the analysis
module load bwa       # Load BWA for sequence alignment
module load samtools  # Load Samtools for handling SAM/BAM files
module load java      # Load Java, needed for VarScan

# Prompt for the output directory name
read -p "Enter the output directory name: " OUTDIR

# Create the output directory and navigate into it
mkdir "$OUTDIR"
cd "$OUTDIR"

# Define sample name and paths for reference and annotation databases
SAMPLE=my_sample
REF_PATH=/home/scw1557/UNIX_5/BWA  # Path to the reference genome
ANNOVAR_DB=/home/scw1557/UNIX_5/annovar/humandb  # Path to ANNOVAR database

# Prompt for R1 and R2 fastq file URLs and names
read -p "Enter the URL for R1: " R1_URL
read -p "Enter the file name for R1 (e.g., R1.fastq.gz): " R1_NAME
read -p "Enter the URL for R2: " R2_URL
read -p "Enter the file name for R2 (e.g., R2.fastq.gz): " R2_NAME

# Create a directory for downloaded data
mkdir data
cd data

# Download R1 and R2 fastq files using user-provided URLs
echo
echo "Downloading R1..."
wget -O "$R1_NAME" "$R1_URL"  # Download R1 and save it as specified by user
echo
echo "Downloading R2..."
wget -O "$R2_NAME" "$R2_URL"  # Download R2 and save it as specified by user
echo

# Change permissions to read-only
chmod 444 "$R1_NAME"
chmod 444 "$R2_NAME"
cd ..

# Create a directory for mapping results
mkdir mapping
cd mapping

# Create symbolic links to the downloaded data
ln -s ../data/"$R1_NAME" R1.fastq.gz  # Link to R1 file
ln -s ../data/"$R2_NAME" R2.fastq.gz  # Link to R2 file

# Count the number of reads in R1 and R2
echo
gunzip -c R1.fastq.gz | grep -cP "^@\w+\.\w+" | xargs echo "Number of reads in R1:"
gunzip -c R2.fastq.gz | grep -cP "^@\w+\.\w+" | xargs echo "Number of reads in R2:"
echo

# Map the reads using BWA
echo
echo "Mapping data..."
bwa mem -R "@RG\tID:$SAMPLE\tPL:illumina\tSM:$SAMPLE" \
	$REF_PATH/hg38.fasta \
	R1.fastq.gz \
	R2.fastq.gz \
	> $SAMPLE\.sam \
	2>> ../log
echo "Done!"

# Convert SAM to BAM format
echo 
echo "Converting SAM to BAM..."
samtools view -bt "$REF_PATH/hg38.fasta" "$SAMPLE.sam" | samtools sort -o "$SAMPLE.bam" -  # Sort and convert to BAM format
samtools index "$SAMPLE.bam"  # Index the BAM file

# Remove the SAM file since it's no longer needed
echo
echo "Removing SAM file..."
rm "$SAMPLE.sam"  # Delete the SAM file to save space

# Summarize mapping statistics
echo
echo "Summarise mapping..."
samtools view "$SAMPLE.bam" | cut -f3 | uniq -c | sort -n  # Count unique mappings by chromosome
cd ..

# Create a directory for variants
mkdir variants
cd variants

# Generate a pileup file for variant calling
echo
echo "Run samtools mpileup..."
samtools mpileup -B -f "$REF_PATH/hg38.fasta" ../mapping/"$SAMPLE.bam" > "$SAMPLE.pileup" 2>> ../log  # Create pileup file
echo "Done!"

# Summarize the number of base positions covered in each chromosome
echo
echo "Summarise number of base positions covered in each chromosome..."
cut -f1 "$SAMPLE.pileup" | uniq -c  # Count positions per chromosome

# Representation of base coverage
echo
echo "Representation according to base..."
cut -f3 "$SAMPLE.pileup" | tr 'acgt' 'ACGT' | sort | uniq -c | sort -nr  # Count occurrences of each base
echo

# Create a coverage depth file
echo
echo -e "depth\tlocations" > coverage_depth.txt 
cut -f4 "$SAMPLE.pileup" | sort | uniq -c | sort -n -k2 | awk ' { t = $1; $1 = $2; $2 = t; print $1 "\t" $2; } ' >> coverage_depth.txt

# Call variants using VarScan
echo
echo "Call variants with varscan..."
java -jar /software/genomics/varscan/2.4.1/VarScan.jar mpileup2snp "$SAMPLE.pileup" --output-vcf 1 > "$SAMPLE.vcf" 2>> ../log  # Call SNPs
echo "Done!"

# DELETE the pileup file
echo
echo "Removing pileup file..."
rm "$SAMPLE.pileup"


# Report total number of SNPs identified
echo 
echo "Total number of SNPs identified:"
grep -vc "^#" "$SAMPLE.vcf"  # Count the number of non-header lines in the VCF file

# Create a simple spreadsheet from the VCF file
echo
echo "Creating a simple spreadsheet from the VCF file..."
grep -v "^##" "$SAMPLE.vcf" | cut -f1,2,4,5 | sed 's/\t/,/g' > snps.csv  # Extract relevant columns and save as CSV
echo "Done!"

# Create a spreadsheet with genotype information
echo "Creating a spreadsheet with genotype information..."
grep -v "^##" "$SAMPLE.vcf" | awk -F'\t' 'BEGIN {OFS=","} !/^#/ {split($10, genotype, ":"); print $1, $2, $4, $5, genotype[1]}' > snps_with_genotype.csv
echo "Done!"

# Report SNP locations
echo
echo "SNP locations:"
grep -v "^#" "$SAMPLE.vcf" | cut -f1 | uniq -c | perl -ne '/\s+(\d+)\s+(\d+)/; print "$1 snps in chromosome $2\n"'  # Count SNPs by chromosome

# Annotate the variants using ANNOVAR
echo
echo "Annotate..."
/home/scw1557/UNIX_5/annovar/convert2annovar.pl "$SAMPLE.vcf" -format vcf4old --includeinfo > "$SAMPLE.av" 2>> ../log  # Convert VCF to ANNOVAR format
/home/scw1557/UNIX_5/annovar/annotate_variation.pl -buildver hg38 -geneanno -dbtype knownGene "$SAMPLE.av" "$ANNOVAR_DB" --outfile annovar_out 2>> ../log  # Annotate variants
echo -e "chr\tposition\tref\talt\tgene\texon\tcDNA pos\taa pos\tmutation type" > final.txt

# DELETING the .av file
echo
echo "Removing the .av file..."
rm "$SAMPLE.av"

# Extract and format annotation results into a final output file
cut -f2 annovar_out.exonic_variant_function | cut -f1 -d' ' > cols2.txt  # Gene names
cut -f3 annovar_out.exonic_variant_function | cut -f1,3,4,5 -d':' --output-delimiter=" " | tr ' ' '\t' > cols3.txt  # cDNA and AA positions
cut -f4,5,7,8 annovar_out.exonic_variant_function > cols4578.txt  # Other details
paste cols4578.txt cols3.txt cols2.txt >> final.txt  # Combine all columns into final output


# Deleting unnecessary cols.txt files
rm cols4578.txt cols3.txt cols2.txt

# Set permissions for the final output file
chmod 444 final.txt

# Display the final output
cat final.txt

# Return to the previous directory
cd ..

# Return to the working directory
cd ..

# Create a tarball of the output directory
echo "Create tarball from working directory..."
tar czf "$OUTDIR.tar.gz" "$OUTDIR"  # Compress the output directory into a tar.gz file

# If you see this, the code has probably worked! Pat yourself on the back.
echo
echo "Complete!"


