# ctDNA Sequencing in HCC-TACE Response Analysis

[![BioProject](https://img.shields.io/badge/BioProject-PRJNA1199049-brightgreen.svg)](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1199049)
[![License](https://img.shields.io/badge/License-Academic%20Use-yellow.svg)](LICENSE)

## Overview

This repository contains the computational analysis pipeline used in our study of circulating tumor DNA (ctDNA) as a biomarker for clinical outcomes following transarterial chemoembolization (TACE) in patients with hepatocellular carcinoma (HCC).


The primary purpose of this repository is to provide the code used for the analyses presented in the associated publication, facilitating the reproducibility of the study's findings.

## Related Publication

**Title:** Deep sequencing of circulating tumour DNA as a biomarker of clinical outcome to transarterial chemoembolisation in hepatocellular carcinoma

**Authors:** Rohini Sharma, Sultan N. Alharbi, Ksenia Ellum, Leila Motedayen-Aval, Andrea Casadei-Gardini, David J. Pinato, Dominik Bettinger, Bertram Bengsch, Rishi Patel, Joanne Evans

## Data Availability

Raw sequencing data generated in this study are available from the NCBI Sequence Read Archive (SRA):

- **BioProject:** [PRJNA1199049](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1199049)

## Pipeline Overview

Our analysis pipeline consists of the following key steps:

1. **UMI Processing & QC**: Extract UMIs from trimmed PE reads and perform quality control
2. **Alignment & Sorting**: Map reads to GRCh38 reference genome and sort alignments
3. **Duplicate Marking**: UMI-aware deduplication
4. **Base Quality Recalibration**: Improve base quality scores
5. **Variant Calling & Filtering**: Identify and filter somatic variants
6. **Variant Annotation**: Functionally annotate variants
7. **Downstream Analysis**: Statistical analysis and visualization

## Prerequisites

### Software Requirements

| Software | Version | Purpose |
|----------|---------|---------|
| [UMI-tools](https://github.com/CGATOxford/UMI-tools) | v1.0.0 | UMI extraction and processing |
| [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) | v0.11.9 | Quality control |
| [BWA](https://github.com/lh3/bwa) | v0.7.17 | Read alignment |
| [Samtools](http://www.htslib.org/) | v1.14+ | SAM/BAM manipulation |
| [GATK](https://github.com/broadinstitute/gatk) | v4.2.6.1 | Variant calling & filtering |
| [R](https://www.r-project.org/) | v4.2+ | Statistical analysis |
| [Python](https://www.python.org/) | v3.9+ | Data processing |

### Data Requirements

- **Reference Genome**: GRCh38 FASTA and BWA index files
- **Known Sites Files**: dbSNP VCF, gnomAD VCF, Panel of Normals (PON) VCF
- **Target Interval BED file**: HCC panel target regions
- **Funcotator Data Sources**: For variant annotation
* **Reference Data:** (Links point to Google Cloud Storage buckets commonly used for GATK resources. Download may require `gsutil` or browser access. Ensure corresponding index files are also downloaded/generated.)
    * Human Reference Genome GRCh38 FASTA file.
    * BWA Index files for GRCh38 - (Must be generated from the FASTA file using `bwa index`)
    * dbSNP VCF file: `Homo_sapiens_assembly38.dbsnp138.vcf` ([Link - GATK Bundle](https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf), requires `.vcf.idx` index file too) for GATK BQSR
    * Panel of Normals (PON) VCF file: `1000g_pon.hg38.vcf.gz` ([Link](https://storage.googleapis.com/gatk-best-practices/somatic-hg38/1000g_pon.hg38.vcf.gz), requires `.tbi` index file too) for GATK Mutect2
    * Germline Resource VCF file: `af-only-gnomad.hg38.vcf.gz` ([Link](https://storage.googleapis.com/gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz), requires `.tbi` index file too) for GATK Mutect2
    * Target Regions BED file (`.bed`) - Included in `HCC_panel_target_regions.bed`
    * Funcotator Data Sources (Version `v1.7.20200521s` used in example GATK command): `funcotator_dataSources.v1.7.20200521s.tar.gz` ([Link](https://console.cloud.google.com/storage/browser/broad-public-datasets/funcotator)) - Required for GATK Funcotator annotation.


### Installation / Setup

### Step 1: Clone the repository

```bash
git clone https://github.com/username/ctDNA-HCC-TACE-Analysis.git
cd ctDNA-HCC-TACE-Analysis
```
### Step 2: Create the environment

```bash
# Using Conda (recommended)
conda env create -f environment.yml

# Activate the environment
conda activate ctdna-hcc
```

### Step 3: Obtain raw data

Download the raw sequencing data from NCBI SRA:

```bash
# Install SRA Toolkit if needed
conda install -c bioconda sra-tools

# Download data
prefetch PRJNA1199049
fasterq-dump --split-files SRR*
```

### Step 4: Configure pipeline

Edit the `config/pipeline_config.yaml` file to specify paths to reference files and adjust parameters.



## Usage

### Workflow Steps:


1.  **UMI Extraction & QC:** `bash scripts/preprocess_reads.sh ...` [TODO: Add specific example]

   
2.  **Alignment & Sorting (BWA + Samtools):**
*This script processes a list of sample prefixes (one per line, provided as input `$1`), aligning corresponding paired-end FASTQ files (`${prefix}_R1.fastq.gz`, `${prefix}_R2.fastq.gz`) and creating sorted BAM files.*
    ```bash
    #!/bin/bash
    
    BWA_INDEX_BASE="/path/to/reference/GRCh38.fa" # BWA index prefix
    OUTPUT_BAM_DIR="/path/to/output/BAM"
    TEMP_DIR="/path/to/temp/sorting_temp" # Temporary directory for sorting
    THREADS=8 # Number of threads for BWA and Samtools
    MEM_PER_THREAD="4G" # Memory per thread for Samtools sort

    mkdir -p $OUTPUT_BAM_DIR $TEMP_DIR

    # --- Input Validation ---
    if [ -z "$1" ]; then
      echo "Usage: $0 <file_with_sample_prefixes>"
      exit 1
    fi
    INPUT_LIST_FILE="$1"

    # --- Processing Loop ---
    while IFS= read -r line || [[ -n "$line" ]]; do
        echo "Processing Sample Prefix: $line"
        
        # Construct FASTQ paths (Assumes files are in the current directory or specified path)
        FASTQ1="${line}_R1.fastq.gz" # [Adjust path/naming if needed]
        FASTQ2="${line}_R2.fastq.gz" # [Adjust path/naming if needed]
        
        # Check if FASTQ files exist
        if [[ ! -f "$FASTQ1" || ! -f "$FASTQ2" ]]; then
            echo "ERROR: FASTQ files not found for prefix $line ($FASTQ1, $FASTQ2)"
            continue # Skip to next sample
        fi

        # Define Read Group string (Essential for downstream tools like GATK)
        # ID: Unique identifier (e.g., flowcell.lane)
        # SM: Sample Name
        # LB: Library ID
        # PL: Platform (e.g., ILLUMINA)
        RG_STRING=$(echo -e "@RG\\tID:FLOWCELL1.LANE1\\tSM:$line\\tLB:TARGTED-SEQ\\tPL:ILLUMINA")
        
        # Define Output BAM path
        OUTPUT_BAM="$OUTPUT_BAM_DIR/${line}.sorted.bam"

        echo "Aligning $FASTQ1 & $FASTQ2 with BWA..."
        # bwa mem: Aligns reads
        # -R: Adds Read Group information
        # -t: Number of threads
        # samtools sort: Sorts alignments by coordinate
        # -@: Number of sorting threads
        # -m: Max memory per thread
        # -T: Temporary file prefix
        # -o: Output BAM file
        # - : Input is read from standard input (pipe from bwa mem)
        bwa mem -R "$RG_STRING" -t $THREADS "$BWA_INDEX_BASE" "$FASTQ1" "$FASTQ2" | \
        samtools sort -@ $THREADS -m $MEM_PER_THREAD -T "$TEMP_DIR/${line}_sort_temp" -o "$OUTPUT_BAM" -
        
        # Check samtools exit status
        if [ $? -eq 0 ]; then
            echo "Successfully created sorted BAM: $OUTPUT_BAM"
            # Optional: Index the BAM file immediately
            # samtools index -@ $THREADS $OUTPUT_BAM
        else
            echo "ERROR during alignment/sorting for $line. Check logs."
        fi
        echo "---"

    done < "$INPUT_LIST_FILE"

    echo "Alignment and sorting complete."
    ```
    *[Note: Replace placeholder paths. Adjust threads (`-t`, `-@`) and memory (`-m`) based on your system. Ensure BWA index exists.]*


3.  **Mark Duplicates (Je Example):**
    *This command processes a single sorted BAM file (path assumed to be in `$line`) to mark duplicates using UMIs.*
    ```bash
    # --- Configuration ---
    JE_JAR_PATH="[path/to/je.jar]" # Path to the Je JAR file
    JAVA_MEM="[e.g., 8g]" # Java heap memory allocation
    INPUT_BAM="$line" # Assumes $line holds the path to the sorted BAM
    OUTPUT_PREFIX="/path/to/output/MarkDuplicated_Je/$(basename ${INPUT_BAM%.sorted.bam})" # Base for output files

    # --- Command ---
    echo "Marking duplicates for $INPUT_BAM..."
    java -Xmx$JAVA_MEM -jar $JE_JAR_PATH markdupes \
      INPUT=$INPUT_BAM \
      OUTPUT="${OUTPUT_PREFIX}MarkDuplicated_Je.bam" \
      METRICS_FILE="${OUTPUT_PREFIX}MarkDuplicated_Je.txt" \
      MM=1 `# Max mismatches allowed in UMI` \
      SLOTS=-1 `# Parameter specific to Je` \
      ASSUME_SORTED=true `# Input BAM is coordinate sorted` \
      REMOVE_DUPLICATES=true `# Remove duplicate reads instead of just marking them`
      
    
    if [ $? -eq 0 ]; then
        echo "Successfully marked/removed duplicates for $INPUT_BAM"
    else
        echo "ERROR during duplicate marking for $INPUT_BAM"
    fi
    ```
    
   *[Note: Replace `[path/to/je.jar]`, `[e.g., 8g]`, and output paths. Ensure the input `$line` variable correctly points to the BAM file from the previous step.]*


4.  **Merge BAMs (if needed):**
    ```bash
    # Example using samtools merge if multiple BAMs per sample exist
    samtools merge -@ $THREADS /path/to/output/merged/${sample_id}.merged.bam /path/to/input/${sample_id}.run1.bam /path/to/input/${sample_id}.run2.bam
    ```

5.  **Downstream Analysis (GATK Example Script):**
    *The following script demonstrates BQSR, Mutect2 variant calling, filtering, and annotation steps using GATK 4.2.6.1. It expects a file containing a list of input BAM file paths as its first argument (`$1`). Adjust paths to tools, reference 
    files, and output directories.*
```bash
#!/bin/bash
set -e  # Exit on error

# --- Configuration ---
GATK_PATH="/path/to/gatk-4.2.6.1/gatk"
GENOME="/path/to/reference/GRCh38.fa"
DBSNP="/path/to/reference/Homo_sapiens_assembly38.dbsnp138.vcf"
BED_FILE="/path/to/HCC_panel_target_regions.bed"
PON="/path/to/reference/1000g_pon.hg38.vcf.gz"
GERMLINE_RESOURCE="/path/to/reference/af-only-gnomad.hg38.vcf.gz"
FUNCOTATOR_DATASOURCES="/path/to/reference/funcotator_dataSources.v1.7.20200521s"
JAVA_OPTS="-Xmx8g" # Example Java options for GATK

# --- Output Directories ---
RECAL_TABLE_DIR="/path/to/output/Recalibrated_dataTable"
RECAL_BAM_DIR="/path/to/output/Recalibrated_BAM"
MUTECT_BAM_DIR="/path/to/output/Final_VCF_Mutect2/OutputBAM"
MUTECT_VCF_DIR="/path/to/output/Final_VCF_Mutect2/RawData_VCF"
F1R2_DIR="/path/to/output/Final_VCF_Mutect2/F1R2"
ORIENT_MODEL_DIR="/path/to/output/Final_VCF_Mutect2/F1R2"
FILTERED_VCF_DIR="/path/to/output/Final_VCF_Mutect2/OptimizinngmutectFiltering"
ANNOTATED_VCF_DIR="/path/to/output/Final_VCF_Mutect2/AnnotatedVCF"
MAF_DIR="/path/to/output/Final_VCF_Mutect2/MAF_output"

# --- Error handling function ---
error_exit() {
    echo "Error: ${1:-"Unknown Error"}" 1>&2
    exit 1
}

# --- Logging function ---
log_message() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1"
}

# --- Initial Validation ---
# Check if GATK exists
command -v $GATK_PATH >/dev/null 2>&1 || error_exit "GATK not found at $GATK_PATH"

# Check required input files
for file in "$GENOME" "$DBSNP" "$BED_FILE" "$PON" "$GERMLINE_RESOURCE" "$FUNCOTATOR_DATASOURCES"; do
    [ -f "$file" ] || error_exit "Required file not found: $file"
done

# Create output directories
for dir in "$RECAL_TABLE_DIR" "$RECAL_BAM_DIR" "$MUTECT_BAM_DIR" "$MUTECT_VCF_DIR" \
           "$F1R2_DIR" "$ORIENT_MODEL_DIR" "$FILTERED_VCF_DIR" "$ANNOTATED_VCF_DIR" "$MAF_DIR"; do
    mkdir -p "$dir" || error_exit "Cannot create directory: $dir"
done

# --- Input Validation ---
if [ -z "$1" ]; then 
    echo "Usage: $0 <file_with_bam_paths>"
    exit 1
fi

INPUT_LIST_FILE="$1"
[ -f "$INPUT_LIST_FILE" ] || error_exit "Input list file not found: $INPUT_LIST_FILE"

# --- Processing Loop ---
while IFS= read -r line || [[ -n "$line" ]]; do
    if [ -z "$line" ]; then continue; fi  # Skip empty lines
    
    log_message "Processing GATK steps for BAM: $line"
    
    # Validate input BAM exists
    [ -f "$line" ] || error_exit "Input BAM file not found: $line"
    
    # Extract base filename for output naming
    base_name=$(basename "$line" .bam)
    base_name=${base_name%_MarkDuplicated_Je}

    # --- GATK Steps ---
    # Step 5.1: Generate recalibration table
    log_message "Step 5.1: BaseRecalibrator - Generate recalibration table"
    $GATK_PATH --java-options "$JAVA_OPTS" BaseRecalibrator \
        -I "$line" \
        -R "$GENOME" \
        --known-sites "$DBSNP" \
        -L "$BED_FILE" \
        -O "$RECAL_TABLE_DIR/${base_name}_recal_data.table" || error_exit "BaseRecalibrator failed"

    # Step 5.2: Apply recalibration to reads
    log_message "Step 5.2: ApplyBQSR - Apply recalibration to reads"
    $GATK_PATH --java-options "$JAVA_OPTS" ApplyBQSR \
        -I "$line" \
        -R "$GENOME" \
        --bqsr-recal-file "$RECAL_TABLE_DIR/${base_name}_recal_data.table" \
        -L "$BED_FILE" \
        -O "$RECAL_BAM_DIR/${base_name}_recalibrated.bam" || error_exit "ApplyBQSR failed"
    
    recal_bam="$RECAL_BAM_DIR/${base_name}_recalibrated.bam"
    
    # Step 5.3: Call somatic variants
    log_message "Step 5.3: Mutect2 - Call somatic variants"
    $GATK_PATH --java-options "$JAVA_OPTS" Mutect2 \
        -R "$GENOME" \
        -I "$recal_bam" \
        -L "$BED_FILE" \
        --panel-of-normals "$PON" \
        --germline-resource "$GERMLINE_RESOURCE" \
        --f1r2-tar-gz "$F1R2_DIR/${base_name}_f1r2.tar.gz" \
        --bam-output "$MUTECT_BAM_DIR/${base_name}_mutect2.bam" \
        -O "$MUTECT_VCF_DIR/${base_name}_unfiltered.vcf" || error_exit "Mutect2 failed"

    # Step 5.4: Model orientation bias
    log_message "Step 5.4: LearnReadOrientationModel - Model orientation bias"
    $GATK_PATH LearnReadOrientationModel \
        -I "$F1R2_DIR/${base_name}_f1r2.tar.gz" \
        -O "$ORIENT_MODEL_DIR/${base_name}_read-orientation-model.tar.gz" || error_exit "LearnReadOrientationModel failed"

    # Step 5.5: Filter variant calls
    log_message "Step 5.5: FilterMutectCalls - Filter variant calls"
    $GATK_PATH FilterMutectCalls \
        -R "$GENOME" \
        -V "$MUTECT_VCF_DIR/${base_name}_unfiltered.vcf" \
        -L "$BED_FILE" \
        --ob-priors "$ORIENT_MODEL_DIR/${base_name}_read-orientation-model.tar.gz" \
        -O "$FILTERED_VCF_DIR/${base_name}_filtered.vcf" || error_exit "FilterMutectCalls failed"

    # Optional: VariantAnnotator
    # log_message "Step 5.6: VariantAnnotator - Add standard annotations"
    # $GATK_PATH VariantAnnotator \
    #     --reference "$GENOME" \
    #     -I "$recal_bam" \
    #     --variant "$FILTERED_VCF_DIR/${base_name}_filtered.vcf" \
    #     --output "$ANNOTATED_VCF_DIR/${base_name}_annotated.vcf" \
    #     --enable-all-annotations \
    #     --dbsnp "$DBSNP" || error_exit "VariantAnnotator failed"

    # Step 5.7: Add functional annotations (MAF format)
    log_message "Step 5.7: Funcotator - Add functional annotations (MAF format)"
    $GATK_PATH --java-options "$JAVA_OPTS" Funcotator \
        --variant "$FILTERED_VCF_DIR/${base_name}_filtered.vcf" \
        --reference "$GENOME" \
        --ref-version hg38 \
        --data-sources-path "$FUNCOTATOR_DATASOURCES" \
        --output "$MAF_DIR/${base_name}.maf" \
        --output-file-format MAF \
        -L "$BED_FILE" \
        --remove-filtered-variants true || error_exit "Funcotator failed"

    log_message "Finished GATK processing for: $line"
    log_message "---"

done < "$INPUT_LIST_FILE"

log_message "GATK pipeline steps complete."
```
*[Note: This script is an example. You **must** replace placeholder paths (e.g., `/path/to/...`) with actual paths on your system. Ensure GATK, reference files, and input BAM list are correctly specified. The script includes basic error handling for the input file argument.]*

## Disclaimer and Support


This repository provides supplementary code to accompany the publication mentioned above. **This is not the release of a validated software package.** We are providing this information and code solely in addition to a description of methods for making it easier to reproduce our analyses as presented in the paper.

**We are *not* providing any user support, troubleshooting, or maintenance for these scripts.** Due to resource limitations, we cannot respond to queries regarding the execution of the code, adaptation to other datasets, integration into other workflows, or general bioinformatics assistance. Please refer to the publication for detailed methods.

## Contributing

While this repository primarily serves to document the analysis for our publication, we welcome issue reports for critical errors that prevent reproduction of our results.
## License

This project is licensed under a custom Academic Use License - see the [LICENSE](LICENSE) file for details.

This license allows free use for academic and research purposes only, while prohibiting commercial use without explicit permission. We encourage academic collaboration and only request proper citation when using this work in research or publications.


## Contact / Maintainers

This code was developed by the authors of the publication:

* Rohini Sharma, Sultan N. Alharbi.

For questions about the scientific findings, please refer to the publication or contact the corresponding author listed therein.

**Please remember the "no support" policy outlined in the Disclaimer section regarding the code itself.**
