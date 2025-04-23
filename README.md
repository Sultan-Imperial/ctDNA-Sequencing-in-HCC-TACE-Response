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

   
3.  **Alignment & Sorting (BWA + Samtools):**
*This script processes a list of sample prefixes (one per line, provided as input `$1`), aligning corresponding paired-end FASTQ files (`${prefix}_R1.fastq.gz`, `${prefix}_R2.fastq.gz`) and creating sorted BAM files.*
    ```bash
    #!/bin/bash
    # --- Configuration (Consider moving these to a config file) ---
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
        FASTQ1="${line}_R1.fastq.gz" # [TODO: Adjust path/naming if needed]
        FASTQ2="${line}_R2.fastq.gz" # [TODO: Adjust path/naming if needed]
        
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


5.  **Mark Duplicates (Je Example):**
    *The following command uses `je markdupes`. It assumes the input BAM filename is stored in a shell variable `$line`. Adjust paths and memory (`-Xmx`) as needed.*
    ```bash
    # Example assuming input file path is in variable $line
    # Adjust java path and memory (-Xmx) if needed
    java -Xmx[Mem] -jar [path/to/je.jar] markdupes \
      INPUT=$line \
      OUTPUT=/path/to/output/MarkDuplicated_Je/${line%.sorted.bam}MarkDuplicated_Je.bam \
      MM=1 SLOTS=-1 ASSUME_SORTED=true \
      METRICS_FILE=/path/to/output/MarkDuplicated_Je/${line%.sorted.bam}MarkDuplicated_Je.txt \
      REMOVE_DUPLICATES=true
      # Add other necessary parameters if used
    ```
    *[Note: Replace `[Mem]` (e.g., `4g`), `[path/to/je.jar]`, and output paths with appropriate values.]*

6.  **Merge BAMs (if needed):** `bash scripts/merge_bams.sh ...` [TODO: Add specific example]

7.  **Downstream Analysis (GATK Example Script):**
    *The following script demonstrates BQSR, Mutect2 variant calling, filtering, and annotation steps using GATK 4.2.6.1. It expects a file containing a list of input BAM file paths as its first argument (`$1`). Adjust paths to tools, reference files, and output directories.*
    ```bash
    #!/bin/bash

    # --- Configuration ---
    GATK_PATH="/path/to/gatk-4.2.6.1/gatk" # Path to GATK executable
    GENOME="/path/to/BWA_GRCh38.fa"
    DBSNP="/path/to/Homo_sapiens_assembly38.dbsnp138.vcf"
    BED_FILE="/path/to/Rohini_target_merged.bed"
    PON="/path/to/1000g_pon.hg38.vcf.gz" # Panel of Normals
    GERMLINE_RESOURCE="/path/to/af-only-gnomad.hg38.vcf.gz"
    FUNCOTATOR_DATASOURCES="/path/to/funcotator_dataSources.v1.7.20200521s"

    # Output Directories (Create if they don't exist)
    RECAL_TABLE_DIR="/path/to/output/Recalibrated_dataTable"
    RECAL_BAM_DIR="/path/to/output/Recalibrated_BAM"
    MUTECT_BAM_DIR="/path/to/output/Final_VCF_Mutect2/OutputBAM"
    MUTECT_VCF_DIR="/path/to/output/Final_VCF_Mutect2/RawData_VCF"
    F1R2_DIR="/path/to/output/Final_VCF_Mutect2/F1R2"
    ORIENT_MODEL_DIR="/path/to/output/Final_VCF_Mutect2/F1R2" # Often same as F1R2
    FILTERED_VCF_DIR="/path/to/output/Final_VCF_Mutect2/OptimizinngmutectFiltering"
    ANNOTATED_VCF_DIR="/path/to/output/Final_VCF_Mutect2/AnnotatedVCF" # Example dir for VariantAnnotator
    MAF_DIR="/path/to/output/Final_VCF_Mutect2/MAF_output" # Example dir for Funcotator

    mkdir -p $RECAL_TABLE_DIR $RECAL_BAM_DIR $MUTECT_BAM_DIR $MUTECT_VCF_DIR $F1R2_DIR $ORIENT_MODEL_DIR $FILTERED_VCF_DIR $ANNOTATED_VCF_DIR $MAF_DIR

    # --- Processing Loop ---
    # Input: File containing one BAM path per line (passed as $1)
    if [ -z "$1" ]; then
      echo "Usage: $0 <file_with_bam_paths>"
      exit 1
    fi

    while IFS= read -r line || [[ -n "$line" ]]; do
        echo "Processing BAM: $line"
        
        # Extract base filename for output naming
        base_name=$(basename "${line%MarkDuplicated_Je.bam}") # Assumes input ends with MarkDuplicated_Je.bam
        recal_base_name=$(basename "${line%_recalibrated.bam}") # For steps after ApplyBQSR if naming changes

        # --- GATK Steps ---
        echo "Step 1: BaseRecalibrator"
        $GATK_PATH BaseRecalibrator \
          -I "$line" \
          --reference "$GENOME" \
          --known-sites "$DBSNP" \
          --intervals "$BED_FILE" \
          --output "$RECAL_TABLE_DIR/${base_name}_recal_data.table"

        echo "Step 2: ApplyBQSR"
        $GATK_PATH ApplyBQSR \
          -I "$line" \
          --reference "$GENOME" \
          --bqsr-recal-file "$RECAL_TABLE_DIR/${base_name}_recal_data.table" \
          --intervals "$BED_FILE" \
          -O "$RECAL_BAM_DIR/${base_name}_recalibrated.bam"
        
        # Define recalibrated BAM path for subsequent steps
        recal_bam="$RECAL_BAM_DIR/${base_name}_recalibrated.bam"
        recal_base_name=$(basename "${recal_bam%_recalibrated.bam}") # Update base name for VCF outputs

        echo "Step 3: Mutect2"
        $GATK_PATH Mutect2 \
          --reference "$GENOME" \
          --panel-of-normals "$PON" \
          --germline-resource "$GERMLINE_RESOURCE" \
          --intervals "$BED_FILE" \
          --input "$recal_bam" \
          --bam-output "$MUTECT_BAM_DIR/${recal_base_name}_output.bam" \
          --output "$MUTECT_VCF_DIR/${recal_base_name}_unfiltered.vcf" \
          --f1r2-tar-gz "$F1R2_DIR/${recal_base_name}_f1r2.tar.gz" \
          --force-active true --initial-tumor-lod 0.0 --tumor-lod-to-emit 0.0 -genotype-pon-sites true

        echo "Step 4: LearnReadOrientationModel"
        $GATK_PATH LearnReadOrientationModel \
          --input "$F1R2_DIR/${recal_base_name}_f1r2.tar.gz" \
          --output "$ORIENT_MODEL_DIR/${recal_base_name}_read-orientation-model.tar.gz"

        echo "Step 5: FilterMutectCalls"
        $GATK_PATH FilterMutectCalls \
          --reference "$GENOME" \
          --intervals "$BED_FILE" \
          --ob-priors "$ORIENT_MODEL_DIR/${recal_base_name}_read-orientation-model.tar.gz" \
          --variant "$MUTECT_VCF_DIR/${recal_base_name}_unfiltered.vcf" \
          --output "$FILTERED_VCF_DIR/${recal_base_name}_filtered.vcf" \
          --lenient true --min-median-base-quality 10 --f-score-beta 1
          # Note: The original example had --input $line here, which seems incorrect; should likely use --variant

        # Optional: VariantAnnotator (Example - adapt input/output as needed)
        # echo "Step 6: VariantAnnotator"
        # $GATK_PATH VariantAnnotator \
        #  --reference "$GENOME" \
        #  -I "$recal_bam" \
        #  --variant "$FILTERED_VCF_DIR/${recal_base_name}_filtered.vcf" \
        #  --output "$ANNOTATED_VCF_DIR/${recal_base_name}_annotated.vcf" \
        #  --enable-all-annotations \
        #  --dbsnp "$DBSNP"

        echo "Step 7: Funcotator (MAF Output)"
        $GATK_PATH Funcotator \
          --reference "$GENOME" \
          --ref-version hg38 \
          --intervals "$BED_FILE" \
          --data-sources-path "$FUNCOTATOR_DATASOURCES" \
          --output "$MAF_DIR/${recal_base_name}.maf" \
          --output-file-format MAF \
          --variant "$FILTERED_VCF_DIR/${recal_base_name}_filtered.vcf" \
          --remove-filtered-variants true
          # Note: The original example had --input $line here, which seems incorrect for Funcotator variant input

        echo "Finished processing: $line"
        echo "---"

    done < "$1"

    echo "All BAM files processed."

    ```
    *[Note: This script is an example. You **must** replace placeholder paths (e.g., `/path/to/...`) with actual paths on your system. Ensure GATK, reference files, and input BAM list are correctly specified. The script includes basic error handling for the input file argument.]*

    **Alternative Mutect2 Examples:**
    *The following examples show alternative ways to run `gatk Mutect2` directly, assuming the input BAM file path is in the `$line` variable (e.g., within a loop reading from a file list passed as `$1`). Adjust paths to GATK, reference, intervals, and output directories.*

    *Example 1: Basic Mutect2 call*
    ```bash
    #!/bin/bash
    GATK_PATH="/path/to/gatk-4.2.6.1/gatk"
    GENOME="/path/to/BWA_GRCh38.fa"
    BED_FILE="/path/to/Rohini_target_merged.bed"
    OUTPUT_DIR="/path/to/output/Mutect2_LastOne"
    mkdir -p $OUTPUT_DIR

    if [ -z "$1" ]; then echo "Usage: $0 <file_with_bam_paths>"; exit 1; fi

    while IFS= read -r line || [[ -n "$line" ]]; do
        echo "Processing BAM (Basic Mutect2): $line"
        base_name=$(basename "${line%.marked_duplicates.bam}") # Assumes input ends differently here
        $GATK_PATH Mutect2 \
          --reference "$GENOME" \
          --input "$line" \
          --output "$OUTPUT_DIR/${base_name}.unfiltered.vcf" \
          --intervals "$BED_FILE"
    done < "$1"
    ```

    *Example 2: Mutect2 with enhanced annotations*
    ```bash
    #!/bin/bash
    GATK_PATH="/path/to/gatk-4.2.6.1/gatk"
    GENOME="/path/to/BWA_GRCh38.fa"
    BED_FILE="/path/to/Rohini_target_merged.bed"
    OUTPUT_DIR="/path/to/output/VCF_Mutect2"
    mkdir -p $OUTPUT_DIR

    if [ -z "$1" ]; then echo "Usage: $0 <file_with_bam_paths>"; exit 1; fi

    while IFS= read -r line || [[ -n "$line" ]]; do
        echo "Processing BAM (Annotated Mutect2): $line"
        base_name=$(basename "${line%.marked_duplicates.bam}") # Assumes input ends differently here
        $GATK_PATH Mutect2 \
          --input "$line" \
          --output "$OUTPUT_DIR/${base_name}.vcf" \
          --reference "$GENOME" \
          --annotation Coverage \
          --f1r2-max-depth 200000000 \
          --genotype-germline-sites true \
          --intervals "$BED_FILE" \
          --max-reads-per-alignment-start 0 \
          --enable-all-annotations true
    done < "$1"
    ```

## Disclaimer and Support


This repository provides supplementary code to accompany the publication mentioned above. **This is not the release of a validated software package.** We are providing this information and code solely in addition to a description of methods for making it easier to reproduce our analyses as presented in the paper.

**We are *not* providing any user support, troubleshooting, or maintenance for these scripts.** Due to resource limitations, we cannot respond to queries regarding the execution of the code, adaptation to other datasets, integration into other workflows, or general bioinformatics assistance. Please refer to the publication for detailed methods.

## Contributing

While this repository primarily serves to document the analysis for our publication, we welcome issue reports for critical errors that prevent reproduction of our results.
## License

[TODO: Choose an open-source license (e.g., MIT, Apache 2.0, GPLv3) for your code. Create a `LICENSE` file in your repository containing the full license text and state the chosen license here.]

*Example: This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.*

*Note: If no license is specified, default copyright laws apply, restricting reuse.*

## Contact / Maintainers

This code was developed by the authors of the publication:

* Rohini Sharma, Sultan N. Alharbi, Ksenia Ellum, Leila Motedayen-Aval, Andrea Casadei-Gardini, David J. Pinato, Dominik Bettinger, Bertram Bengsch, Rishi Patel, Joanne Evans.

For questions about the scientific findings, please refer to the publication or contact the corresponding author listed therein.

**Please remember the "no support" policy outlined in the Disclaimer section regarding the code itself.**
