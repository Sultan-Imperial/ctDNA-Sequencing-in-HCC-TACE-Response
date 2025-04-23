# ctDNA Sequencing in HCC-TACE Response Analysis

[![BioProject](https://img.shields.io/badge/BioProject-PRJNA1199049-brightgreen.svg)](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1199049)
[![License](https://img.shields.io/badge/License-Academic%20Use-yellow.svg)](LICENSE)

## Overview

This repository contains the computational analysis pipeline used in our study of circulating tumor DNA (ctDNA) as a biomarker for clinical outcomes following transarterial chemoembolization (TACE) in patients with hepatocellular carcinoma (HCC).


The primary purpose of this repository is to provide the code used for the analyses presented in the associated publication, facilitating the reproducibility of the study's findings.

## Related Publication

The methods and results using these scripts are detailed in the following publication:

**Title:** Deep sequencing of circulating tumour DNA as a biomarker of clinical outcome to transarterial chemoembolisation in hepatocellular carcinoma.
**Authors:** Rohini Sharma, Sultan N. Alharbi, Ksenia Ellum, Leila Motedayen-Aval, Andrea Casadei-Gardini, David J. Pinato, Dominik Bettinger, Bertram Bengsch, Rishi Patel, Joanne Evans.

## Data Availability

The raw sequencing data generated and analysed during this study have been submitted to the National Center for Biotechnology Information (NCBI) Sequence Read Archive (SRA) and are available under the BioProject accession number: **PRJNA1199049**.

You can access the data here: [https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1199049](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1199049)

[TODO: Specify if any essential intermediate data files needed to run the scripts *directly* are included in this repository. If not, state clearly that only scripts are provided and raw data must be obtained from NCBI and processed.]

## Repository Contents

This repository includes the following main scripts and components:

* `[script_name_1.py/R/sh]`: [Briefly describe what this script does, e.g., Script for aligning raw sequencing reads]
* `[script_name_2.py/R/sh]`: [Briefly describe what this script does, e.g., Script for variant calling and filtering]
* `[script_name_3.py/R/sh]`: [Briefly describe what this script does, e.g., Script for statistical analysis and generating figures]
* `[config_file.yaml/json]` (if applicable): [Describe the purpose of the configuration file]
* `[environment.yml / requirements.txt]` (if applicable): [File listing software dependencies]
* `data/` (if applicable): [Describe any small example, metadata, or annotation files included]

[TODO: List the key scripts, configuration files, and directories in your repository. Provide a concise description for each.]

## Primary data processing of sequencing data

## Getting Started

### Prerequisites

* Java Runtime Environment (for Je & GATK)
* Python [Version, e.g., 3.9+]
* R [Version, e.g., 4.2+]
* UMI-tools (v1.0.0)
* FastQC (v0.11.9)
* BWA (v.0.7.17)
* Je-MarkDuplicates (`org.embl.gbcs.je.jeduplicates.MarkDuplicatesWithMolecularCode`) - [TODO: Specify Je version]
* GATK (Genome Analysis Toolkit, v4.2.6.1 used in example script)
* [Reference Genome Files: GRCh38 FASTA]
* [Known Sites Files: dbSNP VCF, gnomAD VCF, PON VCF]
* [Target Interval BED file]
* [Funcotator Data Sources]
* [Other essential libraries/tools, e.g., pandas, Bioconductor pkgs]
* [TODO: Ensure all required tools and data files are listed]

*Recommendation: Use the provided environment file (`environment.yml` or `requirements.txt`).*
### Installation / Setup

1.  **Clone the repository:**
    ```bash
    git clone [https://docs.github.com/repositories/creating-and-managing-repositories/about-repositories](https://docs.github.com/repositories/creating-and-managing-repositories/about-repositories)
    cd [repository-folder-name]
    ```
    *[Note: The URL `https://docs.github.com/...` in the previous version was incorrect placeholder, replace with your actual repository URL]*
2.  **Set up the environment (if applicable):**
    * Using Conda: `conda env create -f environment.yml`
    * Using pip: `pip install -r requirements.txt`
    * [Add any other manual setup steps required, e.g., compiling tools]
3.  **Download Data:** Obtain the raw sequencing data from NCBI SRA using the accession **PRJNA1199049**. You may need tools like `sra-tools` for this. [Add specific instructions or links if helpful]. Ensure the data is placed in the expected directory structure for the scripts.

[TODO: Adjust these steps based on your actual project setup and how the scripts expect data.]

## Analysis Pipeline Overview

1.  **UMI Processing:** Extract UMIs from trimmed PE reads (`UMI-tools` v1.0.0).
2.  **QC:** Assess quality (`FastQC` v0.11.9).
3.  **Alignment:** Align reads to GRCh38 (`BWA` v.0.7.17).
4.  **Duplicate Marking:** Mark/remove PCR duplicates using UMIs (`Je MarkDuplicatesWithMolecularCode`).
5.  **Base Quality Score Recalibration (BQSR):** Recalibrate base quality scores using known variants (`GATK BaseRecalibrator`, `ApplyBQSR`).
6.  **Variant Calling:** Identify somatic variants (`GATK Mutect2`).
7.  **Variant Filtering & Annotation:** Filter variants based on various metrics and annotate them (`GATK FilterMutectCalls`, `VariantAnnotator`, `Funcotator`).
8.  **BAM Merging:** Merge BAMs per sample across lanes/runs (if applicable, may occur earlier).

## Usage

[TODO: Provide concise, step-by-step execution instructions with example commands for each pipeline stage.]

*Example Workflow Steps:*

1.  **UMI Extraction & QC:** `bash scripts/preprocess_reads.sh ...` [TODO: Add specific example]
2.  **Alignment:** `bash scripts/run_alignment.sh ...` [TODO: Add specific example]
3.  **Mark Duplicates (Je Example):**
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

4.  **Merge BAMs (if needed):** `bash scripts/merge_bams.sh ...` [TODO: Add specific example]

5.  **Downstream Analysis (GATK Example Script):**
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
