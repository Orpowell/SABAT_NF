# SABAT Exon Annotation Pipeline

The SABAT Exon Annotation Pipeline is workflow written using Nextflow for the automatic 
prediction of homologous genes using SABAT. 

## Dependencies

- Nextflow (v23.0.0+)
- SABAT (prepackaged with the pipeline)
- BLAST+ (2.12.0+) (Required for SABAT)
- Python (v3.12.1) (Required for SABAT)
- Biopython (v1.78) (Required for SABAT)
- pandas(v2.1.4) (Required for SABAT)


We recommend creating a conda environment containing all dependencies using the follow commands:

    conda create -n SABAT python=3.12.1
    conda install -c bioconda biopython blast nextflow
    conda install pandas

For devices using apple silicon (M1 or later) add the following parameter when creating the environment:

    CONDA_SUBDIR=osx-64 conda create -n SABAT python=3.12.1

## Usage

    conda activate SABAT_NF

    nextflow run https://github.com/Orpowell/SABAT_NF -r main \
    -resume \
    --input samplesheet.csv \
    --sequence cds.fasta \
    --outdir path/to/output/
    --b2b_exons 4 \
    --b2b_coverage 0.95 \
    --b2b_locus_size 4000 

### Parameters

| Parameter | Description |
| --- | --- |
| --input   | Path to the CSV sample sheet (See: Preparing a Sample Sheet).  |
| --sequence| Path to the FASTA file containing the CDS of the gene. |
| --outdir | Path to the output directory. Note: A individual directory is made for each sample in the output directory. |
| --b2b_exons| SABAT *blast2bed* parameter (See: [SABAT](https://github.com/Orpowell/SABAT)). |
| --b2b_coverage| SABAT *blast2bed* parameter (See: [SABAT](https://github.com/Orpowell/SABAT)). |
| --b2b_locus_size| SABAT *blast2bed* parameter (See: [SABAT](https://github.com/Orpowell/SABAT)). |

## Preparing a Sample Sheet

The pipeline requires a CSV file containing 2 columns (see: test.csv):
1. sample: The name of the sample (used name the output directory).
2. genome: Path to the genome to analysed.

Please note all sample names must be unique or results will be overwritten by the pipeline! 

## Pipeline Outputs

The Pipeline creates a directory within the specified output directory (--outdir) for each sample. Within the directory for each sample will be the BED file, the output from SABAT, and 
a second directory containing the BLASTDB for the sample. All directories and files are named using the sample's name. To save space, the BLASTDB is symbollic link to the work directory generated
by NextFlow.

output
    |
    |-sample_1
    |   |-sample_1.bed
    |   |-sample_1
    |        |-BLASTDB files for sample_1
    |
    |-sample_2
    |    |-sample_2.bed
    |    |-sample_2
    |        |-BLASTDB files for sample_2

## How does the pipeline work?

The pipeline performs the following steps for each sample:
1. Generate a BLASTDB of the genome (makeblastdb)
2. BLAST target gene CDS against the genome BLASTDB (dc-MEGABLAST)
3. Predict identified gene loci from BLAST output (SABAT blast2bed)



