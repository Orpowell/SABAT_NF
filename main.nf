workflow{
    input | MakeBlASTDB | dc_megaBLAST | SABAT_blast2bed

    
}

input = Channel.fromPath(params.input)
               .splitCsv(header:true, sep:",")
               .map{row -> [row.sample, file("${row.genome}")]}

process MakeBlASTDB {

    publishDir "${params.outdir}/${sample}"

    maxForks 5

    input:
    tuple val(sample), path(genome)

    output:

    tuple val(sample), path("${sample}")

    script:
    """
    mkdir -p $sample

    makeblastdb \
    -in $genome \
    -parse_seqids \
    -blastdb_version 5 \
    -title TA10171_blast_db \
    -dbtype nucl \
    -out $sample/$sample  

    """
}

process dc_megaBLAST {

    maxForks 10

    input:

    tuple val(sample), path(blastdb)

    output:

    tuple val(sample), path("${sample}.blastn")

    script:
    """
    blastn \
    -task dc-megablast \
    -query $params.sequence \
    -db $blastdb/$sample \
    -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen" \
    -out ${sample}.blastn \
    -num_threads 8 \
    -template_type coding \
    -template_length 16 
    """
}

process SABAT_blast2bed{
    publishDir "${params.outdir}/${sample}", mode: 'copy'

    input:

        tuple val(sample), path(blastn)

    output:

        path "${sample}.bed"

    script:
    """
    __main__.py blast2bed -i $blastn -e $params.b2b_exons -c $params.b2b_coverage -l $params.b2b_locus_size > ${sample}.bed
    """
}