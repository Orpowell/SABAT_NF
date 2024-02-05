workflow{
    input.view()

    blastdb = MakeBlASTDB(input)
    blastdb.view()
    blasthit = dc_megaBLAST(blastdb)
    blasthit.view()
    
}

input = Channel.fromPath(params.input)
               .splitCsv(header:true, sep:",")
               .map{row -> [row.sample, file("${row.genome}")]}

process MakeBlASTDB {
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
    -num_threads 100 \
    -template_type coding_and_optimal \
    -template_length 16 
    """
}

process SABAT_blast2bed{
    input:

    output:

    script:
    """

    """
}