workflow{
    input.view()

    blastdb = MakeBlASTDB(input)
    blastdb.view()

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

    output:

    script:
    """

    """
}

process SABAT_blast2bed{
    input:

    output:

    script:
    """

    """
}