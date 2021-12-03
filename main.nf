#!/usr/bin/env nextflow
/*
Simply run nextflow pipeline as:
nextflow run main.nf --input_bed input_bed --reference_fasta tair10.fasta --db_genomes *fasta --outdir output_folder
*/
/*
 * SET UP CONFIGURATION VARIABLES
*/
params.input_bed = false
params.reference_fasta = false
params.db_genomes = false


/// blast parameters
params.perc_identity = 80




ch_input_bed = Channel.fromPath ( params.input_bed )
    // .splitCsv( header: false, sep: "\t" )
    // .map{ row -> [ "${row[0]}", "${row[1]}", "${row[2]}"  ]  }

ch_ref_fasta = Channel.fromPath ( params.reference_fasta )

ch_db_genomes = Channel
    .fromPath ( params.db_genomes )
    .map{ ["$it.baseName", "${it.getParent()}", file("$it")] }
// need to do makeblastdb for each of these 

process createBlastDB {
    tag { "$db_id" }
    storeDir "$db_parent_dir"

    input:
    set val(db_id), val(db_parent_dir), file(db_genome) from ch_db_genomes

    output:
    set val(db_id), file(db_genome), file("${db_id}.*") into ch_blast_db

    script:
    """
    makeblastdb -dbtype nucl -in ${db_genome} -out ${db_id}
    """
}

process getSequence {
    // tag { "ref_${chr},${start},${end}" }
    // publishDir "", mode: 'symlink'

    input:
    file query_bed from ch_input_bed
    // set val(chr), val(start), val(end) from ch_input_bed
    file ref_fasta from ch_ref_fasta.collect()

    output:
    // file "ref_${chr},${start},${end}.fasta" into ch_query_fasta
    file "query_sequences.fasta" into ch_query_fasta

    script:
    """
    bedtools getfasta -fi $ref_fasta -bed $query_bed -name -fo query_sequences.fasta
    """
    // ref_${chr},${start},${end}.fasta
    // cat '$chr\t$start\t$end\tref_$chr,$start,$end\n' > req_region.bed
}

// combine query fasta to all the genomes
// ch_each_fasta = ch_query_fasta.splitFasta(  )
// ch_inputs_query_genomes = ch_blast_db.combine( ch_query_fasta )

process doBlast {
    tag { "ref--$db_genome" }
    publishDir "$params.outdir", mode: 'copy'

    input:
    set val(db_id), file(db_genome), file(db_genome_blasts) from ch_blast_db
    file(query_fasta) from ch_query_fasta.collect()

    output:
    file("query_${db_id}.bed") into ch_output_hits

    script:
    // query_id = $query_fasta.baseName
    """
    blastn -query $query_fasta \
    -db $db_id \
    -task blastn \
    -evalue 1 \
    -out query_${db_id}.raw.bed \
    -perc_identity $params.perc_identity \
    -outfmt "7 qseqid qacc sacc evalue qstart qend sseqid sstart send" \
    -penalty -2 \
    -reward 1 \
    -word_size 28 \
    -gapopen 1 \
    -gapextend 1
    echo 'query_id\tquery_acc\tsubject_acc\tevalue\tquery_start\tquery_end\tsub_id\tsub_start\tsub_end' > query_${db_id}.bed
    grep -v "^#" query_${db_id}.raw.bed >> query_${db_id}.bed
    """
}