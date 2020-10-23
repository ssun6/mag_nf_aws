#!/usr/bin/env nextflow

params.reads = 's3://genomics-workflows-c4d2c590-123b-11eb-a1d7-066f8b7a8772/mag_project1/rawData/*_R{1,2}*.fastq.gz'
params.kneaddatabase = 's3://genomics-workflows-c4d2c590-123b-11eb-a1d7-066f8b7a8772/mag_project1/database/knead_database/'
params.outdir = '/Users/shansun/git/mag_nf_aws/results'

log.info """\
         M A G - N F   P I P E L I N E    
         ===================================
         reads        : ${params.reads}
         kneaddatabase    : ${params.kneaddatabase}
         outdir       : ${params.outdir}
         """
         .stripIndent()


/* 
 * Kneaddata
 */
Channel 
    .fromFilePairs( params.reads, checkIfExists:true )
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
    .set { read_pairs_ch } 
    
process kneaddata {
    
    input:
    tuple val(pair_id), path(reads) from read_pairs_ch
    path kneaddatabase from params.kneaddatabase
     
    output:
    tuple val(pair_id), path ('*kneaddata') into filtered_reads_ch

    script:       
    """
    kneaddata --input ${reads[0]} --input ${reads[1]} -db $kneaddatabase -t 10 -p 10 --max-memory 80g -v --output ${pair_id}_kneaddata --bypass-trf --trimmomatic /usr/local/bin/Trimmomatic-0.39
    gzip ${pair_id}_kneaddata/*_paired_1.fastq
    gzip ${pair_id}_kneaddata/*_paired_2.fastq
    """
}

/* 
 * Assembly with metaspades 
 */

process metaspades {
    
    input:  
    tuple val(pair_id), path(kneaddata) from filtered_reads_ch
     
    output:
    tuple val(pair_id), path ('*contigs') into contigs_ch

    script:       
    """
    metaspades.py -1 $kneaddata/*_paired_1.fastq.gz  -2 $kneaddata/*_paired_2.fastq.gz -o ${pair_id}_contigs -m 150
    """
}

/*
 * Bin with metabat2
 */
process metabat2 {
    
    publishDir params.outdir, mode:'copy'
    
    input:
    tuple val(pair_id), path ('contigs') from contigs_ch
 
    output:
    tuple val(pair_id), path ('*bins') into bins_ch
 
    script:
    """
    metabat2 -i $contigs/contigs.fasta -o ${pair_id}_bins/${pair_id}_bin -v -m 2000
    """
}

/*
 * Run checkM to check quality of bins
 */
process checkM {

    publishDir params.outdir, mode:'copy'
    
    input:
    tuple val(pair_id), path ('bins') from bins_ch

    output:
    tuple val(pair_id), path ('*CheckM.txt')
    tuple val(pair_id), path ('*SCG')

    script:
    """
    checkm lineage_wf -f ${pair_id}_CheckM.txt -x fa -t 10 $bins/ ${pair_id}_SCG
    """  
}
