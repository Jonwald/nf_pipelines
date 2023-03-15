#!/usr/bin/env nextflow

params.reads = "data/final/*_R{1,2}_001.fastq.gz"

params.adapters = "refs/adapters.fa"
adapters = file(params.adapters)

params.genome_index = "refs/bbmap_reduced/ref/"
ref_genome = file(params.genome_index)

Channel
    .fromFilePairs( params.reads )
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
    .set { read_pairs }


read_pairs.into{trim_input}

/* process FASTQC {

        publishDir "results/pl10/fastqc"

        input:
        set pair_id, file(reads) from fastqc_input

        output:
        path "*.html"

        script:
        """
        fastqc ${pair_id}*.fastq.gz
        """

}

*/
process BBDUK {

        publishDir "results/pl10/trim"

        input:
        set pair_id, file(reads) from trim_input
        file "adapters.fa" from adapters

        output:
        set val(pair_id), file('*_R{1,2}_trimmed.fastq.gz') into mapping_input

        script:
        """
        bbduk.sh in1=${reads[0]} \
        in2=${reads[1]} \
        out1="${pair_id}_R1_trimmed.fastq.gz" \
        out2="${pair_id}_R2_trimmed.fastq.gz" \
        ref=adapters.fa \
        minlength=30 \
        qtrim=t \
        trimq=10 \
        outs="${pair_id}_singletons.fastq.gz" \
        ktrim=r \
        k=19 \
        tbo=t
    """
}

process BBMAP {

        publishDir "results/pl10/bam"

        input:
        set pair_id, file(reads) from mapping_input
        path ref_genome

        output:
        set val(pair_id), file("*.bam") into bam_file

        script:
        """
        bbmap.sh in1=${reads[0]} in2=${reads[1]} out=${pair_id}.bam trimreaddescriptions=t

        """
}

process BCFTOOLS {

        publishDir "results/pl10/vcf"

        input:
        set pair_id, file(bam) from bam_file

        output:
        path "*.vcf"

        script:
        """
        samtools sort ${bam[0]} > ${pair_id}".sorted.bam"
        samtools index ${pair_id}".sorted.bam"
        bcftools mpileup --redo-BAQ --min-BQ 10 --per-sample-mF -a DP,AD -d 100000 -f ../../../refs/il33_reduced_ref_fixed.fasta ${pair_id}".sorted.bam" |  bcftools call --multiallelic-caller -Ov > ${pair_id}".vcf"

        """
}