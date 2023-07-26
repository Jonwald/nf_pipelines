#!/usr/bin/env nextflow

params.reads = "fastq/*_R{1,2}.fastq.gz"
params.genome_index = "refs/bbmap_reduced/ref/"
ref_genome = file(params.genome_index)

Channel
    .fromFilePairs( params.reads )
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
    .set { read_pairs }

read_pairs.into{fastqc_input; rrna_input; star_input}


process FASTQC {

        publishDir "fastqc"

        input:
        set pair_id, file(reads) from fastqc_input

        output:
        path "*.html"

        script:
        """
        fastqc ${pair_id}*.fastq.gz
        """
}


process BWAMEM {
        publishDir "bwamem"

        input:
        set pair_id, file(reads) from rrna_input

        output:
        path "*.txt"

        script:
        """
        /data/jyoung/rnaseq_nf/bin/bwa-0.7.17/./bwa mem -t 8 -M /data/jyoung/rnaseq_nf/refs/human_rrna/genome.fa ${reads[0]} ${reads[1]} | samtools view -u -S - | samtools sort -m 256M -@ 8 - | samtools flagstat - > ${pair_id}_rrna_stats.txt

        """
}


process STAR {

        publishDir "star"

        input:
        set pair_id, file(reads) from star_input

        output:
        set pair_id, file("*.bam") into star_output_a
        set pair_id, file("*.bam") into star_output_b

        script:
        """
        /data/jyoung/rnaseq_nf/bin/STAR-2.6.1a/source/./STAR --runThreadN 8 --genomeDir \
        /data/jyoung/rnaseq_nf/refs/Homo_sapiens.GRCh38.dna.primary_assembly.ERCC.Homo_sapiens.GRCh38.91  \
        --readFilesIn ${reads[0]} ${reads[1]} \
        --readFilesCommand zcat \
        --outFilterMultimapNmax 1 \
        --outFilterType BySJout \
        --outSAMattributes NH HI AS NM MD \
        --outFilterMismatchNmax 3 \
        --outFilterMismatchNoverReadLmax 0.04 \
        --alignSJoverhangMin 8 \
        --alignSJDBoverhangMin 1 \
        --outFilterIntronMotifs RemoveNoncanonical \
        --genomeLoad NoSharedMemory \
        --outSAMtype BAM SortedByCoordinate \
        --quantMode GeneCounts \
        --outFileNamePrefix ${pair_id}_
        """
}



/*picard markduplicates */

process markduplicates {

        publishDir "markduplicates"

        input:
        set pair_id, file(bam) from star_output_a

        output:
        set pair_id, file("*.bam") into markduplicates_output

        script:
        """
        java -jar bin/picard.jar MarkDuplicates I=${bam[0]} \
        O=${pair_id}_Aligned.sortedByCoord.out.marked.bam \
        M=${pair_id}_marked_dup_metrics.txt \
        REMOVE_DUPLICATES=False \

        """

}

process removeduplicates {

        publishDir "removeduplicates"

        input:
        set pair_id, file(bam) from star_output_b

        output:
        set pair_id, file("*.bam") into removeduplicates_output

        script:
        """
        java -jar bin/picard.jar MarkDuplicates I=${bam[0]} \
        O=${pair_id}_Aligned.sortedByCoord.out.removed.bam \
        M=${pair_id}_marked_dup_metrics.txt \
        REMOVE_DUPLICATES=True \

        """
}

/*stringtie*/

process stringtie {

        publishDir "stringtie"

        input:
        set pair_id, file(bam) from markduplicates_output

        output:
        set pair_id, file("*.gtf") into stringtie_all_output

        script:
        """
        stringtie -e -B -p 4 --rf -c 0.001 -G refs/gencode.v29.annotation.gtf -o ${pair_id}_stringtie.gtf ${bam[0]}

        """
}

process stringtie_dedup {

        publishDir "stringtie"

        input:
        set pair_id, file(bam) from removeduplicates_output

        output:
        set pair_id, file("*.gtf") into stringtie_dedup_output

        script:
        """
        stringtie -e -B -p 4 --rf -c 0.001 -G refs/gencode.v29.annotation.gtf -o ${pair_id}_stringtie.gtf ${bam[0]}

        """
}

/*generate gene / count matrix*/

process generate_counts {

        publishDir "matrices"

        input:
        set pair_id, file(reads) from stringtie_all_output

        output:
        set pair_id, file("*.txt") into generate_matrix_output

        script:
        """
        python scripts/generate_matrix.py ${pair_id}_featurecounts.txt

        """

}

/*regular QC script*/

process qc_script {

        publishDir "qc_script"

        input:
        set pair_id, file(reads) from generate_matrix_output

        output:
        set pair_id, file("*.txt") into qc_script_output

        script:
        """
        python scripts/QC_script.py ${pair_id}_featurecounts.txt

        """

}


/*stringtie 1 and 2*/

/*generate gene / count matrix*/

/*regular QC script*/

/*expanded QC*/





