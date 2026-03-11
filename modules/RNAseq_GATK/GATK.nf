process GATK_RNASeq_Trim {
        tag "$meta.lib"

        input:
        tuple val(meta),
        path(bam),
        path(index),
        path(genome),
        path(genome_fai),
        path(genome_dict)

        output:
        tuple val(meta),
        path("${meta.lib}.trim.bam"),
        path("${meta.lib}.trim.bai")


        script:

        """
        java -jar -Xmx40g \$GATK_JAR -T SplitNCigarReads -R $genome -I $bam -o ${meta.lib}.trim.bam -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS

        """
}


process GATK_RTC_IR {
        tag "$meta.lib"

        input:
	tuple val(meta),
        path(bam),
        path(index),
        path(genome),
        path(genome_fai),
        path(genome_dict),
        path(phase1_1000g),
        path(Mills_and_1000g)

        output:
        tuple val(meta),path("${meta.lib}.Ir.bam"), emit: Ir_bam
        tuple val(meta),path("${meta.lib}.Ir.bai"), emit: Ir_bai
        path "versions.yml"             , emit: versions


        script:
        """

	java -jar -Xmx80g \$GATK_JAR -T RealignerTargetCreator -nt 10 -R $genome -known $phase1_1000g -known $Mills_and_1000g -I $bam -o ${meta.lib}.realignment.intervals

        java -jar -Xmx80g \$GATK_JAR -T IndelRealigner -R $genome -known $phase1_1000g -known $Mills_and_1000g -I $bam --targetIntervals ${meta.lib}.realignment.intervals -o ${meta.lib}.Ir.bam --maxReadsInMemory 1500000

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            GATK: \$(java -jar \$GATK_JAR -T RealignerTargetCreator --version)
        END_VERSIONS
        """
}



process GATK_BR_PR {
        tag "$meta.lib"

        publishDir "${params.resultsdir}/${meta.id}/${meta.casename}/${meta.lib}", mode: "${params.publishDirMode}"

        input:
	tuple val(meta),
        path(bam),
        path(index),
        path(genome),
        path(genome_fai),
        path(genome_dict),
        path(phase1_1000g),
        path(Mills_and_1000g)

        output:
        tuple val(meta),path("${meta.lib}.final.bam"), emit: final_bam
        tuple val(meta),path("${meta.lib}.final.bam.bai"), emit: final_bai
        path "versions.yml"             , emit: versions


        script:
        def prefix = task.ext.prefix ?: "${meta.lib}"
        """

        java -jar -Xmx40g \$GATK_JAR -T BaseRecalibrator -R $genome -knownSites $phase1_1000g -knownSites $Mills_and_1000g -I $bam -o ${prefix}.recalibration.matrix.txt

        java -jar -Xmx40g \$GATK_JAR -T PrintReads -R $genome -I $bam -o ${prefix}.final.bam -BQSR ${prefix}.recalibration.matrix.txt

        mv ${prefix}.final.bai ${prefix}.final.bam.bai

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            GATK: \$(java -jar \$GATK_JAR -T PrintReads --version)
        END_VERSIONS
        """
}


process RNAseq_HaplotypeCaller {

     tag "$meta.lib"

     publishDir "${params.resultsdir}/${meta.id}/${meta.casename}/${meta.lib}/calls", mode: "${params.publishDirMode}"

     input:
     tuple val(meta),
        path(bam),
        path(index),
        path(genome),
        path(genome_fai),
        path(genome_dict),
        path(dbsnp)

     output:
     tuple val(meta),path("${meta.lib}.HC_RNASeq.raw.vcf"), emit: rna_HC
     path "versions.yml"             , emit: versions

     stub:
     """
     touch "${meta.lib}.HC_RNASeq.raw.vcf"
     """

     script:
     def prefix = task.ext.prefix ?: "${meta.lib}"
     """
     set -exo pipefail
     java -jar -Xmx40g \$GATK_JAR -T HaplotypeCaller -R ${genome} -I ${bam} -o ${prefix}.vcf --dbsnp ${dbsnp} -dontUseSoftClippedBases -stand_call_conf 30 -nct ${task.cpus}
     java -jar -Xmx40g \$GATK_JAR -T VariantFiltration -R ${genome} -V ${prefix}.vcf -window 35 -cluster 3 --filterExpression "FS > 30.0 || QD < 2.0" -filterName "RNASeqFilters_FS_QD" -o ${prefix}.HC_RNASeq.raw.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        GATK: \$(java -jar \$GATK_JAR -T HaplotypeCaller --version)
    END_VERSIONS

     """
}

process Exome_HaplotypeCaller {

     tag "$meta.lib"

     publishDir "${params.resultsdir}/${meta.id}/${meta.casename}/${meta.lib}/calls", mode: "${params.publishDirMode}"

     input:
     tuple val(meta),
        path(bam),
        path(index),
        path(Capture_bed),
        path(genome),
        path(genome_fai),
        path(genome_dict),
        path(dbsnp)


     output:
     tuple val(meta),path("${meta.lib}.HC_DNASeq.raw.vcf"), emit: exome_HC
     path "versions.yml"             , emit: versions

     stub:
     """
     touch "${meta.lib}.HC_DNASeq.raw.vcf"
     """

     script:
     def prefix = task.ext.prefix ?: "${meta.lib}"
     """
     set -exo pipefail
     java -jar -Xmx100g \$GATK_JAR -T HaplotypeCaller -R ${genome} -I ${bam} -o ${prefix}.HC_DNASeq.raw.vcf --dbsnp ${dbsnp} -L ${Capture_bed} -mbq 20 -mmq 30 -nct ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        GATK: \$(java -jar \$GATK_JAR -T HaplotypeCaller --version)
    END_VERSIONS
     """
}
