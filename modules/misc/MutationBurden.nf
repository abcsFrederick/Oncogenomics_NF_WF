process MutationBurden {

    tag "$meta.id"

    publishDir "${params.resultsdir}/${meta.id}/${meta.casename}/${meta.lib}/qc", mode: "${params.publishDirMode}"

    input:
    tuple val(meta),path(Mutect_annotationfull),
    path(strelka_indels_annotationfull),
    path(strelka_snvs_annotationfull),
    path(tumor_target_capture),
    val(normal),val(tumor),val(vaf),
    val(mutect_ch),
    val(strelka_indelch),
    val(strelka_snvsch)

    output:
    tuple val(meta),path("${meta.lib}.${mutect_ch}.mutationburden.txt"),       emit: mutect
    tuple val(meta),path("${meta.lib}.${strelka_indelch}.mutationburden.txt"), emit: strelka_indels
    tuple val(meta),path("${meta.lib}.${strelka_snvsch}.mutationburden.txt"),  emit: strelka_snvs


    script:
    """
    mutationBurden.py  ${strelka_indels_annotationfull} ${tumor_target_capture} ${tumor} ${normal} ${vaf} > ${meta.lib}.${strelka_indelch}.mutationburden.txt
    mutationBurden.py  ${strelka_snvs_annotationfull} ${tumor_target_capture} ${tumor} ${normal} ${vaf} > ${meta.lib}.${strelka_snvsch}.mutationburden.txt
    mutationBurden.py  ${Mutect_annotationfull} ${tumor_target_capture} ${tumor} ${normal} ${vaf} > ${meta.lib}.${mutect_ch}.mutationburden.txt
    """

}

process combineTMB {

    tag "$meta.id"

    publishDir "${params.resultsdir}/${meta.id}/${meta.casename}/${meta.lib}/qc", mode: "${params.publishDirMode}"

    input:
    tuple val(meta),
    path(mutect_mutationburden),
    path(strelka_indels_mutationburden)

    output:
    tuple val(meta),path("${meta.lib}.combined.mutationburden.txt")

    script:
    """
    combineTMB_oc.py ${mutect_mutationburden} ${strelka_indels_mutationburden} > ${meta.lib}.combined.mutationburden.txt
    """

}
