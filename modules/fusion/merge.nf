process Mergefusion {
    tag { dataset_id }

    publishDir "${params.resultsdir}/${dataset_id}/${params.casename}/Actionable", mode: "${params.publishDirMode}"

    input:
    tuple val(dataset_id),
        path(arriba), 
        path(FC),
        path(FC_summary),
        path(SF)
    
    output:
    tuple val("${dataset_id.split('_')[0]}"),
        path("${dataset_id.split('_')[0]}.actionable.fusion.txt")

    stub:
    """
    touch "${dataset_id.split('_')[0]}.actionable.fusion.txt"
    """

    shell:
    '''
    ActionableFusion.v1.pl  !{dataset_id} !{FC} !{SF} !{arriba} $PWD | awk 'NR<2{print $0;next}{print $0| "sort "}' > !{dataset_id.split('_')[0]}.actionable.fusion.txt

    '''
}
