process CNVkitPooled {

 tag "$meta.lib"
    publishDir "${params.resultsdir}/${meta.id}/${meta.casename}/${meta.lib}/cnvkit", mode: "${params.publishDirMode}"

  input:
    tuple val(meta),
        path(bam),
        path(index),
        path(cnv_ref)
 
  output:
    tuple val(meta),path("${meta.lib}.cns"), emit: cnvkit_cns
    tuple val(meta),path("${meta.lib}.cnr"), emit: cnvkit_cnr
    tuple val(meta),path("${meta.lib}.pdf"), emit: cnvkit_pdf
       //path("${meta.lib}.png")
  stub:
     """
     touch "${meta.lib}.cns"
     touch "${meta.lib}.cnr"
     touch "${meta.lib}.pdf"
     
     """
  script:
     def prefix = task.ext.prefix ?: "${meta.lib}"
   """
   cnvkit.py batch -p ${task.cpus} ${bam} --reference ${cnv_ref} --output-dir .
   mv ${prefix}.final.cns ${prefix}.cns
   mv ${prefix}.final.cnr ${prefix}.cnr
   cnvkit.py scatter -s ${prefix}.cn{s,r} -o ${prefix}.pdf
   ##convert -density 150 ${prefix}.pdf ${prefix}.png   --- install imagemagick to the docker
   """
}

process CNVkit_png {

tag "$meta.lib"
    publishDir "${params.resultsdir}/${meta.id}/${meta.casename}/${meta.lib}/cnvkit", mode: "${params.publishDirMode}"

input:

    tuple val(meta),path(pdf)

output:
    tuple val(meta),path("${meta.lib}.png")

stub:
   """
   touch "${meta.lib}.png"
   """

script:
def prefix = task.ext.prefix ?: "${meta.lib}"
"""
pdftoppm -png ${pdf} -singlefile  ${prefix}  

"""

}