process RECONCNV {

    tag "$meta.lib"
    publishDir "${params.resultsdir}/${meta.id}/${meta.casename}/${meta.lib}/cnvkit", mode: "${params.publishDirMode}", pattern: "${meta.lib}*"

    input:
    tuple val(meta), path(cnr), path(cns)
    path recon_config_file
    path recon_data_dir

    output:
    tuple val(meta), path("${meta.lib}*.html"), emit: html
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def bn = task.ext.prefix ?: "${meta.lib}"

    """
    in=${cnr}
    cns=${cns}


    if [ ! -f ${bn}.html ];then
    	cleaned_cnr=${bn}.nochrM_Y.cnr
    	grep -v '^chrM' \$in | grep -v '^chrY' > \${cleaned_cnr}
    	reconCNV.py -r \${cleaned_cnr} -g \${cleaned_cnr} -x ${recon_data_dir}/hg19_genome_length_chr.txt -a ${recon_data_dir}/hg19_COSMIC_genes_model_chr.txt -s \$cns -c ${recon_config_file} -o ${bn}.html -d ./
    	chrs=()
    	for chr in {1..22};do chrs+=("chr\$chr");done
    	chrs+=('chrX')
    	for chr in "\${chrs[@]}";do
    		chr_cnr=${bn}.\${chr}.cnr
    		chr_cns=${bn}.\${chr}.cns
        chr_len="${recon_data_dir}/hg19_genome_length_\${chr}.txt"

        genome_len=${recon_data_dir}/hg19_genome_length_chr.txt

    		awk -v chr="\$chr" '\$1==chr || \$1=="chromosome"' \$in > \${chr_cnr}
    		awk -v chr="\$chr" '\$1==chr || \$1=="chromosome"' \$cns > \${chr_cns}
    		awk -v chr="\$chr" '\$1==chr || \$1=="chromosome"' \${genome_len} > \${chr_len}
            reconCNV.py -r \${chr_cnr} -g \${chr_cnr} -x \${chr_len} -a ${recon_data_dir}/annotation/hg19/\${chr}.annotation.txt -s \${chr_cns} -c ${recon_config_file} -o ${bn}.\${chr}.html -d ./
    	done
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(reconCNV.py --version 2>&1 | sed 's/Python //g')
    END_VERSIONS
    """
}
