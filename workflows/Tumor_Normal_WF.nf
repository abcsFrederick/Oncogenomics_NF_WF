include {Exome_common_WF} from './Exome_common_WF.nf'
include {Mutect} from '../modules/Variant_analysis/Mutect'
include {Mutect_order} from '../modules/Variant_analysis/Mutect'

workflow Tumor_Normal_WF {

    genome                  = Channel.of(file(params.genome, checkIfExists:true))
    dbsnp_138_b37_vcf       = Channel.of(file(params.dbsnp, checkIfExists:true))
    cosmic_v67_hg19_vcf     = Channel.of(file(params.cosmic_v67_hg19_vcf, checkIfExists:true))
    genome_fai              = Channel.of(file(params.genome_fai, checkIfExists:true))
    genome_dict             = Channel.of(file(params.genome_dict, checkIfExists:true))


samples_exome = Channel.fromPath("Tumor_Normal.csv")
.splitCsv(header:true)
.filter { row -> row.type == "Tumor" || row.type == "Normal" }
.map { row ->
    def meta = [:]
    meta.id    =  row.sample
    meta.lib   =  row.library
    meta.sc    =  row.sample_captures
    meta.casename  = row.casename 
    meta.type     = row.type
    def fastq_meta = []
    fastq_meta = [ meta,  file(row.read1), file(row.read2)  ]

    return fastq_meta
}
//samples_exome.view()

Exome_common_WF(samples_exome)

//Exome_common_WF.out.exome_final_bam.view()
bam_target_ch = Exome_common_WF.out.exome_final_bam.combine(Exome_common_WF.out.target_capture_ch,by:[0])
//bam_target_ch.view()
tumor_bam_channel = bam_target_ch.branch { 
    Tumor: it[0].type == "Tumor"
    Normal: it[0].type == "Normal"
}
//tumor_bam_channel.Normal.view()

Mutect(
   tumor_bam_channel.Tumor,
   tumor_bam_channel.Normal.map { tuple -> tuple.take(tuple.size() - 1) },
   genome,
   genome_fai,
   genome_dict,
   dbsnp_138_b37_vcf,
   cosmic_v67_hg19_vcf
)

Mutect_order(Mutect.out.mutect_raw_vcf)

}