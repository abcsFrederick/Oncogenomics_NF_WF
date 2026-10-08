
include {BWA_picard} from '../subworkflows/Bwa_picard_bam_processing'
include {Exome_GATK} from '../subworkflows/Exome_GATK'
include {QC_exome_bam} from '../subworkflows/QC_exome_bam'
include {HLA_calls_exome} from '../subworkflows/HLA_calls_exome'
include {Kraken
        Fastqc
        Krona
        Fastq_screen} from '../modules/qc/qc'


workflow Exome_common_WF {

kraken_bacteria = Channel.of(file(params.kraken_bacteria, checkIfExists:true))
fastq_screen_config         = Channel.of(file(params.fastq_screen_config, checkIfExists:true))

take:
     samples_exome_ch

main:


BWA_picard(samples_exome_ch)


Fastq_screen_input = samples_exome_ch.map{ meta, r1, r2 -> [meta, [r1,r2]] }
                        .combine(fastq_screen_config)
                        .combine(Channel.fromPath(params.fastq_screen_db,type: 'dir', checkIfExists: true))

Fastq_screen(Fastq_screen_input)

Kraken(samples_exome_ch
    .combine(kraken_bacteria)
)

ch_versions = BWA_picard.out.ch_versions.mix(Kraken.out.versions)

Krona(Kraken.out.kraken_out)

HLA_calls_exome(samples_exome_ch)

ch_versions = ch_versions.mix(HLA_calls_exome.out.ch_versions)

capture_ch = BWA_picard.out.picard_MD
    .map { tuple ->
        def meta = tuple[0]
        def bam = tuple[1]
        def bai = tuple[2]
        def target_file = ''

        if (meta.sc == 'clin.ex.v1') {
            target_file = params.clin_ex_v1
        } else if (meta.sc == 'seqcapez.hu.ex.v3') {
            target_file = params.seqcapez.hu.ex.v3
        } else if (meta.sc == 'seqcapez.rms.v1') {
            target_file = params.seqcapez.rms.v1_target
        } else if (meta.sc == 'agilent.v7') {
            target_file = params.agilent.v7_target
        } else if (meta.sc == 'idt_v2_plus') {
            target_file = params.idt_v2_plus
        } else if (meta.sc == 'xgen-hyb-panelv2') {
            target_file = params.xgen_hyb_panelv2_target
        } else if (meta.sc == 'comp_ex_v1') {
            target_file = params.comp_ex_v1_target
        } else if (meta.sc == 'seqcapez.hu.ex.utr.v1') {
            target_file = params.seqcapez.hu.ex.utr.v1_target
        } else if (meta.sc == 'hybrid_selection') {
            target_file = params.hybrid_selection
        }
        return [meta,target_file]
    }



Exome_GATK(BWA_picard.out.picard_MD,
           capture_ch
)


ch_versions = ch_versions.mix(Exome_GATK.out.ch_versions)

//update design channel
design_ch =  Exome_GATK.out.GATK_Exome_bam
   .map { tuple ->
        def meta = tuple[0]
        def bam = tuple[1]
        def bai = tuple[2]
        def design_file = ''

        if (meta.sc == 'clin.ex.v1') {
            design_file = params.clin_ex_v1_design
        } else if (meta.sc == 'agilent.v7') {
            design_file = params.agilent.v7_design
        } else if (meta.sc == 'idt_v2_plus') {
            design_file = params.idt_v2_plus_design
        } else if (meta.sc == 'seqcapez.rms.v1') {
            design_file = params.seqcapez.rms.v1_design
        } else if (meta.sc == 'xgen-hyb-panelv2') {
            design_file = params.xgen_hyb_panelv2_design
        } else if (meta.sc == 'comp_ex_v1') {
            design_file = params.comp_ex_v1_design
        } else if (meta.sc == 'seqcapez.hu.ex.v3') {
            design_file = params.seqcapez.hu.ex.v3_design
        } else if (meta.sc == 'seqcapez.hu.ex.utr.v1') {
            design_file = params.seqcapez.hu.ex.utr.v1_design
        } else if (meta.sc == 'hybrid_selection') {
            design_file = params.hybrid_selection_design
        }

        return [meta,design_file]
    }


QC_exome_bam(
    Exome_GATK.out.GATK_Exome_bam,
    BWA_picard.out.bwa_bam,
    capture_ch,
    design_ch
)

ch_versions = ch_versions.mix(QC_exome_bam.out.ch_versions)

emit:
coverage = QC_exome_bam.out.coverage
pileup = QC_exome_bam.out.hotspot_pileup
snpeff_vcf =Exome_GATK.out.SnpEff_vcf
exome_final_bam = Exome_GATK.out.GATK_Exome_bam
loh = QC_exome_bam.out.loh
gt = QC_exome_bam.out.gt
target_capture_ch = capture_ch
HC_snpeff_snv_vcf2txt = Exome_GATK.out.HC_snpeff_snv_vcf2txt
mergehla_exome = HLA_calls_exome.out.mergehla_exome
exome_qc = QC_exome_bam.out.Exome_qc
markdup_txt = BWA_picard.out.markdup_txt
hsmetrics = QC_exome_bam.out.hsmetrics
flagstat = QC_exome_bam.out.flagstat
verifybamid = QC_exome_bam.out.verifybamid
Fastqc_out = BWA_picard.out.Fastqc_out
krona = Krona.out
kraken = Kraken.out.kraken_out
ch_versions = ch_versions
design_ch = design_ch
fastq_screen = Fastq_screen.out
conpair_pileup = QC_exome_bam.out.conpair_pileup
hotspot_depth = QC_exome_bam.out.hotspot_depth
}
