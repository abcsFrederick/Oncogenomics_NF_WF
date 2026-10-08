
include {Read_depth} from  '../modules/qc/plots'
include {Bam2tdf} from '../modules/qc/plots'
include {Flagstat} from '../modules/qc/plots'
include {Bamutil} from '../modules/qc/plots'
include {HotspotPileup} from '../modules/qc/plots'
include {Genotyping} from '../modules/qc/qc'
include {CircosPlot_lib} from '../modules/qc/qc'
include {VerifyBamID} from '../modules/qc/plots'
include {FailedExons_Genes} from '../modules/qc/plots'
include {Coverage} from  '../modules/qc/plots'
include {TargetIntervals} from  '../modules/qc/plots'
include {HSMetrics} from  '../modules/qc/plots'
include {Conpair_pile} from  '../modules/qc/qc'
include {Exome_QC} from '../modules/qc/qc.nf'
include {Hotspot_Coverage} from '../modules/qc/plots'

workflow QC_exome_bam {

    genome                  = Channel.of(file(params.genome, checkIfExists:true))
    genome_fai              = Channel.of(file(params.genome_fai, checkIfExists:true))
    genome_dict             = Channel.of(file(params.genome_dict, checkIfExists:true))
    chrom_sizes             = Channel.of(file(params.chrom_sizes, checkIfExists:true))
    hg19_hotspot_pos        = Channel.of(file(params.hg19_hotspot_pos, checkIfExists:true))
    clin_ex_v1              = Channel.of(file(params.clin_ex_v1, checkIfExists:true))
    Sites1000g4genotyping   = Channel.of(file(params.Sites1000g4genotyping, checkIfExists:true))
    recode_vcf              = Channel.of(file(params.recode_vcf, checkIfExists:true))
    conpair_refbed          = Channel.of(file(params.conpair_refbed, checkIfExists:true))
    access_hotspot           = Channel.of(file(params.access_hotspot, checkIfExists:true))
    sorted_chr_order        = Channel.of(file(params.sorted_chr_order, checkIfExists:true))
    genomelength             = Channel.of(file(params.genomelength, checkIfExists:true))

    take:
         GATK_exome_bam
         bwa_picard_bam
         capture_ch
         design_ch
    main:
         Genotyping(
          bwa_picard_bam
            .combine(Sites1000g4genotyping)
            .combine(genome)
            .combine(genome_fai)
            .combine(genome_dict)
         )
         CircosPlot_lib(Genotyping.out.loh)
         Read_depth(GATK_exome_bam
            .combine(capture_ch, by:[0])
            .combine(sorted_chr_order)
            .combine(genomelength)
         )

         ch_versions = Genotyping.out.versions.mix(Read_depth.out.versions)


        failed_exon_input = Read_depth.out.read_depth_output
          .map {tuple ->
                def meta = tuple[0]
                def depth = tuple[1]
                def threshold = ''

                if ((meta.sc == 'clin.ex.v1'|| meta.sc == 'idt_v2_plus'|| meta.sc == 'hybrid_selection') && (meta.type == 'tumor_DNA' )) {
                  threshold = params.failed_exon['clin_ex_v1']['Tumor']
                } else if ((meta.sc == 'clin.ex.v1'|| meta.sc == 'idt_v2_plus' || meta.sc == 'hybrid_selection') && (meta.type == 'normal_DNA' || meta.type == 'cell_line_DNA' || meta.type == 'blood_DNA')) {
                  threshold = params.failed_exon['clin_ex_v1']['Normal']
                } else if (meta.sc == 'wholegenome' && meta.type == 'tumor_DNA' ) {
                  threshold = params.failed_exon['wholegenome']['Tumor']
                } else if (meta.sc == 'wholegenome' && (meta.type == 'normal_DNA' || meta.type == 'cell_line_DNA' || meta.type == 'blood_DNA')) {
                  threshold = params.failed_exon['wholegenome']['Normal']
                } else if (meta.sc == 'seqcapez.rms.v1' && meta.type == 'tumor_DNA' ) {
                  threshold = params.failed_exon['seqcapez_rms_v1']['Tumor']
                } else if (meta.sc == 'seqcapez.rms.v1' && (meta.type == 'normal_DNA' || meta.type == 'cell_line_DNA' || meta.type == 'blood_DNA')) {
                  threshold = params.failed_exon['seqcapez_rms_v1']['Normal']
                } else if (meta.sc == 'comp_ex_v1' && meta.type == 'tumor_DNA' ) {
                  threshold = params.failed_exon['comp_ex_v1']['Tumor']
                } else if (meta.sc == 'comp_ex_v1' && (meta.type == 'normal_DNA' || meta.type == 'cell_line_DNA' || meta.type == 'blood_DNA')) {
                  threshold = params.failed_exon['comp_ex_v1']['Normal']
                } else if ((meta.sc == 'clin.snv.v1'|| meta.sc == 'clin.snv.v2') && (meta.type =='tumor_DNA' )) {
                  threshold = params.failed_exon['clin_snv']['Tumor']
                } else if ((meta.sc == 'clin.snv.v1'|| meta.sc == 'clin.snv.v2') && (meta.type =='normal_DNA' || meta.type == 'cell_line_DNA' || meta.type == 'blood_DNA')) {
                  threshold = params.failed_exon['clin_snv']['Normal']
                } else if (meta.sc == 'xgen-hyb-panelv2' && meta.type == 'tumor_DNA' ) {
                  threshold = params.failed_exon['xgen_hyb_panelv2']['Tumor']
                } else if (meta.sc == 'xgen-hyb-panelv2' && (meta.type == 'normal_DNA' || meta.type == 'cell_line_DNA' || meta.type == 'blood_DNA')) {
                  threshold = params.failed_exon['xgen_hyb_panelv2']['Normal']
                } else if ((meta.sc == 'seqcapez.hu.ex.v3'|| meta.sc == 'seqcapez.hu.ex.utr.v1') && (meta.type == 'tumor_DNA' )) {
                  threshold = params.failed_exon['seqcapez_hu_ex_v3']['Tumor']
                } else if ((meta.sc == 'seqcapez.hu.ex.v3'|| meta.sc == 'seqcapez.hu.ex.utr.v1') && (meta.type == 'normal_DNA' || meta.type == 'cell_line_DNA' || meta.type == 'blood_DNA')) {
                  threshold = params.failed_exon['seqcapez_hu_ex_v3']['Normal']
                } else if (meta.sc == 'agilent.v7' && meta.type == 'tumor_DNA' ) {
                  threshold = params.failed_exon['agilent_v7']['Tumor']
                } else if (meta.sc == 'agilent.v7' && (meta.type == 'normal_DNA' || meta.type == 'cell_line_DNA' || meta.type == 'blood_DNA')) {
                  threshold = params.failed_exon['agilent_v7']['Normal']
                }
                return [meta,depth,threshold]
          }
         failed_exon_input|FailedExons_Genes
         //FailedExons_Genes(Read_depth.out.read_depth_output)
         /*
         Bam2tdf(
          GATK_exome_bam
             .combine(genome)
             .combine(genome_fai)
             .combine(genome_dict)
         )
        */
         //ch_versions = ch_versions.mix(Bam2tdf.out.versions)

         Flagstat(GATK_exome_bam)

         ch_versions = ch_versions.mix(Flagstat.out.versions)

         Bamutil(
          GATK_exome_bam
             .combine(genome)
             .combine(genome_fai)
             .combine(genome_dict)
         )

         ch_versions = ch_versions.mix(Bamutil.out.versions)

        HotspotPileup(
          GATK_exome_bam
             .combine(genome)
             .combine(genome_fai)
             .combine(genome_dict)
             .combine(hg19_hotspot_pos)
        )
        VerifyBamID(
          GATK_exome_bam
             .combine(recode_vcf)
        )

        ch_versions = ch_versions.mix(VerifyBamID.out.versions)

        Coverage(
        GATK_exome_bam
            .combine(capture_ch, by:[0])
            .combine(sorted_chr_order)
            .combine(genomelength)
        )

        ch_versions = ch_versions.mix(Coverage.out.versions)

        TargetIntervals(
          GATK_exome_bam.combine(capture_ch, by:[0]).combine(design_ch, by:[0])
        )

        ch_versions = ch_versions.mix(TargetIntervals.out.versions)

        Conpair_pile(
          GATK_exome_bam
          .combine(genome)
          .combine(genome_fai)
          .combine(genome_dict)
          .combine(conpair_refbed)
        )

        HSMetrics(GATK_exome_bam
          .combine(TargetIntervals.out.probe_intervals, by:[0])
          .combine(TargetIntervals.out.target_intervals, by:[0])
          .combine(genome)
          .combine(genome_fai)
          .combine(genome_dict)
          )

        ch_versions = ch_versions.mix(HSMetrics.out.versions)

        Exome_QC(
          GATK_exome_bam.combine(capture_ch, by:[0]).combine(HSMetrics.out.hsmetrics_out, by:[0])
        )
        Hotspot_Coverage(
        GATK_exome_bam
             .combine(chrom_sizes)
             .combine(access_hotspot)
        )

   emit:
         hotspot_pileup = HotspotPileup.out
         coverage = Coverage.out.coverage_out
         loh = Genotyping.out.loh
         gt = Genotyping.out.gt
         Exome_qc = Exome_QC.out
         hsmetrics = HSMetrics.out.hsmetrics_out
         flagstat = Flagstat.out.flagstat
         verifybamid = VerifyBamID.out.verifybamid_out
         ch_versions = ch_versions
         conpair_pileup = Conpair_pile.out
         hotspot_depth = Hotspot_Coverage.out

}
