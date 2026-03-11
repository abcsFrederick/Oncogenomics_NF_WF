include {Common_RNAseq_WF} from './Common_RNAseq_WF'
include {Exome_common_WF} from './Exome_common_WF.nf'
include {MakeHotSpotDB
        Hotspot_Boxplot
        CoveragePlot} from '../modules/qc/plots'
include {Exome_QC} from '../modules/qc/qc.nf'
include {Vcf2txt} from '../modules/misc/snpEff'
include {FormatInput
        AddAnnotation as AddAnnotation_exome
        AddAnnotation as AddAnnotation_rnaseq} from '../modules/annotation/annot'
include {Annotation} from '../subworkflows/Annotation'
include {DBinput_exome_rnaseq} from '../modules/misc/DBinput'
include {CNVkitPooled
        CNVkit_png} from '../modules/cnvkit/CNVkitPooled'
include {RECONCNV} from '../modules/cnvkit/reconcnv/main'
include {QC_summary_Patientlevel} from '../modules/qc/qc'
include {Genotyping_Sample
        Multiqc
        CircosPlot} from '../modules/qc/qc'
include {Actionable_fusion} from '../modules/Actionable.nf'
include {Fusion_Annotation
        Merge_fusion_annotation} from '../modules/annotation/Fusion_Annotation'
include {Combine_customRNAQC
        RNAqc_TrancriptCoverage} from '../modules/qc/picard'
include {CUSTOM_DUMPSOFTWAREVERSIONS} from '../modules/nf-core/dumpsoftwareversions/main.nf'
include {Allstepscomplete} from '../modules/misc/Allstepscomplete'


def combinelibraries(inputData) {
    def processedData = inputData.map { meta, file ->
        meta2 = [
            id: meta.id,
            casename: meta.casename,
            diagnosis: meta.diagnosis
        ]
        [meta2, file]
    }.groupTuple()
     .map { meta, files -> [meta, [*files]] }
}

/*
def metadatareducer(inputChannel) {
    return inputChannel.map { meta, file ->
            [ meta.id, meta.casename, meta.diagnosis, file ]}
        .reduce([[:], []]) { result, item ->
            def (patient, casename, diagnosis, filePath) = item

            // Dynamically set the id from the meta data
            result[0] = [id: patient, casename: casename, diagnosis:diagnosis]

            // Append file paths
            result[1].add(filePath)


            return result
        }
}
*/
def metadatareducer(inputData) {
    return inputData
        .map { meta, file ->
            [ meta.id, meta.casename, meta.diagnosis, file ]
        }
        .groupTuple()
        .map { patient, casenames, diagnoses, files ->
            def metadata = [
                id: patient,
                casename: casenames[0],
                diagnosis: diagnoses[0]
            ]
            [metadata, files]
        }
}


    genome                  = Channel.of(file(params.genome, checkIfExists:true))
    genome_fai              = Channel.of(file(params.genome_fai, checkIfExists:true))
    genome_dict             = Channel.of(file(params.genome_dict, checkIfExists:true))
    pfamdb  = Channel.of(file(params.pfamdb, checkIfExists:true))
    genome_version_fusion_annotation =  Channel.from(params.genome_version_fusion_annotation)
    genome_version = Channel.from(params.genome_version)
    Pipeline_version = Channel.from(params.Pipeline_version)
    recon_config_ch = Channel.fromPath(params.reconcnv_config_file)
    recon_data_ch   = Channel.fromPath(params.reconcnv_data_dir)


workflow Tumor_RNAseq_WF {

take:
    tumor_rnaseq_samplesheet


main:
// Parse the samplesheet to generate fastq tuples
//samples = Channel.fromPath("Tumor_RNAseq.csv")
samples = tumor_rnaseq_samplesheet
.splitCsv(header:true)
.filter { row -> row.type == "tumor_DNA" || row.type == "cell_line_DNA" || row.type == "tumor_RNA" || row.type == "cell_line_RNA" || row.type == "normal_DNA" }
.map { row ->
    def meta = [:]
    meta.id    =  row.sample
    meta.lib   =  row.library
    meta.sc    =  row.sample_captures
    meta.casename  = row.casename
    meta.type     = row.type
    meta.diagnosis =row.Diagnosis
    def fastq_meta = []
    fastq_meta = [ meta,  file(row.read1), file(row.read2)  ]

    return fastq_meta
}
ch_allcomplete = Channel.empty()
samples_branch = samples.branch{
        exome: it[0].type == "tumor_DNA" || it[0].type == "cell_line_DNA" || it[0].type == "normal_DNA"
        rnaseq: it[0].type == "tumor_RNA" || it[0].type == "cell_line_RNA"
}

samples_branch.rnaseq|Common_RNAseq_WF
samples_branch.exome|Exome_common_WF

ch_versions = Exome_common_WF.out.ch_versions.mix(Common_RNAseq_WF.out.ch_versions)

pileup_Channel = Exome_common_WF.out.pileup.concat(Common_RNAseq_WF.out.pileup)
pileup_pair = metadatareducer(pileup_Channel)


MakeHotSpotDB(pileup_pair)

snpeff_vcf2txt_ch = Exome_common_WF.out.HC_snpeff_snv_vcf2txt.concat(Common_RNAseq_WF.out.snpeff_vcf)
Combined_snpeff_vcf2txt_ch = metadatareducer(snpeff_vcf2txt_ch)
format_input_ch = Combined_snpeff_vcf2txt_ch.join(MakeHotSpotDB.out,by:[0])
FormatInput(format_input_ch)
Annotation(FormatInput.out)
ch_versions = ch_versions.mix(Annotation.out.version)
add_annotation_exome_input = Exome_common_WF.out.HC_snpeff_snv_vcf2txt.combine(Annotation.out.rare_annotation.map{meta, file -> [file]})
add_annotation_rnaseq_input = Common_RNAseq_WF.out.snpeff_vcf.combine(Annotation.out.rare_annotation.map{meta, file -> [file]})
add_annotation_exome_input|AddAnnotation_exome
add_annotation_rnaseq_input|AddAnnotation_rnaseq
snpeff_vcf2txt_ch = Exome_common_WF.out.HC_snpeff_snv_vcf2txt.concat(Common_RNAseq_WF.out.snpeff_vcf)
annot_ch_dbinput = AddAnnotation_exome.out.concat(AddAnnotation_rnaseq.out)
/*
Combined_dbinput_ch = annot_ch_dbinput.map { meta, file ->
            [ meta.id, meta.casename, meta, file ]
        }
        .map { patient, casename, meta, file ->
            def meta2 = [
                lib: meta.lib,
                sc: meta.sc,
                type: meta.type
            ]
            [patient, casename, meta2, file]
        }
    .reduce([[:], [], [], [], []]) { result, item ->
        def (patient, casename, meta, filePath) = item

        // Populate result[0] with patient and casename
        result[0] = [id: patient, casename: casename]

        // Populate result[1] with lib and type
        result[1].add(meta.lib)
        result[1].add(meta.type)

        // Populate result[2] with lib and sc
        result[2].add(meta.lib)
        result[2].add(meta.sc)

        if (filePath.toString().contains('DNA')) {result[3] << filePath.toString()}

        if (filePath.toString().contains('RNA')) {result[4] << filePath.toString()}

        return result
    }
*/

Combined_dbinput_ch = annot_ch_dbinput.map { meta, file ->
            [ meta.id, meta.casename, meta.lib, meta.sc, meta.type, file ]
        }
        .groupTuple()
        .map { patient, casenames, libs, scs, types, files ->
            def meta = [id: patient, casename: casenames[0]]
            def libsAndTypes = libs.collect { [it, types[libs.indexOf(it)]] }.flatten()
            def libsAndScs = libs.collect { [it, scs[libs.indexOf(it)]] }.flatten()
            def dnaFiles = []
            def rnaFiles = []

            files.each { filePath ->
                if (filePath.toString().contains('DNA')) {dnaFiles << filePath}
                if (filePath.toString().contains('RNA')) {rnaFiles << filePath}
            }

            [meta, libsAndTypes, libsAndScs, dnaFiles, rnaFiles]
        }
DBinput_exome_rnaseq(Combined_dbinput_ch,
                    Combined_snpeff_vcf2txt_ch.map{meta, file -> file})
ch_allcomplete = ch_allcomplete.mix( DBinput_exome_rnaseq.out.map { all -> all[1..-1] }.flatten())

cnvkit_input_bam = Exome_common_WF.out.exome_final_bam.branch{
        tumor: it[0].type == "tumor_DNA" || it[0].type == "cell_line_DNA"  }
        .map { tuple ->
        def meta = tuple[0]
        def bam = tuple[1]
        def bai = tuple[2]
        def cnv_ref = ''

        if (meta.sc == 'clin.ex.v1') {
            cnv_ref = params.cnvkit_clin_ex_v1
        } else if (meta.sc == 'clin.snv.v1') {
            cnv_ref = params.cnvkit_clin_snv_v1
        } else if (meta.sc == 'clin.cnv.v2') {
            cnv_ref = params.cnvkit_clin_cnv_v2
        } else if (meta.sc == 'clin.snv.v2') {
            cnv_ref = params.cnvkit_clin_snv_v2
        } else if (meta.sc == 'seqcapez.rms.v1') {
            cnv_ref = params.cnvkit_seqcapez_rms_v1
        } else if (meta.sc == 'seqcapez.hu.ex.v3') {
            cnv_ref = params.cnvkit_seqcapez_hu_ex_v3
        } else if (meta.sc == 'agilent.v7') {
            cnv_ref = params.cnvkit_agilent.v7
        } else {
            return [meta, bam, bai]
        }
        return [meta, bam, bai, cnv_ref]
}
.filter { tuple ->
    tuple.size() == 4
}

cnvkit_input_bam|CNVkitPooled

cnvkit_recon_input = CNVkitPooled.out.cnvkit_cnr
    .join(CNVkitPooled.out.cnvkit_cns, by:[0])
    .map { meta, cnr, cns -> [meta, cnr, cns] }

RECONCNV(
    cnvkit_recon_input,
    recon_config_ch,
    recon_data_ch
)

CNVkitPooled.out.cnvkit_pdf|CNVkit_png
ch_allcomplete = ch_allcomplete.mix(
    CNVkit_png.out.map { meta, file -> file }.ifEmpty([]) )

hotspot_ch = Exome_common_WF.out.hotspot_depth.concat(Common_RNAseq_WF.out.hotspot_depth)
combined_hotspot_ch = metadatareducer(hotspot_ch)
Hotspot_Boxplot(combined_hotspot_ch)
ch_allcomplete = ch_allcomplete.mix( Hotspot_Boxplot.out.map { meta, file -> file } )

genotyping_ch = Exome_common_WF.out.gt.concat(Common_RNAseq_WF.out.gt)
combined_genotyping_ch = metadatareducer(genotyping_ch)
Genotyping_Sample(combined_genotyping_ch,
                Pipeline_version)
ch_allcomplete = ch_allcomplete.mix( Genotyping_Sample.out.map { all -> all[1..-1] }.flatten())

loh_ch = Exome_common_WF.out.loh.concat(Common_RNAseq_WF.out.loh)
combined_loh_ch = metadatareducer(loh_ch)
CircosPlot(combined_loh_ch)
ch_allcomplete = ch_allcomplete.mix( CircosPlot.out.map { meta, file -> file } )

coverage_ch = Exome_common_WF.out.coverage.concat(Common_RNAseq_WF.out.coverage)
combined_coverage_ch = metadatareducer(coverage_ch)
CoveragePlot(combined_coverage_ch)
ch_allcomplete = ch_allcomplete.mix( CoveragePlot.out.map { meta, file -> file } )


 //RNA lib processing steps
actionable_fusion_input = Common_RNAseq_WF.out.fusion_calls.map{ meta, fusion -> [meta, [fusion]] }
Actionable_fusion(actionable_fusion_input)
ch_allcomplete = ch_allcomplete.mix( Actionable_fusion.out.map { all -> all[1..-1] }.flatten())

Fusion_Annotation_input = Common_RNAseq_WF.out.rsem_isoforms
                        .join(Common_RNAseq_WF.out.fusion_calls, by:[0])
                        .combine(pfamdb)
                        .combine(genome)
                        .combine(genome_version_fusion_annotation)
                        .combine(genome_version)
Fusion_Annotation(Fusion_Annotation_input)

merge_fusion_anno_input = combinelibraries(Fusion_Annotation.out)
Merge_fusion_annotation(merge_fusion_anno_input.combine(genome_version))
ch_allcomplete = ch_allcomplete.mix( Merge_fusion_annotation.out.map { all -> all[1..-1] }.flatten())



qc_summary_input_ch = combinelibraries(Exome_common_WF.out.exome_qc)
QC_summary_Patientlevel(qc_summary_input_ch)
ch_allcomplete = ch_allcomplete.mix( QC_summary_Patientlevel.out.map { meta, file -> file } )

customRNAqc_ch = combinelibraries(Common_RNAseq_WF.out.rnalib_custum_qc)
Combine_customRNAQC(customRNAqc_ch)
ch_allcomplete = ch_allcomplete.mix( Combine_customRNAQC.out.map { meta, file -> file } )

transcriptcovRNAqc_ch = combinelibraries(Common_RNAseq_WF.out.picard_rnaseqmetrics)
RNAqc_TrancriptCoverage(transcriptcovRNAqc_ch)
ch_allcomplete = ch_allcomplete.mix( RNAqc_TrancriptCoverage.out.map { meta, file -> file } )



rnaseq_qc_ch = Common_RNAseq_WF.out.Fastqc_out.join(Common_RNAseq_WF.out.pileup, by: [0])
                      .join(Common_RNAseq_WF.out.chimeric_junction, by: [0])
                      .join(Common_RNAseq_WF.out.rsem_genes, by: [0])
                      .join(Common_RNAseq_WF.out.rnaseqc, by: [0])
                      .join(Common_RNAseq_WF.out.circos_plot, by: [0])
                      .join(Common_RNAseq_WF.out.strandedness, by: [0])
                      .join(Common_RNAseq_WF.out.rnalib_custum_qc, by: [0])
                      .join(Common_RNAseq_WF.out.picard_rnaseqmetrics, by: [0])
                      .join(Common_RNAseq_WF.out.picard_rnaseqmetrics_pdf, by: [0])
                      .join(Common_RNAseq_WF.out.picard_alignmetrics, by: [0])
                      .join(Common_RNAseq_WF.out.markdup_txt, by: [0])
                      .join(Common_RNAseq_WF.out.flagstat, by: [0])
                      .join(Common_RNAseq_WF.out.fastq_screen, by: [0])
.map { meta, fastqc, pileup, chimeric, rsem, rnaseqc, circos, strand, rna_qc, picardqc, picardqc_pdf, picardqc_metric, markdup, flagstat, fastq_screen ->
    meta2 = [
        id: meta.id,
        casename: meta.casename
    ]
    [ meta2, fastqc, pileup, chimeric, rsem, rnaseqc, circos, strand, rna_qc, picardqc, picardqc_pdf, picardqc_metric, markdup, flagstat, fastq_screen ]
  }.groupTuple()
   .map { meta, fastqcs, pileups, chimerics, rsems, rnaseqcs, circoss, strands, rna_qcs, picardqcs, picardqc_pdfs, picardqc_metrics, markdups, flagstats, fastq_screens -> [ meta, *fastqcs, *pileups, *chimerics, *rsems, *rnaseqcs, *circoss, *strands, *rna_qcs, *picardqcs, *picardqc_pdfs, *picardqc_metrics, *markdups, *flagstats, *fastq_screens ] }

exome_qc_ch = Exome_common_WF.out.Fastqc_out
            .join(Exome_common_WF.out.verifybamid)
            .join(Exome_common_WF.out.flagstat)
            .join(Exome_common_WF.out.exome_final_bam)
            .join(Exome_common_WF.out.hsmetrics)
            .join(Exome_common_WF.out.krona)
            .join(Exome_common_WF.out.kraken)
            .join(Exome_common_WF.out.exome_qc)
            .join(Exome_common_WF.out.markdup_txt)
            .join(Exome_common_WF.out.fastq_screen)
.map { meta, fastqc, verifybamid, flagstat, bam, bai, hsmetrics, krona, kraken, exome_qc, markdup, fastq_screen ->
    meta2 = [
        id: meta.id,
        casename: meta.casename
    ]
    [ meta2, fastqc, verifybamid, flagstat, bam, bai, hsmetrics, krona, kraken, exome_qc, markdup, fastq_screen ]
  }.groupTuple()
   .map { meta, fastqcs, verifybamids, flagstats, bams, bais, hsmetricss, kronas, krakens, exome_qcs, markdups, fastq_screens -> [ meta, *fastqcs, *verifybamids, *flagstats, *bams, *bais, *hsmetricss, *kronas, *krakens, *exome_qcs, *markdups, *fastq_screens ] }


Combined_multiqc_input = exome_qc_ch.merge(rnaseq_qc_ch) { item1, item2 ->
    if (item1[0].id == item2[0].id && item1[0].casename == item2[0].casename) {
        return [[id: item1[0].id, casename: item1[0].casename]] + [item1[1..-1] + item2[1..-1]]
    } else {
        return null
    }

}

Multiqc(Combined_multiqc_input)

ch_versions = ch_versions.mix(Multiqc.out.versions)

combine_versions  = ch_versions.unique().collectFile(name: 'collated_versions.yml')

custom_versions_input = Multiqc.out.multiqc_report
        .combine(combine_versions).map{ meta, multiqc, version -> [meta, version] }
        .combine(Pipeline_version)

CUSTOM_DUMPSOFTWAREVERSIONS(custom_versions_input)

Allstepscomplete(CUSTOM_DUMPSOFTWAREVERSIONS.out.config,
                ch_allcomplete)


}
