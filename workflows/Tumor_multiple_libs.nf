include {Exome_common_WF} from './Exome_common_WF.nf'
include {MakeHotSpotDB
        Hotspot_Boxplot
        CoveragePlot} from  '../modules/qc/plots'
include {FormatInput} from '../modules/annotation/annot'
include {Annotation} from '../subworkflows/Annotation'
include {AddAnnotation} from '../modules/annotation/annot'
include {CNVkitPooled} from '../modules/cnvkit/CNVkitPooled'
include {RECONCNV} from '../modules/cnvkit/reconcnv/main'
include {CNVkit_png} from '../modules/cnvkit/CNVkitPooled'
include {DBinput_multiple_new} from '../modules/misc/DBinput'
include {Genotyping_Sample
        Multiqc
        CircosPlot
        QC_summary_Patientlevel} from '../modules/qc/qc'
include {TcellExtrect} from '../modules/misc/TcellExtrect'
include {CUSTOM_DUMPSOFTWAREVERSIONS} from '../modules/nf-core/dumpsoftwareversions/main.nf'
include {Allstepscomplete} from '../modules/misc/Allstepscomplete'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

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


workflow Tumor_multiple_libs {


    genome_version_tcellextrect         = Channel.of(params.genome_version_tcellextrect)
    Pipeline_version = Channel.from(params.Pipeline_version)
    recon_config_ch = Channel.fromPath(params.reconcnv_config_file)
    recon_data_ch   = Channel.fromPath(params.reconcnv_data_dir)


take:
    multiple_tumor_samplesheet

main:

samples_exome = multiple_tumor_samplesheet
.splitCsv(header:true)
.filter { row -> row.type == "tumor_DNA" || row.type == "cell_line_DNA" }
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
Exome_common_WF(samples_exome)

makehotspotdb_input = combinelibraries(Exome_common_WF.out.pileup)
MakeHotSpotDB(makehotspotdb_input)

combined_snpefflibs = combinelibraries(Exome_common_WF.out.HC_snpeff_snv_vcf2txt)

format_input = combined_snpefflibs.join(MakeHotSpotDB.out,by:[0])
FormatInput(format_input)
Annotation(FormatInput.out)
ch_versions = Exome_common_WF.out.ch_versions.mix(Annotation.out.version)
add_annotation_input = Exome_common_WF.out.HC_snpeff_snv_vcf2txt.combine(Annotation.out.rare_annotation.map{meta, file -> [file]})
AddAnnotation(add_annotation_input)
annot_ch_dbinput = AddAnnotation.out.map{ meta, file -> [ meta.id, meta.casename, meta, file ] }
            .map { patient, casename, meta, file ->
            meta2 = [
                lib: meta.lib,
                sc: meta.sc,
                type: meta.type
            ]
            [patient,casename, meta2, file]
    }
    .reduce([[:], [], [], []]) { result, item ->
        def (patient, casename, meta, filePath) = item

        // Dynamically set the id from the meta data
        result[0] = [id: patient, casename: casename]

        // Append library and type to the result
        result[1].add(meta.lib)
        result[1].add(meta.type)

        // Append library and source class to the result
        result[2].add(meta.lib)
        result[2].add(meta.sc)

        // Append file paths
        result[3].add(filePath)

        return result
    }

//combined_annotationlibs = combinelibraries(AddAnnotation.out)
//dbinput_input = combined_annotationlibs.join(combined_snpefflibs,by:[0])
DBinput_multiple_new(annot_ch_dbinput,
                    combined_snpefflibs.map{meta, file -> file})
ch_allcomplete = ch_allcomplete.mix( DBinput_multiple_new.out.map { all -> all[1..-1] }.flatten())


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

CNVkitPooled.out.cnvkit_pdf|CNVkit_png
ch_allcomplete = ch_allcomplete.mix(
    CNVkit_png.out.map { meta, file -> file }.ifEmpty([]) )
cnvkit_recon_input = CNVkitPooled.out.cnvkit_cnr
    .join(CNVkitPooled.out.cnvkit_cns, by:[0])
    .map { meta, cnr, cns -> [meta, cnr, cns] }

RECONCNV(
    cnvkit_recon_input,
    recon_config_ch,
    recon_data_ch
)


TcellExtrect(
    Exome_common_WF.out.exome_final_bam
    .join(Exome_common_WF.out.target_capture_ch,by:[0])
    .combine(genome_version_tcellextrect)
)
ch_allcomplete = ch_allcomplete.mix( TcellExtrect.out.naive_txt.map { all -> all[1..-1] }.flatten())

genotyping_input = combinelibraries(Exome_common_WF.out.gt)

Genotyping_Sample(genotyping_input,
                Pipeline_version)
ch_allcomplete = ch_allcomplete.mix( Genotyping_Sample.out.map { all -> all[1..-1] }.flatten())

circos_input = combinelibraries(Exome_common_WF.out.loh)
CircosPlot(circos_input)
ch_allcomplete = ch_allcomplete.mix( CircosPlot.out.map { meta, file -> file } )


hotspot_depth_input = combinelibraries(Exome_common_WF.out.hotspot_depth)
Hotspot_Boxplot(hotspot_depth_input)
ch_allcomplete = ch_allcomplete.mix( Hotspot_Boxplot.out.map { meta, file -> file } )


coverage_plot_input = combinelibraries(Exome_common_WF.out.coverage)
CoveragePlot(coverage_plot_input)
ch_allcomplete = ch_allcomplete.mix( CoveragePlot.out.map { meta, file -> file } )


qc_summary_Patientlevel_input = combinelibraries(Exome_common_WF.out.exome_qc)
QC_summary_Patientlevel(qc_summary_Patientlevel_input)
ch_allcomplete = ch_allcomplete.mix( QC_summary_Patientlevel.out.map { meta, file -> file } )



multiqc_input = Exome_common_WF.out.Fastqc_out.join(Exome_common_WF.out.pileup, by: [0])
                .join(Exome_common_WF.out.kraken, by: [0])
                .join(Exome_common_WF.out.verifybamid, by: [0])
                .join(Exome_common_WF.out.hsmetrics, by: [0])
                .join(Exome_common_WF.out.fastq_screen, by: [0])
                .join(Exome_common_WF.out.flagstat,by: [0])
                .join(Exome_common_WF.out.markdup_txt,by: [0])
                .join(Exome_common_WF.out.krona,by: [0])



multiqc_input.map { meta, fastqc, pileup, kraken, verifybamid, hsmetrics, fastq_screen, flagstat, markdup, krona ->
    meta2 = [
        id: meta.id,
        casename: meta.casename,
    ]
    [ meta2, fastqc, pileup, kraken, verifybamid, hsmetrics, fastq_screen, flagstat, markdup, krona ]
  }.groupTuple()
   .map { meta, fastqcs, pileups, krakens, verifybamids, hsmetricss, fastq_screens, flagstats, markdups, kronas -> [ meta, [*fastqcs,*pileups, *krakens, *verifybamids, *hsmetricss, *fastq_screens, *flagstats, *markdups, *kronas ]] }
   .set { multiqc_ch }

Multiqc(multiqc_ch)
ch_versions = ch_versions.mix(Multiqc.out.versions)

combine_versions  = ch_versions.unique().collectFile(name: 'collated_versions.yml')

custom_versions_input = Multiqc.out.multiqc_report
        .combine(combine_versions).map{ meta, multiqc, version -> [meta, version] }
        .combine(Pipeline_version)

CUSTOM_DUMPSOFTWAREVERSIONS(custom_versions_input)
ch_allcomplete.view()
Allstepscomplete(CUSTOM_DUMPSOFTWAREVERSIONS.out.config,
                ch_allcomplete)


}
