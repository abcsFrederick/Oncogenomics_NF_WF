#!/usr/bin/env nextflow
import Utils

nextflow.enable.dsl = 2

// Default parameter values

//OUTDIR = "${params.outdir}"

def patientid_val = Channel.value(null)
def casename_val = Channel.value(null)


log.info """\
         E X O M E - R N A S E Q - N F   P I P E L I N E
         ===================================
         NF version   : $nextflow.version
         runName      : $workflow.runName
         username     : $workflow.userName
         configs      : $workflow.configFiles
         cmd line     : $workflow.commandLine
         start time   : $workflow.start
         projectDir   : $workflow.projectDir
         launchDir    : $workflow.launchDir
         workdDir     : $workflow.workDir
         homeDir      : $workflow.homeDir
         samplesheet  : ${params.samplesheet}
         outdir       : ${params.resultsdir}
         genome_version       : ${params.genome_v}
         platform      : ${params.platform}
         """
         .stripIndent()

Channel
  .fromPath(params.samplesheet)
  .splitCsv(header: true)
  .map { row ->
    [
      sample   : (row.sample ?: '').toString(),
      casename : (row.casename ?: '').toString(),
      library      : (row.library ?: '').toString(),
      read1    : (row.read1 ?: '').toString(),
      read2    : (row.read2 ?: '').toString(),
      bam      : (row.bam   ?: '').toString(),
      sample_captures : (row.sample_captures ?: '').toString(),
      Matched_RNA    : (row.Matched_RNA ?: '').toString(),
      Matched_normal : (row.Matched_normal ?: '').toString(),
      Diagnosis      : (row.Diagnosis ?: '').toString(),
      type           : (row.type ?: '').toString(),
      filetype : ((row.file_type ?: row.filetype ?: '') as String).toLowerCase()
    ]
  }
  .set { ROWS }

// Partition rows
FASTQ_ROWS = ROWS.filter { it.read1 }     // rows providing FASTQs
BAM_ROWS   = ROWS.filter { it.bam }       // rows providing BAM/CRAM

FASTQ_ROWS
  .map { meta ->
    tuple(meta, file(meta.read1), meta.read2 ? file(meta.read2) : null)
  }
  .set { FASTQ_IN }

FASTQ_IN
  .map { meta, r1, r2 -> tuple("${meta.sample}||${meta.casename}||${meta.library}", meta, r1, r2) }
  .groupTuple(by: 0)
  .map { key, metas, r1s, r2s ->
    def meta = metas[0]  // representative; carries all your extra columns
    def sortedR1s = r1s.sort { it.toString() }
    def sortedR2s = r2s.findAll { it != null }.sort { it.toString() }
    tuple(meta, sortedR1s, sortedR2s)
  }
  .set { FASTQ_GROUPED }

// channel of samples with >1 R1 file OR >1 R2 file
TO_MERGE = FASTQ_GROUPED.filter { meta, r1s, r2s ->
  ((r1s?.size() ?: 0) > 1) || ((r2s?.size() ?: 0) > 1)
}

// channel of samples with exactly 1 R1 AND exactly 1 R2
PASSTHRU = FASTQ_GROUPED.filter { meta, r1s, r2s ->
  ((r1s?.size() ?: 0) == 1) && ((r2s?.size() ?: 0) == 1)
}
.map { meta, r1s, r2s ->
  def r1 = r1s[0]
  def r2 = r2s[0]
  tuple(meta, r1, r2)
}

BAM_ROWS
  .map { meta -> tuple(meta, file(meta.bam)) }
  .set { BAM_IN }

process MERGE_FASTQS {
  tag { "${meta.sample}_${meta.library}" }

  //publishDir "${params.resultsdir}/fastq", mode:'copy', overwrite:true

  input:
    tuple val(meta), path(r1_lanes), path(r2_lanes)

  output:
    tuple val(meta),
          path("${meta.sample}_${meta.library}_R1.fastq.gz"),
          path("${meta.sample}_${meta.library}_R2.fastq.gz")

  script:
  """
  # R1: merge fastq
cat ${r1_lanes.join(' ')} > ${meta.sample}_${meta.library}_R1.fastq.gz

# R2: same logic;
cat ${r2_lanes.join(' ')} > ${meta.sample}_${meta.library}_R2.fastq.gz
"""
}


process BAM_TO_FASTQ {
  tag { "${meta.sample}_${meta.library}" }

  //publishDir "${params.resultsdir}/fastq1", mode:'copy', overwrite:true

  input:
    tuple val(meta), path(bam)

  output:
    tuple val(meta),
          path("${meta.sample}_${meta.library}_R1.fastq.gz"),
          path("${meta.sample}_${meta.library}_R2.fastq.gz")

  script:
  """
  set -euo pipefail
#  module -q load samtools || true
#  module load picard

  # Sort the BAM by queryname (required for SamToFastq)
  samtools sort -n -@ ${task.cpus} -m 16G -o ${meta.sample}_${meta.library}.queryname.bam ${bam}

  # Use Picard SamToFastq to split into R1/R2
  java -Xmx40g -jar \$PICARDJAR SamToFastq \
      I=${meta.sample}_${meta.library}.queryname.bam \
      FASTQ=${meta.sample}_${meta.library}_R1.fastq.gz \
      SECOND_END_FASTQ=${meta.sample}_${meta.library}_R2.fastq.gz
  """
}

process PREPARE_SAMPLESHEET {

    input:
    path samplesheet
    val genome_version

    output:
    path("*csv"), emit: csv_files
    env patientid, emit: patientid
    env casename, emit: casename

    script:
    """

    if [ "${genome_version}" == "hg19" ] || [ "${genome_version}" == "hg38" ]; then
        python ${workflow.projectDir}/bin/prepare_samplesheet.py ${samplesheet} .
    elif [ ${genome_version} == "mm39" ]; then
        cp ${samplesheet} mouse_rnaseq.csv
    else
        echo "Error: Unknown genome: ${genome_version}"
        exit 1
    fi
    patientid=\$(awk -F',' 'NR==1 {for (i=1; i<=NF; i++) if (\$i=="sample") s=i} NR>1 {print \$s}' "${samplesheet}" | sort | uniq )
    casename=\$(awk -F',' 'NR==1 {for (i=1; i<=NF; i++) if (\$i=="casename") s=i} NR>1 {print \$s}' "${samplesheet}" | sort | uniq )

    export patientid
    export casename
    """
}

include {RNAseq_only} from './workflows/RNAseq_only.nf'
include {RNAseq_multiple_libs} from './workflows/RNAseq_multiple_libs.nf'
include {Exome_only_WF} from './workflows/Exome_only_WF.nf'
include {Tumor_multiple_libs} from './workflows/Tumor_multiple_libs.nf'
include {Tumor_Normal_WF} from './workflows/Tumor_Normal_WF.nf'
include {Tumor_Normal_RNAseq_WF} from './workflows/Tumor_Normal_RNAseq_WF.nf'
include {Tumor_RNAseq_WF} from './workflows/Tumor_RNAseq_WF.nf'
include {Mouse_RNA} from './workflows/Mouse_RNA.nf'


workflow {

TO_MERGE | MERGE_FASTQS
//FASTQ_MERGED = MERGE_FASTQS.out

BAM_IN|BAM_TO_FASTQ
//FASTQ_FROM_BAM = BAM_TO_FASTQ.out

def HEADER = [
  'sample','casename','library','read1','read2',
  'bam','sample_captures','Matched_RNA','Matched_normal','Diagnosis','type','filetype'
]
// Combine streams and render rows using meta + the produced R1/R2
def data_rows = MERGE_FASTQS.out
  .mix(PASSTHRU)
  .mix(BAM_TO_FASTQ.out)
  .map { meta, r1, r2 ->
    // build row values in header order
    def rowmap = meta + [ read1: r1.toString(), read2: r2.toString() ]
    HEADER.collect { k -> (rowmap[k] ?: '').toString() }.join(',')
  }

data_rows
  .collect()
  .map { rows -> ([HEADER.join(',')] + rows).join('\n') + '\n' }
  .collectFile(name: 'merged_samplesheet.csv', storeDir: workflow.workDir)
  .set { MERGED_SHEET }

  // 4) Splitter on the merged sheet
def prepared_samplesheets = PREPARE_SAMPLESHEET( MERGED_SHEET, params.genome_v )

//def prepared_samplesheets = PREPARE_SAMPLESHEET(params.samplesheet,params.genome_v)
    //prepared_samplesheets.patientid.view()


    // SAFELY extract values
    prepared_samplesheets.patientid.subscribe { patientid_val = it }
    prepared_samplesheets.casename.subscribe { casename_val = it }


    prepared_samplesheets.csv_files.branch {
        rnaseq: it.name == 'RNAseq.csv'
        exome: it.name == 'Exome.csv'
        multiple_exome: it.name == 'Tumor_lib.csv'
        tumor_rnaseq_normal: it.name == 'Tumor_RNAseq_Normal.csv'
        multiple_rna: it.name == 'RNA_lib.csv'
        tumor_normal: it.name == 'Tumor_Normal.csv'
        tumor_rnaseq: it.name == 'Tumor_RNAseq.csv'
        mouse_rna: it.name == 'mouse_rnaseq.csv'


    }.set { branched_samplesheets }
    branched_samplesheets.rnaseq | RNAseq_only
    branched_samplesheets.tumor_rnaseq | Tumor_RNAseq_WF
    branched_samplesheets.exome | Exome_only_WF
    branched_samplesheets.multiple_exome | Tumor_multiple_libs
    branched_samplesheets.tumor_rnaseq_normal | Tumor_Normal_RNAseq_WF
    branched_samplesheets.multiple_rna | RNAseq_multiple_libs
    branched_samplesheets.tumor_normal | Tumor_Normal_WF
    branched_samplesheets.mouse_rna | Mouse_RNA
}


workflow.onComplete {
    println "✅ Done!"
    println "🧬 PatientID: ${patientid_val}"
    println "📁 Casename: ${casename_val}"
    println "📁 Workflow: ${params.platform}"


    def message = Utils.handleWorkflowCompletion(
        workflow,
        params.genome_v,
        params.genome_v == "mm39" ? "completed.txt" : "successful.txt",
        params.platform,
        patientid_val,
        casename_val,
        params.resultsdir
    )
    if ( params.platform == "biowulf" ) {
        sendMail(
            to: "${workflow.userName}@mail.nih.gov",
            cc: workflow.profile == "biowulf_mouse_RNA_slurm" ? "" : "gangalapudiv2@mail.nih.gov",
            subject: workflow.success ? "khanlab ngs-pipeline execution successful" : "khanlab ngs-pipeline execution failed",
            body: message,
            mimeType: 'text/html'
        )
    }
}
