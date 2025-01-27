#!/usr/bin/env nextflow

/*
========================================================================================
    miRPlan Nextflow Workflow
========================================================================================
    Github   :
    Contact  :
----------------------------------------------------------------------------------------
*/


include { samplesheetToList } from 'plugin/nf-schema'

nextflow.enable.dsl=2

// println """\
//          M I R P L A N - N F   P I P E L I N E
//          ===================================
//          genome       : ${params.genome}
//          reads        : ${params.reads}
//          outdir       : ${params.outdir}
//          """
//          .stripIndent()

/*
========================================================================================
    Include Modules
========================================================================================
*/

include { DOWNLOADLIB } from '../modules/local/downloadlib'
include { FASTP       } from '../modules/nf-core/fastp'


/*
========================================================================================
    Include Sub-Workflows
========================================================================================
*/

include { ID_RESOLUTION                           } from "../subworkflows/local/idresolution"
include { QUALITY_CONTROL as QUALITY_CONTROL_RAW  } from "../subworkflows/local/qualitycontrol"
include { QUALITY_CONTROL as QUALITY_CONTROL_TRIM } from "../subworkflows/local/qualitycontrol"
include { VALIDATION                              } from "../subworkflows/local/validation"
include { FILTERING as FILTERING_DB               } from "../subworkflows/local/filtering"
include { FILTERING as FILTERING_GENOME           } from "../subworkflows/local/filtering"
include { QUANTIFICATION                          } from "../subworkflows/local/quantification"
include { DIFFEXPANALYSIS                         } from "../subworkflows/local/diffexpanalysis"
include { ANNOTATION                              } from "../subworkflows/local/annotation"
include { RESULTS_GENERATION                      } from "../subworkflows/local/results"


// Where do we store this function?
// FUNCTIONS

// Define una función para dividir los nombres de archivo
def splitFileName(filename) {
    def parts = filename.split('_') // def se utiliza def para definir una variable local del entorno en el que este
    def name = parts[0]
    def accession = parts[1].split('\\.')[0] // Usa split con regex para evitar el uso de \\.
    return [name, accession]
}

// Esta función comprueba si el fichero es FASTQ 

def validateFastq(file_path) {
    
    def input_file =  Channel.fromPath(file_path)
    
    input_file
        .splitFastq(decompress: true)
        .view()

}

// Esta funcion sirve para validad ficheros de Accession list.
def validateAccessionList(file_path) {

    // List of SRA patterns
    def sra_patterns = ['^SRR', '^ERR', '^DRR']
    def is_accession_list = true

    // Input file
    def input_file =  file(file_path)

    try {
        // Read lines
        input_file.eachLine { line ->

            // Check if the Runs match the SRA patterns
            def matches = sra_patterns.any { pattern ->
                line.trim().matches(pattern + '.*')
            }
            if (!matches) {
                is_accession_list = false
                return
            }
        }
        
        return is_accession_list

    } catch (Exception e) {
        return false
    }

}


/*
========================================================================================
    Create Channels
========================================================================================
*/

// // Canal de entrada para los ficheros con las RUN de cada proyecto
// accession_list_files_path_ch = Channel.fromPath(params.additional_info_dir + '/02-Accession_lists/*.txt')

// Canal con los adaptadores de la secuenciación para el trimming
ch_trimming_adapters_file = Channel.fromPath(params.trimming_adapters).collect() // .collect() devuelve un value channel

// Canales longitud minima y maxima de secuencia (para que el proceso los utilice
// para todos los archivos y no lo ejecute solamente con el primero). Debe ser un
// value channel para que funcione bien.
ch_min_trimming = Channel.value(params.trimming_min_len)
ch_max_trimming = Channel.value(params.trimming_max_len)

// // Canal del directorio de metadatos
// metadata_dir_ch = Channel.fromPath(params.metadata_dir).collect() // value channel

// Canales para los threshold de profundidad de secuenciacion y replicas. Recuerda
// que deben ser canales para que se ejecute para todas las librerias
depth_threshold_ch = Channel.value(params.validation_depth)
rep_threshold_ch = Channel.value(params.validation_rep)

// // Canal del fichero de RNAcentral
// rnacentral_ch = Channel.fromPath(params.f_rnacentral).collect()

// // Numero de mismatches utilizados para el filtrado
// mismatches_filter_ch = Channel.value(params.f_mismatches)

// // Prepare some input channles for QUANTIFICATION
// filt_min_ch_counts = Channel.value(params.q_min_counts_filt)
// filt_min_samples_ch = Channel.value(params.q_min_samples_filt)
// // CHECK these parameters.
// rpm_ch = Channel.value(params.rpm) // CHECK
// project_table_bool_ch = Channel.value(params.q_project)
// project_avg_bool_ch = Channel.value(params.q_project_avg)
// subproject_avg_bool_ch = Channel.value(params.q_subproject_avg)

// // Alpha for differential expression analysis
// alpha_ch = Channel.value(params.e_alpha)

// // Number the mismatches used in the annotation
// mismatches_annot_ch = Channel.value(params.a_mismatches)

// // Rutas de las bases de datos utilizadas para la anotacion de los miRNAs
// // mirbase_ch = Channel.fromPath(params.a_mirbase).collect()
// // pmiren_ch = Channel.fromPath(params.a_pmiren).collect()
// // srnaanno_ch = Channel.fromPath(params.a_srnaanno).collect()

// // Species id info
// species_id_info_ch = Channel.fromPath(params.species_info).collect()

// // Umbral de pvalor utilizado para el test de mww test
// mww_pvalue_thrshld = Channel.value(params.a_mww_pvalue_threshold)

/*
========================================================================================
    WORKFLOW - miRNA GLOBAL ANALYSIS
========================================================================================
*/

workflow MIRPLAN {

    // Fastq files empty channels
    ch_fastq         = Channel.empty()
    pipeline_summary = Channel.empty()
    ch_counts        = channel.empty()

    // Read the input samplesheet and validate it
    Channel
        .fromList(samplesheetToList(params.input, "assets/schema_input.json"))
        .set{ch_input}
    
    /*
    ============================================================================
        Assign a default genome to the input files if necessary.
    ============================================================================
    */

    ch_input
        .map { item ->
            // Get the genome, species name and group
            genome_path = item[4]
            species_name = item[0]

            // Check if a genome has not been provided
            if (genome_path == []) {

                // Get the predetermined genome of this species from the config
                // file. If this species is not in the file, it returns null. 
                predetermined_genome = params.genomes.get(species_name, null)?.fasta
            
            // If a genome has been provided, use that one.
            } else {
                predetermined_genome = genome_path
            }

            // If there is any input without an associated genome, stop the execution.
            if (predetermined_genome == null) {

                // Throw an exception with the list of invalid projects
                throw new RuntimeException("There is no genome associated with the following files:\n${tuple[1]}")
            }

            // Create the new tuple
            [item[0], item[1], item[3], predetermined_genome, item[2], item[5]]

        }.set {ch_input_genome}


    /*
    ============================================================================
        Check the inputs related to counts matrices
    ============================================================================
    */

    // Check the correct usage of the --from_counts parameter.
    ch_input_genome
        .map{item -> item[4]}
        .filter{ it =~ /.*\.tsv$/ }
        .toList()
        .map {files ->
            if (files.isEmpty() && params.from_counts) {
                // Throw an exception with the list of invalid files
                error("\nERROR: The --from_counts parameter must " +
                    "be used only when counts matrix files in TSV " +
                    "format are provided in the samplesheet.\n")
            } else if (!files.isEmpty() && !params.from_counts) {
                error("\nERROR: Counts matrices in TSV format " +
                    "have been provided in the samplesheet, but the " +
                    "--from_counts parameter has not been specified. Please " +
                    "include this parameter to ensure correct processing.\n")
            }
        }

    // Check the correct usage of Group input
    ch_input_genome
        .map { it[5] }
        .toList()
        .map{ item ->

            // Collect valid groups (non-empty elements at index 5)
            def groupsProvided = item.findAll { it != [] }
            
            // Check if all elements have valid groups
            if(groupsProvided.size() != item.size()) {
                if (params.from_counts) {
                    error("\nERROR: Not all input files have an  " +
                        "associated group, even though the --from_counts " +
                        "parameter has been specified. Make sure to use the " +
                        "--from_counts parameter and provide a group only when " +
                        "the input file is a counts matrix.\n")
                }
            }

            if (groupsProvided.size() > 0 && !params.from_counts) {
                // Print a warning if '--from_counts' is false but groups were provided
                log.warn "Groups were detected in the input, but will " +
                    "be ignored because '--from_counts' was not specified. Make " +
                    "sure to use the --from_counts parameter and provide a group " + 
                    "only when the input file is a counts matrix.\n"
            }
        }
    
    /*
    ============================================================================
        SUBWORKFLOW: Download annotation databases
                     and resolve identifier discrepancies
    ============================================================================
    */

    // Get the input species names
    ch_input_genome
        .map{it[0]}
        .unique()
        .collect()
        .set{ch_input_sp_names}

    // Prepare identifiers for the species
    ID_RESOLUTION(ch_input_sp_names)
    
    /*
    ============================================================================
        Separate the fastq files from the accession list files
    ============================================================================
    */

    // Add the assigned identifiers for the species in the meta section of the input channel.
    ch_input_genome.combine(ID_RESOLUTION.out.species_ids, by: 0)
    .map { item ->
            [[species: item[0], species_id: item[6], project: item[1], metadata: item[2], genome: item[3], single_end: true, group_id:item[5]], item[4]]
    }
    .branch { meta, file ->
            
            // Accession lists. It ends with '.txt'
            acclist: file.toString().endsWith('.txt')
                def meta_with_id = [id: meta.project] + meta
                return [meta_with_id, file]

            // Sequencing libraries. They end with '.fastq', '.fastq.gz', '.fq', o '.fq.gz'
            fastq: file.toString() =~ /\.(fastq(\.gz)?|fq(\.gz)?)$/
                def file_wo_extension = file.getName().replaceFirst(/\.(fastq(\.gz)?|fq(\.gz)?)$/, '')
                def meta_with_id = [id: file_wo_extension] + meta
                return [meta_with_id, file]

            // Counts matrices. They end with '.tsv'
            counts: file.toString().endsWith('.tsv')
                def meta_with_id = [id: meta.project] + meta
                return[meta_with_id, file]   
    }
    .set { ch_input_ids }
    

    // Check that the count matrices provided as input are valid.
    if (params.from_counts){

        // Check the counts matrices input files
        VALIDATION(ch_input_ids.counts, 'counts', params.validation_rep, 0)

        // Assign the output channel to the ch_counts channel.
        VALIDATION.out.files
            // Select the valid files
            .filter { tuple ->
                tuple[1] =~ /.*\.valid\.tsv$/
            }.set{ ch_counts }

        // Create the pipeline_summary channel using the meta. OJO. COMO LAS TABLAS DE CONTEOS NO TIENEN LOS PASOS PREVIOS TIENEN NA. ESO HACE QUE LUEGO AL FINAL SE PONGA TODO NA. COMPROBAR
        VALIDATION.out.files
            // Get the sample names from the matrix header 
            .flatMap { meta, file ->
                // Define the new elements of the channel
                def newMeta = meta.findAll { it.key != 'id' }
                def fileHeader = file.withReader { it.readLine() }
                def samples = fileHeader?.split("\t").drop(1) 

                // Crear un nuevo elemento por cada SRR, asignándolo a 'sample'
                samples.collect { srr -> [sample: srr] + newMeta }
            }.set{ pipeline_summary }
    }

    /*
    ============================================================================
        PRE-PROCESSING
    ============================================================================
    */

    // If the input files are count matrices, do not execute the pre-processing.
    if (!params.from_counts){

        /*
        ============================================================================
            SUBWORKFLOW: Download the study libraries
        ============================================================================
        */

        // Validate accession list files
        ch_input_ids.acclist
            .map { tuple ->
                def validity = validateAccessionList(tuple[1])
                [meta: tuple[0], file: tuple[1], validity: validity ]
            }
            .toList() // Collect results as a list
            .map { results ->
                // Check if there are any invalid files
                def invalid_files = results.findAll { it.validity == false }
                if (invalid_files) {
                    // Create a list of invalid file paths
                    def invalid_files_list = invalid_files.collect { it.file.toString() }.join('\n')
                    // Throw an exception with the list of invalid files
                    throw new RuntimeException("Invalid accession list files found:\n${invalid_files_list}")
                }
            }

        // Download the libraries
        DOWNLOADLIB(ch_input_ids.acclist)

        // Split the tuple created in download_libraries
        // [arth, PRJNA157121, [/path/file.fastq.gz,/path/file2.fastq.gz]] ->
        // [arth, PRJNA157121, /path/file.fastq.gz] ... [arth, PRJNA157121, /path/file2.fastq.gz]
        DOWNLOADLIB.out.fastqgz
            .transpose()
            .set{ ch_fastq }

        // Change the meta.id from project to file.
        ch_fastq
            .map { meta, file ->
                def updatedMeta = meta.clone()
                updatedMeta.id = file.getName().replaceFirst(/\.(fastq(\.gz)?|fq(\.gz)?)$/, '')
                return [updatedMeta, file]
            }
            .set {ch_fastq}
        
        // Add the Input information to the summary channel
        ch_input_ids.fastq
            .map { meta, file ->
                [meta.id, [sample: meta.id, species: meta.species, species_id: meta.species_id, project: meta.project, input: 'Local']]
            }
            .set{pipeline_summary}

        ch_fastq
            .map { meta, file ->
                [meta.id, [sample: meta.id, species: meta.species, species_id: meta.species_id, project: meta.project, input: 'Downloaded']]
            }
            .concat(pipeline_summary)
            .set{pipeline_summary}
        
        // Combine the downloaded libraries with those provided by the user in the same channel.
        ch_fastq
            .concat(ch_input_ids.fastq)
            .set{ch_fastq}
        
        /*
        ============================================================================
            SUBWORKFLOW: Perform quality control of RAW data
        ============================================================================
        */

        if (!(params.skip_qc || params.skip_fastqc)) {
            QUALITY_CONTROL_RAW(
                ch_fastq,
                params.skip_multiqc,
                'Raw'
            )
        }

        /*
        ============================================================================
            SUBWORKFLOW: Perform trimming of the libraries using fastp.
        ============================================================================
        */

        // Execute FASTP
        FASTP(
            ch_fastq,
            ch_trimming_adapters_file,
            false,
            false,
            false
        )

        // Add the trimming data to the pipeline_summary channel
        FASTP.out.reads
            .map{meta, file ->
                [meta.id, meta, file]
            }
            .join(pipeline_summary, remainder: true)
            .map { item ->
                // Define new_meta variable
                def new_meta = item.last()
                // Update meta var
                if (item[1] == null){
                    new_meta = new_meta + [Trimming: 'Discarded']
                } else {
                    new_meta = new_meta + [Trimming: 'Trimmed']
                }
                return [item[0], new_meta]
            }
            .set{pipeline_summary}
        
        // Change the meta.id from file to project.
        FASTP.out.reads
            .map { meta, file ->
                def updatedMeta = meta.clone()
                updatedMeta.id = meta.project
                return [updatedMeta, file]
            }
            .groupTuple(by: 0, sort:true)
            .set{ch_fastq}

        /*
        ============================================================================
            SUBWORKFLOW: Perform quality control of RAW data
        ============================================================================
        */

        if (!(params.skip_qc || params.skip_fastqc)) {
            QUALITY_CONTROL_TRIM(
                FASTP.out.reads,
                params.skip_multiqc,
                'Trimmed'
            )
        }

        /*
        ============================================================================
            SUBWORKFLOW: Check if the projects meet the criteria
                        for the number of replicates and sequencing depth.
        ============================================================================
        */

        // Validate the project
        VALIDATION(ch_fastq, 'libraries', params.validation_rep, params.validation_depth)

        // Prepare the projects results channel for the summary channel.
        VALIDATION.out.projects
        .map{ item ->
            [item.project, item]
        }
        .set{ ch_validation_projects }

        // Add the validation information to the summary channel.
        VALIDATION.out.files
            .map{ meta, file -> [meta.id, meta, file]}
            .combine(pipeline_summary, by:0)
            .map { id, meta_lib, file, meta_sum ->
                [meta_sum.project, meta_sum + [
                    depth: meta_lib.depth,
                    depth_validity: meta_lib.depth_validity,
                    replicates_validity: meta_lib.replicates_validity
                ]]
            }
            // Modify the summary channel to be at the subproject level (Now a
            // sample may appear more than once if it belongs to multiple
            // subprojects).
            .combine(ch_validation_projects, by:0)
            .map { lib_sum, lib_meta, projects_sum ->
                lib_meta + [
                    group       : projects_sum.group,
                    group_id    : projects_sum.group_id,
                    num_valid_samples: projects_sum.num_valid_samples,
                    num_notvalid_samples: projects_sum.num_notvalid_samples,
                    group_validity: projects_sum.validity
                ]
            }
            .set { pipeline_summary }
        
        // Select only the valid libraries
        VALIDATION.out.files
            .filter { meta, file -> meta.depth_validity == 'valid' && meta.replicates_validity == 'valid' }
            .set{ ch_fastq }

        /*
        ============================================================================
            SUBWORKFLOW: Remove sequences that are not of interest
        ============================================================================
        */

        if (!params.skip_filt_db) {

            // Remove sequences that are not of interest (rRNA, tRNA, etc.)
            FILTERING_DB(ch_fastq, params.filtering_db_mismatches, "database", params.databases.filtering.rnacentral_non_miRNAs.fasta)

            // Update ch_fastq channel
            FILTERING_DB.out.unaligned.set{ch_fastq}

            // Add the filtering_db data to the pipeline_summary channel
            FILTERING_DB.out.unaligned
                .map{ meta, file -> [meta.id, meta, file]}
                .set{ filt_db_files_ch }
                
            pipeline_summary
                .map{ item -> [item.sample, item]}
                .groupTuple(by:0)
                .join(filt_db_files_ch, remainder:true)
                .flatMap { item ->

                    // Check if the last element of item is null
                    def lastElement = item.last()

                    // Add the filtering summary data to the summary channel
                    def additionalFields = (lastElement == null) ? [
                        filtering_db_total: 'NA',
                        filtering_db_only_align: 'NA',
                        filtering_db_failed: 'NA',
                    ] : [
                        filtering_db_total: item[2].filtering_db_total,
                        filtering_db_only_align: item[2].filtering_db_only_align,
                        filtering_db_failed: item[2].filtering_db_failed,
                    ]
                    def updatedMeta = item[1].collect { elem ->
                        elem + additionalFields
                    }

                    updatedMeta
                }
                .set{ pipeline_summary }
        }

        /*
        ============================================================================
            SUBWORKFLOW: Remove sequences that do not align with the genome
        ============================================================================
        */


        if (!params.skip_filt_genome) {
            
            // Remove those sequences that do not align with the reference genome
            FILTERING_GENOME(ch_fastq, params.filtering_genome_mismatches, "genome", null)

            // Update ch_fastq channel
            FILTERING_GENOME.out.aligned.set{ch_fastq}

            // Add the filtering_genome data to the pipeline_summary channel
            FILTERING_GENOME.out.aligned
                .map{ meta, file -> [meta.id, meta, file]}
                .set{ filt_genome_files_ch }

            pipeline_summary
                .map{ item -> [item.sample, item]}
                .groupTuple(by:0)
                .join(filt_genome_files_ch, remainder:true)
                .flatMap { item ->

                    // Check if the last element of item is null
                    def lastElement = item.last()

                    // Add the filtering summary data to the summary channel
                    def additionalFields = (lastElement == null) ? [
                        filtering_genome_total: 'NA',
                        filtering_genome_only_align: 'NA',
                        filtering_genome_failed: 'NA',
                    ] : [
                        filtering_genome_total: item[2].filtering_genome_total,
                        filtering_genome_only_align: item[2].filtering_genome_only_align,
                        filtering_genome_failed: item[2].filtering_genome_failed,
                    ]
                    def updatedMeta = item[1].collect { elem ->
                        elem + additionalFields
                    }

                    updatedMeta
                }
                .set{ pipeline_summary }
        }

        // Do not run this step when only pre-processing is to be done.
        if (!params.only_preprocessing){

            /*
            ============================================================================
                SUBWORKFLOW: Quantification of small RNA sequences
            ============================================================================
            */

            // Remove those specific fields of the sample except for the ID.
            ch_fastq
                .map { meta, file ->
                    // Filtrar los campos no deseados
                    def filteredMeta = meta.findAll { key, value ->
                        !(key.startsWith('filtering_') || key in ['depth', 'depth_validity', 'replicates_validity'])
                    }
                    [filteredMeta, file]
                }
                .set { ch_fastq }
            
            // Create count matrix
            QUANTIFICATION(ch_fastq)
            
            // Add the quantification data to the pipeline_summary channel
            QUANTIFICATION.out.subproject
                .map { meta, file ->
                    id = file.getName().replaceFirst(/\.counts\.tsv$/, '')
                    return [id, meta, file]
                }
                .set{ quantification_subproject_ch }

            pipeline_summary
                .map{ item -> [item.group, item]}
                .groupTuple(by:0)
                .join(quantification_subproject_ch, remainder:true)
                .flatMap { item ->
                    // Check if the last element of item is null
                    def lastElement = item.last()

                    // Set filteringGenomeValue based on conditions
                    def quantification = (lastElement == null) ? "not-quantified" : "quantified"
                    
                    // Add the "Filtering_genome" value to each map
                    def updatedItem = item[1].collect { element ->
                        element + [quantification: quantification]
                    }
                    
                    return updatedItem
                }
                .set{pipeline_summary}

            // QUANTIFICATION(
            //     FILTERING.out.bowtie_unaligned,
            //     metadata_dir_ch,
            //     filt_min_ch_counts,
            //     filt_min_samples_ch,
            //     rpm_ch,
            //     project_table_bool_ch,
            //     project_avg_bool_ch,
            //     subproject_avg_bool_ch
            // )
            
            // Change the meta.id from project to subproject.
            QUANTIFICATION.out.subproject
                .map { meta, file ->
                    def updatedMeta = meta.clone()
                    updatedMeta.id = file.getName().replaceFirst(/\.counts\.tsv$/, '')
                    return [updatedMeta, file]
                }
                .set {ch_counts}
        }
        
    }

    // Do not run these steps when only pre-processing is to be done.
    if (!params.only_preprocessing){

        /*
        ============================================================================
            UBWORKFLOW: Differential Expression Analysis
        ============================================================================
        */
        
        // Perform exploratory and differential expression analyses.
        DIFFEXPANALYSIS(ch_counts, params.dea_alpha, params.min_counts, params.min_samples)
        
        // Add the EA data to the pipeline_summary channel
        DIFFEXPANALYSIS.out.easum
            .splitCsv( header: true, sep: '\t' )
            .map{ item -> [item.Group, item]}
            .set{ ea_summary_ch }
        
        pipeline_summary
            .map { item -> [item.group, item] }
            .groupTuple(by: 0)
            .join(ea_summary_ch, remainder: true)
            .flatMap { item ->
                // Get the required variables
                def pip_summary = item[1]
                def easum = item.last()

                // Check if there is any information about the subproject in
                // the exploratory analysis.
                def additionalFields = (easum == null) ? [
                    pc1 : 'NA',
                    pc2: 'NA',
                    pc3: 'NA',
                    pc4: 'NA',
                    pc5: 'NA',
                    pc6: 'NA',
                    'p-value(mww)': 'NA'
                ] : [
                    pc1: easum.PC1,
                    pc2: easum.PC2,
                    pc3: easum.PC3,
                    pc4: easum.PC4,
                    pc5: easum.PC5,
                    PC6: easum.PC6,
                    'p-value(mww)': easum.'P-value(MWW)'
                ]

                // Merge each element of pip_summary with additionalFields
                pip_summary.collect { summary ->
                    summary + additionalFields
                }
            }
            .set{pipeline_summary}

        // Add the DEA data to the pipeline_summary channel
        DIFFEXPANALYSIS.out.deasum
            .splitCsv( header: true, sep: '\t' )
            .flatMap { item ->
                def samples = item.Samples.split(',')
                samples.collect { sample ->
                    [sample, item + [sample: sample]]
                }
            }
            .set{ dea_summary_ch }

        pipeline_summary
            .map{ item -> [item.sample, item]}
            .join(dea_summary_ch, remainder:true)
            .map { item ->
                def pip_summary = item[1]
                def dea_summary = item[2]

                // Campos adicionales a añadir
                def additionalFields = dea_summary ? [
                    comparison_id : dea_summary.Group,
                    test: dea_summary.Test,
                    'padj<alpha': dea_summary.'Padj<0.05', // CAMMBIAR LO DE 0.05 POR ALPHA
                    total: dea_summary.Total,
                    coefficient: dea_summary.Coefficient,
                    contrast: dea_summary.Contrast,
                    contrast_coefficient: dea_summary.Contrast_coefficient
                ] : [
                    comparison_id: 'NA',
                    test: 'NA',
                    'padj<0.05': 'NA',
                    total: 'NA',
                    coefficient: 'NA',
                    contrast: 'NA',
                    contrast_coefficient: 'NA'
                ]

                // Combinar los campos originales del segundo elemento con los adicionales
                pip_summary + additionalFields
            }
            .set{ pipeline_summary }
        
        // Change the meta.id from project to file.
        DIFFEXPANALYSIS.out.sig
            .map { meta, file ->
                def updatedMeta = meta.clone()
                updatedMeta.id = file.getName().replaceFirst(/\.dea_sig\.tsv$/, '')
                return [updatedMeta, file]
            }
            .set { dea_sig_ch }


        // Execute the annotation step if params.skip_annotation is false.
        if(!params.skip_annotation){

            /*
            ============================================================================
                SUBWORKFLOW: miRNA Annotation
            ============================================================================
            */

            // Identify which differentially expressed sRNA sequences are miRNAs.
            ANNOTATION(
                dea_sig_ch,
                ID_RESOLUTION.out.mirbase_mature,
                ID_RESOLUTION.out.pmiren_mature,
                ID_RESOLUTION.out.srnaanno_mature,
                params.annotation_mismatches,
                DIFFEXPANALYSIS.out.easum,
                params.ea_p_value
            )

            // Get the number of sequences annotated as miRNAs
            ANNOTATION.out.annotmat
                .map { item ->
                    // Count the number of sequences annotated as miRNAs
                    def filePath = item[1]
                    def numSeq = filePath ? filePath.text.split('\n').drop(1).size() : 0
                    
                    // Get the comparison id
                    comparison_id = filePath.getName().replaceFirst(/\.annot_len\.tsv$/, '')
                    
                    [comparison_id, numSeq]
                }
                .set{ num_annot_miRNAs_ch }

            // Get the number of sequences annotated as miRNAs (filtered)
            ANNOTATION.out.annotfiltmat
                .map { item ->
                    // Count the number of sequences annotated as miRNAs
                    def filePath = item[1]
                    def numSeq = filePath ? filePath.text.split('\n').drop(1).size() : 0
                    
                    // Get the comparison id
                    comparison_id = filePath.getName().replaceFirst(/\.annot_filt\.tsv$/, '')
                    
                    [comparison_id, numSeq]
                }
                .set{ num_annot_miRNAs_filt_ch }

            // Join annotmat and annotfiltmat channels
            num_annot_miRNAs_ch
                .join(num_annot_miRNAs_filt_ch, remainder:true)
                .set { annot_summary_ch }

            // Join the information from the miRNA annotation with the information from family grouping
            ANNOTATION.out.fam_sum
                .splitCsv(sep: '\t' )
                .join(annot_summary_ch, remainder:true)
                .set { annot_and_group_summary_ch }

            // Add the Annotation data to the pipeline_summary channel
            pipeline_summary
                .map{ item -> [item.comparison_id, item]}
                .groupTuple(by:0)
                .join(annot_and_group_summary_ch, remainder:true)
                .flatMap { item ->

                    // Campos adicionales a añadir
                    def additionalFields = item[6] ? [
                        num_annotated_miRNA : item[5],
                        num_annotated_miRNA_filt: item[6],
                        num_miRNA_fam: item[2],
                        num_miRNA_fam_divergent_ExpPatern: item[3],
                        miRNA_fam_divergent_ExpPatern: item[4],
                    ] : [
                        num_annotated_miRNA : 'NA',
                        num_annotated_miRNA_filt: 'NA',
                        num_miRNA_fam: 'NA',
                        num_miRNA_fam_divergent_ExpPatern: 'NA',
                        miRNA_fam_divergent_ExpPatern: 'NA',
                    ]

                    def updatedItem = item[1].collect { element ->
                        element + additionalFields
                    }

                    return updatedItem
                }
                .set{pipeline_summary}
        }
        
    }

    // Specify which samples have been discarded at any stage of the pipeline.
    pipeline_summary
        .map { item ->
            
            // Lista de valores inválidos
            def invalidValues = ['NA', 'not-valid', 'not-quantified']
            
            // Indicar si se ha encontrado un valor inválido
            def stopProcessing = false
            
            // Recorrer las claves del mapa y verificar los valores
            item.each { key, value ->
                // Si encontramos un valor no válido, marcamos todos los siguientes como "NA"
                if (stopProcessing || invalidValues.contains(value)) {
                    item[key] = 'NA'
                    stopProcessing = true
                }
            }
            
            return item
        }
        .set{ pipeline_summary }
    
    pipeline_summary.view()





    // INFORME HTML.
    // def skip_filtering_db = false

    // Channel
    //     .of(
    //         [sample: 'SRR14182749', species: 'Glycine max', species_id: 'gma', project: 'PRJNA720229', input: 'Downloaded', Trimming: 'Trimmed', Depth: 17600922, Depth_validity: 'valid', Replicates_validity: 'valid', Subproject: 'PRJNA720229_3', Subproject_id: 3, Num_valid_samples: 6, Num_notvalid_samples: 12, Subproject_validity: 'valid', Filtering_db_total: 17600922, Filtering_db_only_align: '6699112 (38.06%)', Filtering_db_failed: '10901810 (61.94%)', Filtering_genome_total: 10901810, Filtering_genome_only_align: '9150456 (83.94%)', Filtering_genome_failed: '1751354 (16.06%)'],
    //         [sample: 'ERR4078852', species: 'Amaranthus hypochondriacus', species_id: 'ahp', project: 'PRJEB38055', input: 'Downloaded', Trimming: 'Trimmed', Depth: 318829, Depth_validity: 'not-valid', Replicates_validity: 'not-valid', Subproject: 'PRJEB38055_1', Subproject_id: 1, Num_valid_samples: 0, Num_notvalid_samples: 6, Subproject_validity: 'NA', Filtering_db_total: 'NA', Filtering_db_only_align: 'NA', Filtering_db_failed: 'NA', Filtering_genome_total: 'NA', Filtering_genome_only_align: 'NA', Filtering_genome_failed: 'NA']
    //     )
    //     .map { task ->
    //         [task]
    //     }
    //     .collect()
    //     .map { tasks ->
    //         def template = new File('/home/antonio/Escritorio/Repositorios_GitHub/miRPlan-nf/assets/report-template.html').text

    //         def sampleLevelRows = tasks.collect { task ->
    //             """
    //             <tr>
    //                 <td>${task.sample}</td>
    //                 <td>${task.species}</td>
    //                 <td>${task.species_id}</td>
    //                 <td>${task.Depth}</td>
    //                 <td>${task.Depth_validity}</td>
    //                 <td>${task.Replicates_validity}</td>
    //                 <td>${task.project}</td>
    //                 <td>${task.Subproject}</td>
    //             </tr>
    //             """
    //         }.join("\n")

    //         def subprojectLevelRows = tasks.collect { task ->
    //             """
    //             <tr>
    //                 <td>${task.species}</td>
    //                 <td>${task.species_id}</td>
    //                 <td>${task.project}</td>
    //                 <td>${task.Subproject}</td>
    //                 <td>${task.Num_valid_samples}</td>
    //                 <td>${task.Num_notvalid_samples}</td>
    //                 <td>${task.Subproject_validity}</td>
    //             </tr>
    //             """
    //         }.join("\n")

    //         def dbFilteringRows = tasks.collect { task ->
    //             """
    //             <tr>
    //                 <td>${task.sample}</td>
    //                 <td>${task.species}</td>
    //                 <td>${task.species_id}</td>
    //                 <td>${task.project}</td>
    //                 <td>${task.Subproject}</td>
    //                 <td>${task.Filtering_db_total}</td>
    //                 <td>${task.Filtering_db_only_align}</td>
    //                 <td>${task.Filtering_db_failed}</td>
    //             </tr>
    //             """
    //         }.join("\n")

    //         def genomeFilteringRows = tasks.collect { task ->
    //             """
    //             <tr>
    //                 <td>${task.sample}</td>
    //                 <td>${task.species}</td>
    //                 <td>${task.species_id}</td>
    //                 <td>${task.project}</td>
    //                 <td>${task.Subproject}</td>
    //                 <td>${task.Filtering_genome_total}</td>
    //                 <td>${task.Filtering_genome_only_align}</td>
    //                 <td>${task.Filtering_genome_failed}</td>
    //             </tr>
    //             """
    //         }.join("\n")

    //         def htmlContent = template
    //             .replace('$sample_level_validation', sampleLevelRows)
    //             .replace('$subproject_level_validation', subprojectLevelRows)
    //             .replace('$database_filtering', dbFilteringRows)
    //             .replace('$genome_filtering', genomeFilteringRows)
        
    //         // ESTA PARTE ES PARA MOSTRAR O NO PARTE DEL CONTENIDO DEL HTML. EN
    //         // ESTE CASO, SI SKIP_FILTERING_DB ES true, NO SE MOSTRARÁ LA PARTE
    //         // DEL FILTRADO. CAMBIAR ESTO PARA PONER EL PARAMETRO DE ENTRADA QUE
    //         // ESPECIFICA SI SE SKIPEA O NO.
    //         // Solo insertamos el bloque de filtering_rnacentral si skip_filtering_db es falso
    //         if (skip_filtering_db) {
    //             htmlContent = htmlContent.replace('$filtering_rnacentral', '')  // Eliminar sección de filtrado si es verdadero
    //         } else {
    //             htmlContent = htmlContent.replace('$filtering_rnacentral', 'Contenido de filtrado aquí')  // Mostrar contenido si es falso
    //         }
            
    //     }
    //     .subscribe { htmlContent ->
    //         new File("/home/antonio/Escritorio/Repositorios_GitHub/miRPlan-nf/assets/execution-report.html").text = htmlContent
    //     }




    /// AHORA TENGO QUE VER COMO HACER PARA ESTABLECER LOS NA CUANDO EN ALGUNO
    /// DE LOS PASOS NO HA PASADO LA MUESTA Y GESTIONAR AQUELLOS CASOS EN LOS
    // QUE SE SALTE PASOS DEL PIPELINE. YA QUE AL SALTARSE PASOS NO HABRA SUMMARY.
    /// LUEGO TENGO QUE VER COMO DEVOLVER UNA TABLA Y/O UN INFORME DE RESULTADOS.






    // // #########################################################################
    // // #################### Build miRNA_vs_stress matrix #######################
    // // #########################################################################

    // RESULTS_GENERATION(
    //     DIFFEXPANALYSIS.out.dea_sum,
    //     ANNOTATION.out.fam_annot.collect()
    // )

}