#!/usr/bin/env nextflow

/*
========================================================================================
    miRPlan Nextflow Workflow
========================================================================================
    Github   :
    Contact  :
----------------------------------------------------------------------------------------
*/




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
    IMPORT LOCAL FUNCTIONS/MODULES/SUBWORKFLOWS
========================================================================================
*/


//
// FUNCTIONS
//

include { validateAndAssignGenome    } from "../subworkflows/local/utils_mirplan_pipeline"
include { notTsvFilesError           } from "../subworkflows/local/utils_mirplan_pipeline"
include { validateGroupInputUsage    } from "../subworkflows/local/utils_mirplan_pipeline"
include { validateAccessionList      } from "../subworkflows/local/utils_mirplan_pipeline"


//
// MODULES
//
include { DOWNLOADLIB       } from '../modules/local/downloadlib'
include { DIFFEXPANALYSIS   } from "../modules/local/diffexpanalysis"


//
// SUBWORKFLOWS
//

include { FASTQ_DOWNLOAD_PREFETCH_FASTERQDUMP_SRATOOLS                      } from "../subworkflows/nf-core/fastq_download_prefetch_fasterqdump_sratools"
include { ID_RESOLUTION                           } from "../subworkflows/local/idresolution"
include { QUALITY_CONTROL as QUALITY_CONTROL_RAW  } from "../subworkflows/local/qualitycontrol"
include { QUALITY_CONTROL as QUALITY_CONTROL_TRIM } from "../subworkflows/local/qualitycontrol"
include { VALIDATION                              } from "../subworkflows/local/validation"
include { FILTERING as FILTERING_DB               } from "../subworkflows/local/filtering"
include { FILTERING as FILTERING_GENOME           } from "../subworkflows/local/filtering"
include { QUANTIFICATION                          } from "../subworkflows/local/quantification"
include { ANNOTATION                              } from "../subworkflows/local/annotation"
include { RESULTS_GENERATION                      } from "../subworkflows/local/results"


/*
========================================================================================
    IMPORT NF-CORE FUNCTIONS/MODULES/SUBWORKFLOWS
========================================================================================
*/


//
// FUNCTIONS
//
include { samplesheetToList } from 'plugin/nf-schema'


//
// MODULES: Installed directly from nf-core/modules
//
include { FASTP } from '../modules/nf-core/fastp'


/*
========================================================================================
    WORKFLOW - miRNA GLOBAL ANALYSIS
========================================================================================
*/


workflow MIRPLAN {

    main:

    // Fastq files empty channels
    ch_fastq         = Channel.empty()
    pipeline_summary = Channel.empty()
    ch_counts        = Channel.empty()

    // Create a channel from input file using params.input
    Channel
        .fromList(samplesheetToList(params.input, "${projectDir}/assets/schema_input.json"))
        .set{ch_input}
    
    /*
    ============================================================================
        Assign a default genome to the input files if necessary.
    ============================================================================
    */
    
    ch_input
        // Assign a genome to each library depending on the species.
        .map { item -> validateAndAssignGenome(item) }
        .set {ch_input}

    /*
    ============================================================================
        Check the inputs related to counts matrices
    ============================================================================
    */

    // Check the correct usage of the --from_counts parameter.
    ch_input
        .map{item -> item[4]}
        .filter{ it =~ /.*\.tsv$/ }
        .toList()
        .map { files -> notTsvFilesError(files) }

    // Check the correct usage of Group input
    ch_input
        .map { it[5] }
        .toList()
        .map{ item -> validateGroupInputUsage(item) }
    
    /*
    ============================================================================
        Separate the fastq files from the accession list files
    ============================================================================
    */

    ch_input
        .map { item ->
                [[species: item[0], species_id: 'null', project: item[1], metadata: item[2], genome: item[3], single_end: true, group_id:item[5]], item[4]]
        }
        .set{ ch_input }
    
    ch_input
        .branch { meta, file ->
                
                // Accession lists. It ends with '.txt'
                acclist: file.toString().endsWith('.txt')
                    def meta_acclist_with_id = [id: meta.project] + meta
                    return [meta_acclist_with_id, file]

                // Sequencing libraries. They end with '.fastq', '.fastq.gz', '.fq', o '.fq.gz'
                fastq: file.toString() =~ /\.(fastq(\.gz)?|fq(\.gz)?)$/
                    def file_wo_extension = file.getName().replaceFirst(/\.(fastq(\.gz)?|fq(\.gz)?)$/, '')
                    def meta_fastq_with_id = [id: file_wo_extension] + meta
                    return [meta_fastq_with_id, file]

                // Counts matrices. They end with '.tsv'
                counts: file.toString().endsWith('.tsv')
                    def meta_counts_with_id = [id: meta.project] + meta
                    return[meta_counts_with_id, file]   
        }
        .set { ch_input_files }
    
    // // Check that the count matrices provided as input are valid.
    // if (params.from_counts){

    //     // Check the counts matrices input files
    //     VALIDATION(ch_input_files.counts, 'counts', params.validation_rep, 0)

    //     // Assign the output channel to the ch_counts channel.
    //     VALIDATION.out.files
    //         // Select the valid files
    //         .filter { tuple ->
    //             tuple[1] =~ /.*\.valid\.tsv$/
    //         }.set{ ch_counts }

    //     // Create the pipeline_summary channel using the meta. OJO. COMO LAS TABLAS DE CONTEOS NO TIENEN LOS PASOS PREVIOS TIENEN NA EN PIPELINE_SUMMARY. ESO HACE QUE LUEGO AL FINAL SE PONGA TODO NA. COMPROBAR
    //     VALIDATION.out.files
    //         // Get the sample names from the matrix header 
    //         .flatMap { meta, file ->
    //             // Define the new elements of the channel
    //             def newMeta = meta.findAll { it.key != 'id' }
    //             def fileHeader = file.withReader { it.readLine() }
    //             def samples = fileHeader?.split("\t").drop(1) 

    //             // Crear un nuevo elemento por cada SRR, asignándolo a 'sample'
    //             samples.collect { srr -> [sample: srr] + newMeta }
    //         }.set{ pipeline_summary }
    // }

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
        ch_input_files.acclist
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
        DOWNLOADLIB(ch_input_files.acclist)

        // Split the tuple created in download_libraries
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
        ch_input_files.fastq
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
            .concat(ch_input_files.fastq)
            .set{ch_fastq}

        ///////////////////////// ESTO ES LO BUENO ///////////////////////////////
        // ch_input_files.acclist
        //     .flatMap { meta, file ->

        //         // Leer el archivo línea por línea y generar una nueva tupla por cada fila
        //         return file.readLines().collect { line ->  
        //             def updatedMeta = meta.clone()
        //             updatedMeta.id = line
        //             tuple(updatedMeta, line)  
        //         }
        //     }
        //     .set { ch_samples_sra_id }

        // // Download the libraries
        // FASTQ_DOWNLOAD_PREFETCH_FASTERQDUMP_SRATOOLS(ch_samples_sra_id, [])
        
        // // Add the Input information to the summary channel
        // ch_input_files.fastq
        //     .map { meta, _file ->
        //         [meta.id, [sample: meta.id, species: meta.species, species_id: meta.species_id, project: meta.project, input: 'Local']]
        //     }
        //     .set{pipeline_summary}

        // FASTQ_DOWNLOAD_PREFETCH_FASTERQDUMP_SRATOOLS.out.reads
        //     .map { meta, _file ->
        //         [meta.id, [sample: meta.id, species: meta.species, species_id: meta.species_id, project: meta.project, input: 'Downloaded']]
        //     }
        //     .concat(pipeline_summary)
        //     .set{pipeline_summary}
        
        // // Combine the downloaded libraries with those provided by the user in the same channel.
        // FASTQ_DOWNLOAD_PREFETCH_FASTERQDUMP_SRATOOLS.out.reads
        //     .concat(ch_input_files.fastq)
        //     .set{ch_fastq}
        ////////////////////////////////////////////////////////////////////////

        /*
        ============================================================================
            SUBWORKFLOW: Perform quality control of RAW data
        ============================================================================
        */

        if (!params.skip_fastqc) {
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
            params.trimming_adapters,
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

        if (!(params.skip_fastqc || params.skip_qc_trim)) {
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
            .map { _id, meta_lib, _file, meta_sum ->
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
            .map { _lib_sum, lib_meta, projects_sum ->
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
            .filter { meta, _file -> meta.depth_validity == 'valid' && meta.replicates_validity == 'valid' }
            .set{ ch_fastq }

        /*
        ============================================================================
            SUBWORKFLOW: Remove sequences that are not of interest
        ============================================================================
        */

        if (!params.skip_filt_db) {

            // Remove sequences that are not of interest (rRNA, tRNA, etc.)
            FILTERING_DB(ch_fastq, params.filtering_db_mismatches, "database", params.filtering_db_file)

            // Update ch_fastq channel
            FILTERING_DB.out.unaligned.set{ ch_fastq }

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
                    def filteredMeta = meta.findAll { key, _value ->
                        !(key.startsWith('filtering_') || key in ['depth', 'depth_validity', 'replicates_validity'])
                    }
                    [filteredMeta, file]
                }
                .set { ch_fastq }
            
            // Create count matrix
            QUANTIFICATION(ch_fastq, 'raw')
            
            // // Add the quantification data to the pipeline_summary channel
            // QUANTIFICATION.out.group_matrix
            //     .map { meta, file ->
            //         def id = file.getName().replaceFirst(/\.raw\.tsv$/, '')
            //         return [id, meta, file]
            //     }
            //     .set{ quantification_subproject_ch }

            // pipeline_summary
            //     .map{ item -> [item.group, item]}
            //     .groupTuple(by:0)
            //     .join(quantification_subproject_ch, remainder:true)
            //     .flatMap { item ->
            //         // Check if the last element of item is null
            //         def lastElement = item.last()

            //         // Set filteringGenomeValue based on conditions
            //         def quantification = (lastElement == null) ? "not-quantified" : "quantified"
                    
            //         // Add the "Filtering_genome" value to each map
            //         def updatedItem = item[1].collect { element ->
            //             element + [quantification: quantification]
            //         }
                    
            //         return updatedItem
            //     }
            //     .set{pipeline_summary}
            
            // // Change the meta.id from project to subproject.
            // QUANTIFICATION.out.group_matrix
            //     .map { meta, file ->
            //         def updatedMeta = meta.clone()
            //         updatedMeta.id = file.getName().replaceFirst(/\.raw\.tsv$/, '')
            //         return [updatedMeta, file]
            //     }
            //     .set {ch_counts}
        }
        
    }

    // // Do not run these steps when only pre-processing is to be done.
    // if (!params.only_preprocessing){


    //     /*
    //     ============================================================================
    //         UBWORKFLOW: Differential Expression Analysis
    //     ============================================================================
    //     */
        
    //     // Perform exploratory and differential expression analyses.
    //     DIFFEXPANALYSIS(ch_counts, params.dea_alpha, params.min_counts, params.min_samples)
        
    //     // Add the EA data to the pipeline_summary channel
    //     DIFFEXPANALYSIS.out.easum
    //         .splitCsv( header: true, sep: '\t' )
    //         .map{ item -> [item.Group, item]}
    //         .set{ ea_summary_ch }
        
    //     pipeline_summary
    //         .map { item -> [item.group, item] }
    //         .groupTuple(by: 0)
    //         .join(ea_summary_ch, remainder: true)
    //         .flatMap { item ->
    //             // Get the required variables
    //             def pip_summary = item[1]
    //             def easum = item.last()

    //             // Check if there is any information about the subproject in
    //             // the exploratory analysis.
    //             def additionalFields = (easum == null) ? [
    //                 pc1 : 'NA',
    //                 pc2: 'NA',
    //                 pc3: 'NA',
    //                 pc4: 'NA',
    //                 pc5: 'NA',
    //                 pc6: 'NA',
    //                 'p-value(mww)': 'NA'
    //             ] : [
    //                 pc1: easum.PC1,
    //                 pc2: easum.PC2,
    //                 pc3: easum.PC3,
    //                 pc4: easum.PC4,
    //                 pc5: easum.PC5,
    //                 PC6: easum.PC6,
    //                 'p-value(mww)': easum.'P-value(MWW)'
    //             ]

    //             // Merge each element of pip_summary with additionalFields
    //             pip_summary.collect { summary ->
    //                 summary + additionalFields
    //             }
    //         }
    //         .set{pipeline_summary}

    //     // Add the DEA data to the pipeline_summary channel
    //     DIFFEXPANALYSIS.out.deasum
    //         .splitCsv( header: true, sep: '\t' )
    //         .flatMap { item ->
    //             def samples = item.Samples.split(',')
    //             samples.collect { sample ->
    //                 [sample, item + [sample: sample]]
    //             }
    //         }
    //         .set{ dea_summary_ch }

    //     pipeline_summary
    //         .map{ item -> [item.sample, item]}
    //         .join(dea_summary_ch, remainder:true)
    //         .map { item ->
    //             def pip_summary = item[1]
    //             def dea_summary = item[2]

    //             // Campos adicionales a añadir
    //             def additionalFields = dea_summary ? [
    //                 comparison_id : dea_summary.Group,
    //                 test: dea_summary.Test,
    //                 'padj<alpha': dea_summary.'Padj<0.05', // CAMMBIAR LO DE 0.05 POR ALPHA
    //                 total: dea_summary.Total,
    //                 coefficient: dea_summary.Coefficient,
    //                 contrast: dea_summary.Contrast,
    //                 contrast_coefficient: dea_summary.Contrast_coefficient
    //             ] : [
    //                 comparison_id: 'NA',
    //                 test: 'NA',
    //                 'padj<0.05': 'NA',
    //                 total: 'NA',
    //                 coefficient: 'NA',
    //                 contrast: 'NA',
    //                 contrast_coefficient: 'NA'
    //             ]

    //             // Combinar los campos originales del segundo elemento con los adicionales
    //             pip_summary + additionalFields
    //         }
    //         .set{ pipeline_summary }
        
    //     // Change the meta.id from project to file.
    //     DIFFEXPANALYSIS.out.sig
    //         .map { meta, file ->
    //             def updatedMeta = meta.clone()
    //             updatedMeta.id = file.getName().replaceFirst(/\.dea_sig\.tsv$/, '')
    //             return [updatedMeta, file]
    //         }
    //         .set { dea_sig_ch }

    //     // Execute the annotation step if params.skip_annotation is false.
    //     if(!params.skip_annotation){

    //         /*
    //         ============================================================================
    //             SUBWORKFLOW: miRNA Annotation
    //         ============================================================================
    //         */

    //         // Identify which differentially expressed sRNA sequences are miRNAs.
    //         ANNOTATION(
    //             dea_sig_ch,
    //             params.annotation_mirbase,
    //             params.annotation_mirbase_taxon,
    //             params.annotation_srnaanno,
    //             params.annotation_pmiren,
    //             params.annotation_mismatches,
    //             DIFFEXPANALYSIS.out.easum,
    //             params.ea_p_value,
    //             params.annotation_min_db
    //         )

    //         ANNOTATION.out.annot_sum
    //             .splitCsv(sep: '\t', skip: 1)
    //             .map { [it[1], *it[2..-1]] }
    //             .set{annot_summary_ch}

    //         // Join the general summary and the length summary
    //         ANNOTATION.out.annot_sumlen
    //             .splitCsv(sep: '\t', skip: 1)
    //             .map { [it[1], *it[2..-1]] }
    //             .join(annot_summary_ch, remainder:true)
    //             .set {annot_summary_ch}

    //         // Join the information from the miRNA annotation with the information from family grouping
    //         ANNOTATION.out.fam_sum
    //             .splitCsv(sep: '\t' )
    //             .join(annot_summary_ch, remainder:true)
    //             .set { annot_and_group_summary_ch }

    //         // Add the Annotation data to the pipeline_summary channel
    //         pipeline_summary
    //             .map{ item -> [item.comparison_id, item]}
    //             .groupTuple(by:0)
    //             .join(annot_and_group_summary_ch, remainder:true)
    //             .flatMap { item ->

    //                 // Campos adicionales a añadir
    //                 def additionalFields = item[6] ? [
    //                     num_annotated_miRNA : item[17],
    //                     num_annotated_miRNA_filt: item[18],
    //                     num_annot_20nt: item[5],
    //                     num_annot_21nt: item[6],
    //                     num_annot_22nt: item[7],
    //                     num_annot_23nt: item[8],
    //                     num_annot_24nt: item[9],
    //                     num_annot_25nt: item[10],
    //                     num_annot_20nt_filt: item[11],
    //                     num_annot_21nt_filt: item[12],
    //                     num_annot_22nt_filt: item[13],
    //                     num_annot_23nt_filt: item[14],
    //                     num_annot_24nt_filt: item[15],
    //                     num_annot_25nt_filt: item[16],
    //                     num_miRNA_fam: item[2],
    //                     num_miRNA_fam_divergent_ExpPatern: item[3],
    //                     miRNA_fam_divergent_ExpPatern: item[4],
    //                 ] : [
    //                     num_annotated_miRNA : 'NA',
    //                     num_annotated_miRNA_filt: 'NA',
    //                     num_annot_20nt: 'NA',
    //                     num_annot_21nt: 'NA',
    //                     num_annot_22nt: 'NA',
    //                     num_annot_23nt: 'NA',
    //                     num_annot_24nt: 'NA',
    //                     num_annot_25nt: 'NA',
    //                     num_annot_20nt_filt: 'NA',
    //                     num_annot_21nt_filt: 'NA',
    //                     num_annot_22nt_filt: 'NA',
    //                     num_annot_23nt_filt: 'NA',
    //                     num_annot_24nt_filt: 'NA',
    //                     num_annot_25nt_filt: 'NA',
    //                     num_miRNA_fam: 'NA',
    //                     num_miRNA_fam_divergent_ExpPatern: 'NA',
    //                     miRNA_fam_divergent_ExpPatern: 'NA',
    //                 ]

    //                 def updatedItem = item[1].collect { element ->
    //                     element + additionalFields
    //                 }

    //                 return updatedItem
    //             }
    //             .set{pipeline_summary}
    //     }
        
    // }

    // // Specify which samples have been discarded at any stage of the pipeline.
    // pipeline_summary
    //     .map { item ->
            
    //         // Lista de valores inválidos
    //         def invalidValues = ['NA', 'not-valid', 'not-quantified']
            
    //         // Indicar si se ha encontrado un valor inválido
    //         def stopProcessing = false
            
    //         // Recorrer las claves del mapa y verificar los valores
    //         item.each { key, value ->
    //             // Si encontramos un valor no válido, marcamos todos los siguientes como "NA"
    //             if (stopProcessing || invalidValues.contains(value)) {
    //                 item[key] = 'NA'
    //                 stopProcessing = true
    //             }
    //         }
            
    //         return item
    //     }
    //     .set{ pipeline_summary }





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