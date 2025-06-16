/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Loaded from modules/local/
//

include { TSV_TO_FASTA              } from "../../modules/local/tsv_to_fasta"
include { DIFFEXPANALYSIS           } from "../../modules/local/diffexpanalysis"
include { ANNOTATE_DEA_RESULTS      } from "../../modules/local/annotate_dea_results"
include { BUILD_MIRNA_EVENT_MATRIX  } from "../../modules/local/build_mirna_event_matrix"
include { CONCAT_UNIQUE_GFF3        } from "../../modules/local/concat_unique_gff3"
include { RENAME_FILE_BY_ID         } from "../../modules/local/rename_file_by_id"

//
// SUBWORKFLOW: Loaded from subworkflows/local/
//

include { QUALITY_CONTROL as QUALITY_CONTROL_RAW  } from "../../subworkflows/local/qualitycontrol"
include { QUALITY_CONTROL as QUALITY_CONTROL_TRIM } from "../../subworkflows/local/qualitycontrol"
include { VALIDATION                              } from "../../subworkflows/local/validation"
include { FILTERING as FILTERING_DB               } from "../../subworkflows/local/filtering"
include { FILTERING as FILTERING_GENOME           } from "../../subworkflows/local/filtering"
include { QUANTIFICATION                          } from "../../subworkflows/local/quantification"
include { MIRNOTE as ANNOTATION                   } from "../../subworkflows/local/mirnote"

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//

include { FASTP } from '../../modules/nf-core/fastp'

//
// SUBWORKFLOW: Consisting entirely of nf-core/modules
//

include { FASTQ_DOWNLOAD_PREFETCH_FASTERQDUMP_SRATOOLS } from "../../subworkflows/nf-core/fastq_download_prefetch_fasterqdump_sratools"

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { samplesheetToList          } from 'plugin/nf-schema'
include { validateAndAssignGenome    } from "../../subworkflows/local/utils_mirplan_pipeline"
include { notTsvFilesError           } from "../../subworkflows/local/utils_mirplan_pipeline"
include { validateGroupInputUsage    } from "../../subworkflows/local/utils_mirplan_pipeline"
include { validateAccessionList      } from "../../subworkflows/local/utils_mirplan_pipeline"
include { filterByMwwPvalue          } from "../../subworkflows/local/utils_mirplan_pipeline"
include { writeSampleSheet           } from "../../subworkflows/local/utils_mirplan_pipeline"

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow MIRPLAN {

    take:
    samplesheet             // string: "path/to/sample_sheet.csv"
    ch_versions             // channel: [ path(versions.yml) ]

    main:

    // Empty channels
    ch_fastq            = Channel.empty()
    ch_pipeline_summary = Channel.empty()
    ch_counts           = Channel.empty()
    ch_versions         = Channel.empty()

    // Create a channel from input file using params.input
    Channel
        .fromList(samplesheetToList(samplesheet, "${projectDir}/assets/schema_input.json"))
        .set{ch_input}
    
    /*
    ============================================================================
       1. Assign a default genome to the input files if necessary.
    ============================================================================
    */
    
    ch_input
        // Assign a genome to each library depending on the species.
        .map { item -> validateAndAssignGenome(item) }
        .set {ch_input}

    /*
    ============================================================================
        2. Check the inputs related to counts matrices
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
        3. Separate the fastq files from the accession list files
    ============================================================================
    */

    ch_input
        .map { item ->
                [[species: item[0], project: item[1], metadata: item[2], genome: item[3], single_end: true, group_id:item[5]], item[4]]
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
    
    /*
    ============================================================================
        4. Validate input counts matrices
    ============================================================================
    */

    // Check that the count matrices provided as input are valid.
    if (params.from_counts){

        // Check the counts matrices input files
        VALIDATION(ch_input_files.counts, 'counts', params.validation_rep, 0)

        // Save the software version
        ch_versions = ch_versions.mix(VALIDATION.out.versions)

        // Assign the output channel to the ch_counts channel.
        VALIDATION.out.files
            // Select the valid files
            .filter { tuple ->
                tuple[1] =~ /.*\.valid\.tsv$/
            }.set{ ch_counts }

        // Create the ch_pipeline_summary channel using the meta.
        VALIDATION.out.files
            // Get the sample names from the matrix header 
            .flatMap { meta, file ->
                // Define the new elements of the channel
                def newMeta = meta.findAll { it.key != 'id' }
                def fileHeader = file.withReader { it.readLine() }
                def samples = fileHeader?.split("\t").drop(1) 

                // Crear un nuevo elemento por cada SRR, asignándolo a 'sample'
                samples.collect { srr -> [sample: srr] + newMeta }
            }.set{ ch_pipeline_summary }
            
        // Create FASTA file for the annotation step
        TSV_TO_FASTA(ch_counts, 0, 1, true, true)

        // Save the FASTA files from counts matrix columns in ch_fastq channel
        TSV_TO_FASTA.out.fasta
            .flatMap { meta, fasta_list ->
                fasta_list.collect { file_path ->
                    def file = file_path instanceof Path ? file_path : file(file_path)
                    def id = file.getBaseName()  // nombre del archivo sin extensión
                    def new_meta = meta.clone()
                    new_meta.id = id
                    [new_meta, file]
                }
            }
            .set { ch_fastq }
    }

    /*
    ============================================================================
        5. Pre-processing of the sequencing libraries
    ============================================================================
    */

    // If the input files are count matrices, do not execute the pre-processing.
    if (!params.from_counts) {

        /*
        ============================================================================
            5.1. SUBWORKFLOW: Download the study libraries
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
        
        ch_input_files.acclist
            .flatMap { meta, file ->

                // Leer el archivo línea por línea y generar una nueva tupla por cada fila
                return file.readLines().collect { line ->  
                    def updatedMeta = meta.clone()
                    updatedMeta.id = line
                    tuple(updatedMeta, line)  
                }
            }
            .set { ch_samples_sra_id }

        // Download the libraries
        FASTQ_DOWNLOAD_PREFETCH_FASTERQDUMP_SRATOOLS(ch_samples_sra_id, [])

        // Add the Input information to the summary channel
        ch_input_files.fastq
            .map { meta, _file ->
                [meta.id, [sample: meta.id, species: meta.species, project: meta.project, input: 'Local']]
            }
            .set{ch_pipeline_summary}

        FASTQ_DOWNLOAD_PREFETCH_FASTERQDUMP_SRATOOLS.out.reads
            .map { meta, _file ->
                [meta.id, [sample: meta.id, species: meta.species, project: meta.project, input: 'Downloaded']]
            }
            .concat(ch_pipeline_summary)
            .set{ch_pipeline_summary}
        
        // Combine the downloaded libraries with those provided by the user in the same channel.
        FASTQ_DOWNLOAD_PREFETCH_FASTERQDUMP_SRATOOLS.out.reads
            .concat(ch_input_files.fastq)
            .set{ch_fastq}

        /*
        ============================================================================
            5.2. SUBWORKFLOW: Perform quality control of RAW data
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
            5.3. SUBWORKFLOW: Perform trimming of the libraries using fastp.
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

        // Save the software version
        ch_versions = ch_versions.mix(FASTP.out.versions)

        // Add the trimming data to the ch_pipeline_summary channel
        FASTP.out.reads
            .map{meta, file ->
                [meta.id, meta, file]
            }
            .join(ch_pipeline_summary, remainder: true)
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
            .set{ch_pipeline_summary}
        
        // Change the meta.id from file to project.
        ch_fastq = FASTP.out.reads
            .map { meta, file ->
                def updatedMeta = meta.clone()
                updatedMeta.id = meta.project
                return [updatedMeta, file]
            }
            .groupTuple(by: 0, sort:true)

        /*
        ============================================================================
            5.4. SUBWORKFLOW: Perform quality control of RAW data
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
            5.5. SUBWORKFLOW: Validate the libraries (depth and replicates)
        ============================================================================
        */

        // Validate the project
        VALIDATION(ch_fastq, 'libraries', params.validation_rep, params.validation_depth)

        // Save the software version
        ch_versions = ch_versions.mix(VALIDATION.out.versions)

        // Get the group info from the metadata
        VALIDATION.out.files
            .map{ meta, _file ->  meta.metadata }
            .unique()
            .splitCsv( header: true, sep: '\t' )
            .map{ meta -> [meta.Run, meta]}
            .set { ch_metadata }
        
        // Add sample information to the validation info
        VALIDATION.out.files
            .map{ meta, file -> [meta.id, meta, file]}
            .combine(ch_metadata, by:0)
            .map{ id, meta, file, meta2 ->
                def new_meta = meta.clone()
                new_meta.group = "${meta.project}_${meta2.Group}"
                new_meta.group_id = meta2.Group
                new_meta.run = id
                [ "${meta.project}_${meta2.Group}",new_meta, file ]
            }
            .set { ch_validation_sample_info }

        // Add group information to validation info
        VALIDATION.out.projects
            .map{ item -> [item.group, item] }
            .set{ ch_validation_projects }
        
        // Prepare the validation info for the summary and ch_fastq channels
        ch_validation_sample_info
            .combine(ch_validation_projects, by:0)
            .map { group, meta1, file, meta2 ->
                def merged = meta1 + meta2
                if (meta2.validity == 'not-valid' && meta1.replicates_validity == 'valid') {
                    merged.replicates_validity = 'not-valid'
                }
                merged.validation_check = (meta2.validity == 'valid') ? 'OK' : 'FAIL'
                merged.remove('validity')
                // merged.remove('project')
                [meta1.run, merged, file]
            }
            .set{ ch_validation_group_info }

        // Validation info for the summary channel
        ch_validation_group_info
            .map{ run, meta, _file -> [run, meta]}
            .set{ ch_val_to_summary }
            
        // Add this info to summary channel
        ch_pipeline_summary
            .combine(ch_val_to_summary, by:0)
            .map{ id, meta1, meta2 -> 
                def merged = meta1 + meta2
                ['id', 'run', 'metadata', 'genome', 'single_end',
                'valid_groups', 'notvalid_groups'].each { merged.remove(it) }
                return merged
            }
            .set{ ch_pipeline_summary }

        // Validation info for the the quantification process
        ch_validation_group_info
            .filter { _id, meta, _file -> meta.validation_check == 'OK' }
            .map{ id, meta, file -> 
                def new_meta = meta.clone()
                ['group', 'run', 'num_valid_samples',
                'num_notvalid_samples',
                'validation_check',
                'group_id'].each { new_meta.remove(it) }
                [new_meta, file, meta.group]
            }
            .groupTuple(by: [0, 1])
            .map{ meta, file, groups ->
                def new_meta = meta.clone()
                new_meta.groups = groups
                [new_meta, file]
            }
            .set{ ch_fastq }
            
        /*
        ============================================================================
            5.6. SUBWORKFLOW: Remove sequences that are not of interest
        ============================================================================
        */

        if (params.filt_db) {

            // Remove sequences that are not of interest (rRNA, tRNA, etc.)
            FILTERING_DB(ch_fastq, "database", params.filt_db)

            // Save the software version
            ch_versions = ch_versions.mix(FILTERING_DB.out.versions)

            // Update ch_fastq channel
            ch_fastq = FILTERING_DB.out.unaligned

            // Add the filtering_db data to the ch_pipeline_summary channel
            FILTERING_DB.out.unaligned
                .map{ meta, file -> [meta.id, meta, file]}
                .set{ filt_db_files_ch }
            
            ch_pipeline_summary
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
                .set{ ch_pipeline_summary }
        }

        /*
        ============================================================================
            5.7. SUBWORKFLOW: Remove sequences that do not align with the genome
        ============================================================================
        */

        if (params.filt_genome) {
            
            // Remove those sequences that do not align with the reference genome
            FILTERING_GENOME(ch_fastq, "genome", null)

            // Save the software version
            ch_versions = ch_versions.mix(FILTERING_GENOME.out.versions)

            // Update ch_fastq channel
            ch_fastq = FILTERING_GENOME.out.aligned

            // Add the filtering_genome data to the ch_pipeline_summary channel
            FILTERING_GENOME.out.aligned
                .map{ meta, file -> [meta.id, meta, file]}
                .set{ filt_genome_files_ch }

            ch_pipeline_summary
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
                .set{ ch_pipeline_summary }
        }

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
        
        /*
        ============================================================================
            6. SUBWORKFLOW: Quantification of small RNA sequences
        ============================================================================
        */

        // Do not run this step when only pre-processing is to be done.
        if (!params.only_preprocessing){
            
            // Create count matrix
            QUANTIFICATION(ch_fastq, 'raw', params.counts_not_memory)

            // Save the software version
            ch_versions = ch_versions.mix(QUANTIFICATION.out.versions)
            
            // Add the quantification data to the ch_pipeline_summary channel
            QUANTIFICATION.out.group_matrix
                .map { meta, file ->
                    def id = file.getName().replaceFirst(/\.raw\.tsv$/, '')
                    return [id, meta, file]
                }
                .set{ quantification_subproject_ch }

            ch_pipeline_summary
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
                .set{ch_pipeline_summary}
            
            // Change the meta.id from project to subproject.
            QUANTIFICATION.out.group_matrix
                .map { meta, file ->
                    def updatedMeta = meta.clone()
                    updatedMeta.id = file.getName().replaceFirst(/\.raw\.tsv$/, '')
                    return [updatedMeta, file]
                }
                .set {ch_counts}
            
            // Create a samplesheet with intermediate results
            writeSampleSheet(
                ch_counts,
                "${params.outdir}/02-Results/02-Counts_matrix",
                "${params.outdir}/00-Additional_data/02-Samplesheets/Samplesheet_counts.csv"
            )

        }
    }

    // Do not run these steps when only pre-processing is to be done.
    if (!params.only_preprocessing){

        /*
        ============================================================================
            7. SUBWORKFLOW: Differential Expression Analysis
        ============================================================================
        */

        // Prepare the input channel for the differential expression analysis.
        ch_counts
            .map{ meta, file ->
                return [meta, file, meta.metadata, meta.group_id]
            }
            .set{ ch_counts_to_dea }

        // Perform exploratory and differential expression analyses.
        DIFFEXPANALYSIS(ch_counts_to_dea, params.dea_alpha, params.min_counts, params.min_samples)

        // Save the software version
        ch_versions = ch_versions.mix(DIFFEXPANALYSIS.out.versions)
        
        // Add the EA data to the ch_pipeline_summary channel
        DIFFEXPANALYSIS.out.easum
            .map{it -> it[1]}
            .splitCsv( header: true, sep: '\t' )
            .map{ item -> [item.Group, item]}
            .set{ ea_summary_ch }
        
        ch_pipeline_summary
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
            .set{ch_pipeline_summary}

        // Add the DEA data to the ch_pipeline_summary channel
        DIFFEXPANALYSIS.out.deasum
            .map{it -> it[1]}
            .splitCsv( header: true, sep: '\t' )
            .flatMap { item ->
                def samples = item.Samples.split(',')
                samples.collect { sample ->
                    [sample, item + [sample: sample]]
                }
            }
            .set{ dea_summary_ch }

        ch_pipeline_summary
            .map{ item -> [item.sample, item]}
            .join(dea_summary_ch, remainder:true)
            .map { item ->
                def pip_summary = item[1]
                def dea_summary = item[2]

                // Campos adicionales a añadir
                def additionalFields = dea_summary ? [
                    comparison_id : dea_summary.Group,
                    test: dea_summary.Test,
                    'padj<alpha': dea_summary.'Padj<alpha', // CAMMBIAR LO DE 0.05 POR ALPHA
                    total: dea_summary.Total,
                    coefficient: dea_summary.Coefficient,
                    contrast: dea_summary.Contrast,
                    contrast_coefficient: dea_summary.Contrast_coefficient
                ] : [
                    comparison_id: 'NA',
                    test: 'NA',
                    'padj<alpha': 'NA',
                    total: 'NA',
                    coefficient: 'NA',
                    contrast: 'NA',
                    contrast_coefficient: 'NA'
                ]

                // Combinar los campos originales del segundo elemento con los adicionales
                pip_summary + additionalFields
            }
            .map { record ->
                def new_record = record.clone()
                def isNumeric = { val -> 
                    try { val as Double; return true } catch (e) { return false }
                }

                def mww = record['p-value(mww)']
                new_record.ea_check = (isNumeric(mww) && (mww as Double) < params.ea_p_value) ? 'OK' : 'FAIL'

                def padj_str = record['padj<alpha'] as String
                new_record.dea_check = (padj_str != 'NA' && padj_str != '0') ? 'OK' : 'FAIL'

                return new_record
            }
            .set{ ch_pipeline_summary }
                    
        // Prepare the channel for the next steps of the workflow.
        DIFFEXPANALYSIS.out.deasum
            .map{it -> it[1]}
            .splitCsv( header: true, sep: '\t' )
            .map{ item -> [item.Group, item]}
            .set{ ch_dea_summ_to_sig }

        // Change the meta.id from project to file and remove EMPTY files
        DIFFEXPANALYSIS.out.sig
            .filter { _meta, file ->
                !file.getName().endsWith('dea_sig_EMPTY.tsv')
            }
            .map { meta, file ->
                def updatedMeta = meta.clone()
                updatedMeta.id = file.getName().replaceFirst(/\.dea_sig\.tsv$/, '')
                return [updatedMeta.id, updatedMeta, file]
            }
            .combine(ch_dea_summ_to_sig, by:0)
            .map{ _id, meta, file, dea_sig ->
                return [meta + [coefficient:dea_sig.Coefficient, samples:dea_sig.Samples], file]
            }
            .set { ch_dea_sig }

        /*
        ============================================================================
            8. SUBWORKFLOW: miRNA Annotation
        ============================================================================
        */

        // Execute the annotation step if params.skip_annotation is false.
        if(!params.skip_annotation){


            // Filter dea results using a mww p-value threshold
            if (params.ea_p_value != null) {

                // Prepare ea summary channel
                DIFFEXPANALYSIS.out.easum
                    .map{ meta, file -> [meta.id, meta, file] }
                    .set{ ch_easum }
                
                // Add the ea summay data to the dea data channel
                ch_dea_sig
                .map { meta, file ->
                    def key = "${meta.project}_${meta.group_id}"
                    return [key, meta, file]
                }
                .combine(ch_easum, by:0)
                .map { _project, dea_meta, dea_file, _ea_meta, ea_file ->
                    return [dea_meta, dea_file, ea_file]
                }
                .set{ ch_dea_ea_sig}

                // Filter the groups using the exploratory analysis results
                ch_dea_ea_sig = filterByMwwPvalue(ch_dea_ea_sig, params.ea_p_value)

                // Remove the ea data from the channel dea data channel
                ch_dea_ea_sig
                    .map{item -> [item[0], item[1]]}
                    .set{ch_dea_sig}
            }

            // Get the input species names
            ch_input
                .map{it[0].species}
                .unique()
                .collect()
                .set{ch_species_names}  

            // Only reference (canonical) miRNAs will be considered
            if (params.mirna_classes == 'ref_miRNA') {

                // Identify miRNA sequences
                ANNOTATION(
                    ch_fastq,
                    ch_species_names,
                    params.databases,
                    params.substitutions,
                    params.five_add,
                    params.three_add,
                    params.ends_modification,
                    true
                )

            // Other classes of miRNAs will be considered
            } else {

                // Identify miRNA sequences
                ANNOTATION(
                    ch_fastq,
                    ch_species_names,
                    params.databases,
                    params.substitutions,
                    params.five_add,
                    params.three_add,
                    params.ends_modification,
                    false
                )
            }

            // Save the software version
            ch_versions = ch_versions.mix(ANNOTATION.out.versions)

            // Prepare the summary channel
            ANNOTATION.out.summary
                .map{item -> [item.id, item]}
                .set{ch_annot_summary}

            // Add the Annotation data to the ch_pipeline_summary channel
            ch_pipeline_summary
                .map{ item -> [item.sample, item]}
                .groupTuple(by:0)
                .join(ch_annot_summary, remainder:true)
                .flatMap { item ->
                    // Additional fields to add
                    def additionalFields = item[2] ? [
                        species_db : item[2].species_db,
                        database : item[2].database,
                        num_pot_isomirs : item[2].num_pot_isomirs,
                        num_isomirs : item[2].num_isomirs,
                        ref_miRNA : item[2].ref_miRNA,
                        iso_5p : item[2].iso_5p,
                        iso_3p : item[2].iso_3p,
                        iso_add5p : item[2].iso_add5p,
                        iso_add3p : item[2].iso_add3p,
                        iso_snv_seed : item[2].iso_snv_seed,
                        iso_snv_central : item[2].iso_snv_central,
                        iso_snv_central_offset : item[2].iso_snv_central_offset,
                        iso_snv_central_supp : item[2].iso_snv_central_supp,
                        mixed : item[2].mixed,
                        mixed_shift : item[2].mixed_shift,
                        undefined : item[2].undefined
                    ] : [
                        species_db : 'NA',
                        database : 'NA',
                        num_pot_isomirs : 'NA',
                        num_isomirs : 'NA',
                        ref_miRNA : 'NA',
                        iso_5p : 'NA',
                        iso_3p : 'NA',
                        iso_add5p : 'NA',
                        iso_add3p : 'NA',
                        iso_snv_seed : 'NA',
                        iso_snv_central : 'NA',
                        iso_snv_central_offset : 'NA',
                        iso_snv_central_supp : 'NA',
                        mixed : 'NA',
                        mixed_shift : 'NA',
                        undefined : 'NA'
                    ]

                    def updatedItem = item[1].collect { element ->
                        element + additionalFields
                    }

                    return updatedItem
                }
                .map { record ->
                    def new_record = record.clone()
                    def num_pot_iso = record.num_pot_isomirs?.toString()
                    new_record.annotation_check = (num_pot_iso != 'NA' && num_pot_iso != '0') ? 'OK' : 'FAIL'
                    return new_record
                }
                .set{ch_pipeline_summary}
            
            // Prepare the input chennel for CONCAT_UNIQUE_GFF3 process
            if (!params.from_counts){

                // Create a single GFF3 file for each group with the annotated sequences
                ANNOTATION.out.annotation
                    .flatMap { tuple -> 
                        def meta = tuple[0]
                        def file = tuple[1]

                        // For each valid group, we generate a new map.
                        meta.valid_groups.collect { group_num ->
                            def new_meta = meta.clone()
                            new_meta.id = "${meta.project}_${group_num}"
                            return [[id:new_meta.id, genome:new_meta.genome], file]
                        }
                    }
                    .groupTuple(by: 0)
                    .set{ ch_gff3_by_group}
            } else {

                // Each counts matrix is a group
                ANNOTATION.out.annotation
                    .map { meta, file ->
                        def new_meta = meta.clone()
                        new_meta.id = "${meta.project}_${meta.group_id}"
                        new_meta.remove('run')
                        return [new_meta, file]
                    }
                    .groupTuple(by: 0)
                    .set{ ch_gff3_by_group }
            }

            // Concat gff3 files by group and remove duplicates
            CONCAT_UNIQUE_GFF3(ch_gff3_by_group)

            // Prepare annotation channel
            CONCAT_UNIQUE_GFF3.out.gff3
                .map{meta, file -> return[meta.id, meta, file]}
                .set{ch_mirna_annot}

            // Prepate dea_files channel
            ch_dea_ea_sig
                .map{meta, file, _ea_file ->
                    return["${meta.project}_${meta.group_id}", meta, file]}
                .set{ dea_ea_files }

            // Combine both channels
            dea_ea_files
                .combine(ch_mirna_annot, by: 0)
                .map { item ->
                    // Add the species id to the output meta
                    def new_meta = item[1].clone()
                    return [new_meta, item[2], item[4]]
                }
                .set { ch_group_miRNAs_input }

            // Add annotation to DEA results dataframe
            ANNOTATE_DEA_RESULTS(ch_group_miRNAs_input, params.mirna_classes)

            // Save the software version
            ch_versions = ch_versions.mix(ANNOTATE_DEA_RESULTS.out.versions)

            // Prepare the summary channel
            ANNOTATE_DEA_RESULTS.out.fam_sum
                .splitCsv(sep: '\t')
                .set{ch_annotate_dea_results_sum}
            
            // Add the Annotation data to the ch_pipeline_summary channel
            ch_pipeline_summary
                .map{ item -> [item.comparison_id, item]}
                .groupTuple(by:0)
                .join(ch_annotate_dea_results_sum, remainder:true)
                .flatMap { item ->

                    // Additional fields to add
                    def additionalFields = item[2] ? [
                        fam_members_same_pattern : item[2],
                        fam_members_diff_pattern : item[3],
                        names_members_diff_pattern : item[4]
                    ] : [
                        fam_members_same_pattern : 'NA',
                        fam_members_diff_pattern : 'NA',
                        names_members_diff_pattern : 'NA'
                    ]

                    def updatedItem = item[1].collect { element ->
                        element + additionalFields
                    }

                    return updatedItem
                }
                .map { record ->
                    def new_record = record.clone()
                    def fam = record.fam_members_same_pattern
                    new_record.dea_annotation_check = (fam != 'NA') ? 'OK' : 'FAIL'
                    return new_record
                }
                .set{ch_pipeline_summary}

            /*
            ============================================================================
                9. Create global matrices
            ============================================================================
            */

            // Create global matrices
            if (params.global_matrix){

                // Create a channel to rename metadata files
                ANNOTATE_DEA_RESULTS.out.unique
                    .filter { _meta, file ->
                        !file.getName().endsWith('_EMPTY.unique.tsv')
                    }
                    .map{ meta, _file -> [meta, meta.metadata]}
                    .set{ ch_to_rename_metadata }

                // Rename metadata files 
                RENAME_FILE_BY_ID(ch_to_rename_metadata)

                // Renamed metadata channel
                RENAME_FILE_BY_ID.out.renamed
                    .map{ meta, file -> [meta.id, meta, file]}
                    .set{ ch_renamed_metadata }

                // Add the new metadata to main channel
                ANNOTATE_DEA_RESULTS.out.unique
                    .filter { _meta, file ->
                        !file.getName().endsWith('_EMPTY.unique.tsv')
                    }
                    .map{ meta, file -> [meta.id, meta, file]}
                    .combine(ch_renamed_metadata, by:0)
                    .map{ _id, meta1, file, _meta2, metadata_file ->
                        [file, metadata_file, meta1.samples]
                    }
                    .collect()
                    .map { list_elements ->

                        // Agrupar en sublistas de 3 elementos
                        def groups = list_elements.collate(3)

                        // Extraer listas separadas por posición
                        def files = groups.collect { it[0] }
                        def metas = groups.collect { it[1] }
                        def samples = groups.collect { it[2] }

                        return [files, metas, samples]
                    }
                    .set{ ch_to_create_pa_matrix }

                // Create the both presence-absence and log2fc matrices  
                BUILD_MIRNA_EVENT_MATRIX(ch_to_create_pa_matrix, params.global_fields)

                // Save the software version
                ch_versions = ch_versions.mix(BUILD_MIRNA_EVENT_MATRIX.out.versions)
            }
        }        
    }

    // Move the check fields to the end
    ch_pipeline_summary
        .map { record ->
            def reordered = new LinkedHashMap()

            // Añadir primero las claves que NO terminan en '_check'
            record.keySet().findAll { !it.endsWith('_check') }.each { key ->
                reordered[key] = record[key]
            }

            // Luego añadir las claves que SÍ terminan en '_check'
            record.keySet().findAll { it.endsWith('_check') }.each { key ->
                reordered[key] = record[key]
            }

            return reordered
        }
        .set { ch_pipeline_summary }


    emit:
    summary        = ch_pipeline_summary                     // channel: [id:, sample:, etc]
    versions       = ch_versions                             // channel: [ path(versions.yml) ]



}