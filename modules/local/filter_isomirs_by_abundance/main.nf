process FILTER_ISOMIRS_BY_ABUNDANCE {

    tag "$meta.id"

    // Conda or container environment
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/a1/a125c778baf3865331101a104b60d249ee15fe1dca13bdafd888926cc5490a34/data' :
        'community.wave.seqera.io/library/gawk:5.3.1--e09efb5dfc4b8156' }"

    input:
    tuple val(meta), path(gff_files)
    val filter_type
    val min_value
    val min_samples

    output:
    tuple val(meta), path("*.filtered.gff3"), emit: filt_gff3
    path "summary.tsv"                      , emit: sum

    script:
    """
    set -euo pipefail

    ############################################################################
    ##                     Create file with isomiRs data                      ##  
    ############################################################################

    # Output TSV file
    tsv="isomirs_data.tsv"
    tsv_filt="isomirs_data_filt.tsv"

    # Create header
    echo -e "sample\tsequence\traw\trpm\$([[ "${filter_type}" == "Relative_abundance" ]] && echo -e '\tmirna_ref\trpm_ref\trel_abundance')" > "\$tsv"

    # Extract data from each GFF file
    for gff in ${gff_files.join(' ')}; do
        
        # Skip "filtered" files to avoid errors
        [[ "\$gff" == *.filtered.gff3 ]] && continue

        # Get the file name
        sample=\$(basename "\$gff" .gff3)

        # Create isomiRs_data TSV file
        gawk -v sample="\$sample" -v filter_type="${filter_type}" '
        BEGIN { FS=OFS="\\t"; }

        /^#/ { next }

        {
            # Get the atrributes from the gff3 file
            split(\$9, fields, "; ");
            for (i in fields) {
            split(fields[i], kv, "=");
            attr[kv[1]] = kv[2];
            }

            # Attributes required
            seq = attr["Read"];
            raw = attr["Expression"];
            rpm = attr["Norm"];

            # Execute only if the filter type is Relative abundance
            if (filter_type == "Relative_abundance"){

                # Required attributes when filter_type = Relative_abundance
                class = attr["Class"];
                mirna_ref = attr["miRNA_seq"];

                # Save reference RPM values
                if (class == "ref_miRNA") {
                    rpm_ref_map[seq] = rpm;
                }

                # Create arrays with some attributes for END section
                row_count++;
                raw_rows[row_count] = sample "\\t" seq "\\t" raw "\\t" rpm "\\t" mirna_ref;
                classes[row_count] = class;
            
            # Execute when filter_type != Relative_abundance
            } else {
                # Print the required attributes for Raw or RPM filtering
                print sample, seq, raw, rpm;
            }
        }

        END {
            # Execute only if the filter type is Relative abundance
            if (filter_type == "Relative_abundance"){
                for (i = 1; i <= row_count; i++) {
                    split(raw_rows[i], data, "\\t");
                    sample = data[1]; seq = data[2]; raw = data[3]; rpm = data[4]; mirna_ref = data[5];

                    if (classes[i] == "ref_miRNA") {
                        rel_abund = 1;
                        rpm_ref = rpm;
                    } else {
                        rpm_ref = rpm_ref_map[mirna_ref] + 0;
                        rel_abund = (rpm_ref > 0) ? rpm / rpm_ref : "inf";
                    }

                    print sample, seq, raw, rpm, mirna_ref, rpm_ref, rel_abund;
                }
            }
        }
        ' "\$gff" >> "\$tsv"
    done

    ############################################################################
    ##                     Filter isomiRs data table                          ##  
    ############################################################################

    #  Execute only if the filter type is Relative abundance
    if [[ "${filter_type}" == "Relative_abundance" ]]; then

        # Select only the Maximum value of Relative abundance for each sample-sequence pair
        gawk -F"\\t" '
        NR == 1 { header = \$0; next }
        {
            key = \$1 FS \$2;
            val = \$7 + 0;
            if (!(key in max) || val > max[key]) {
                max[key] = val;
                row[key] = \$0;
            }
        }
        END {
            print header;
            for (k in row) print row[k];
        }
        ' "\$tsv" > "\$tsv_filt"

    # Execute when filter_type != Relative_abundance
    else
        # Get unique rows
        awk -F"\\t" '
        NR == 1 { print; next }
        !seen[\$1 FS \$2]++' "\$tsv" > "\$tsv_filt"
    fi



    ############################################################################
    ##                       Select valid sequences                           ##  
    ############################################################################

    # Determine filter column
    case "${filter_type}" in
        Raw) col=3; reject="lowRawCounts" ;;
        RPM) col=4; reject="lowRPM" ;;
        Relative_abundance) col=7; reject="lowRelAbundance" ;;
        *) echo "Invalid filter type"; exit 1 ;;
    esac

    # Select the sequences considered as valid depending on the min_value and min_samples values
    valid_seqs=\$(gawk -v col="\$col" -v min_val="${min_value}" -v min_s="${min_samples}" -F'\\t' '
    NR == 1 { next }

    {
        seq = \$2;
        sample = \$1;
        val = \$col + 0;

        # Solo contamos una vez por muestra si cumple el umbral
        if (val >= min_val && val != "NA" && val != "inf") {
            key = seq SUBSEP sample;
            if (!(key in seen)) {
                seen[key] = 1;
                count[seq]++;
            }
        }
    }

    END {
        for (seq in count)
            if (count[seq] >= min_s)
                print seq;
    }
    ' "\$tsv")

    # Save the valid sequences into a temporary file
    valid_tmp="__valid.tmp"
    echo "\$valid_seqs" > "\$valid_tmp"

    ############################################################################
    ##                     Filter isomiRs by abundance                        ##  
    ############################################################################

    # Iterate over each GFF file and filter based on the valid sequences
    for gff in ${gff_files.join(' ')}; do

        # Skip "filtered" files to avoid errors
        [[ "\$gff" == *.filtered.gff3 ]] && continue

        # Get the sample name from the GFF file
        sample=\$(basename "\$gff" .gff3)

        # Create the output GFF file name
        out_gff="\${gff%.gff3}.filtered.gff3"

        gawk -v FS="\\t" -v OFS="\\t" \
            -v sample="\$sample" \
            -v tsv="\$tsv" \
            -v val="${min_value}" \
            -v reject="\$reject" \
            -v filter="${filter_type}" \
            -v valid_file="\$valid_tmp" '
        BEGIN {
            # Read valid sequences into an array
            while ((getline line < valid_file) > 0) {
                valid[line] = 1;
            }

            # If filter is Relative_abundance, read the TSV file
            if (filter == "Relative_abundance") {
                while ((getline tsv_line < tsv) > 0) {
                    if (tsv_line ~ /^sample/) continue;           # Skip header
                    split(tsv_line, fields, "\\t");
                    key = fields[1] "|" fields[2] "|" fields[5];  # sample|sequence|mirna_ref
                    rel_ab[key] = fields[7];                      # rel_abundance
                }
            }
        }

        /^#/ { print; next }  # Imprimir cabeceras sin modificar

        {
            attr = \$9;

            # Get the sequence and miRNA reference from attributes
            match(attr, /Read=([^;]+)/, a);        seq = a[1];
            match(attr, /miRNA_seq=([^;]+)/, b);   mirna_ref = b[1];
            match(attr, /Class=([^;]+)/, c);       class = c[1];

            # If the filter is Relative_abundance, get the RPM value
            if (filter == "Relative_abundance") {
                key_rel = sample "|" seq "|" mirna_ref;
                rel_val = (key_rel in rel_ab) ? rel_ab[key_rel] : "NA";
                if (attr ~ /RelAbundance=/) {
                    attr = gensub(/RelAbundance=[^;]+/, "RelAbundance=" rel_val, 1, attr);
                } else {
                    attr = attr ";RelAbundance=" rel_val;
                }
            }

            # Apply Filter:
            #  - If class is ref_miRNA: PASS
            #  - Otherwise: PASS if seq is in valid, else REJECT:<reason>
            if (class == "ref_miRNA") {
                if (attr ~ /Filter=/)
                    attr = gensub(/Filter=[^;]+/, "Filter=PASS", 1, attr);
                else
                    attr = attr ";Filter=PASS";
            } else if (seq in valid) {
                if (attr ~ /Filter=/)
                    attr = gensub(/Filter=[^;]+/, "Filter=PASS", 1, attr);
                else
                    attr = attr ";Filter=PASS";
            } else {
                reason = "REJECT:" reject;
                if (attr ~ /Filter=/)
                    attr = gensub(/Filter=[^;]+/, "Filter=" reason, 1, attr);
                else
                    attr = attr ";Filter=" reason;
            }
                
            \$9 = attr;
            print;
        }
        ' "\$gff" > "\$out_gff"

    done

    # Remove temporary file
    rm -f "\$valid_tmp"

    ############################################################################
    ##                 Create summary file for each GFF3 file                 ##  
    ############################################################################

    # Output file with combined summary
    final_summary="summary.tsv"
    first=1

    # Process each GFF3 file
    for gff_file in *.filtered.gff3; do
        sample="\${gff_file%.filtered.gff3}"

        echo "Processing \$gff_file → \$sample.summary.tsv"

        gawk -F'\\t' -v sample="\$sample" -v first="\$first" '
            BEGIN {
                split("ref_miRNA iso_5p iso_3p iso_add3p iso_add5p iso_snv_seed iso_snv_central_offset iso_snv_central iso_snv_central_supp iso_snv mixed mixed_shift undefined", classes, " ");
            }
            # Extract Read-Class pairs from valid lines
            /^#/ { next }
            \$9 ~ /Filter=PASS/ {
                match(\$9, /Read=([^;]+)/, r)
                match(\$9, /Class=([^;]+)/, c)
                if (r[1] && c[1]) {
                    key = r[1] "\t" c[1]
                    if (!(key in seen)) {
                        seen[key] = 1
                        class_counts[c[1]]++
                    }
                }
            }
            END {
                if (first == 1) {
                    printf("sample\tnum_isomirs")
                    for (i = 1; i <= length(classes); i++)
                        printf("\\t%s", classes[i])
                    printf("\\n")
                }

                total = 0
                for (i = 1; i <= length(classes); i++) {
                    cls = classes[i]
                    val = (cls in class_counts) ? class_counts[cls] : 0
                    if (cls != "ref_miRNA")
                        total += val
                    output[i] = val
                }

                printf("%s\t%d", sample, total)
                for (i = 1; i <= length(classes); i++)
                    printf("\\t%d", output[i])
                printf("\\n")
            }
        ' "\$gff_file" >> "\$final_summary"

        first=0
    done

    """
}