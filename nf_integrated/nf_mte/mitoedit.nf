#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.bystander = "$projectDir/NO_FILE"
params.maf_to_alt_script = params.maf_to_alt_script ?: 'maf_to_alt_ref.py'

process PREP_MAF_ALT {
    publishDir "${params.outdir}/prep", mode: 'copy', overwrite: true

    input:
    path input_tab
    path ref_fasta

    output:
    path "prep_output_trimmed.tsv"
    path "prep_output_alt_genome.fa"
    path "ref_seq.txt"
    path "alt_seq.txt"

    script:
    """
    set -euo pipefail

    # Run maf_to_alt_ref.py to create trimmed table and alt_genome.fa
    python3 ${params.maf_to_alt_script} --maf ${input_tab} --fasta ${ref_fasta} --out_prefix prep_output

    # Convert original reference fasta to single-line sequence text (no header, unwrapped)
    awk '/^>/{next} {printf "%s", \$0} END{printf "\\n"}' ${ref_fasta} > ref_seq.txt

    # Convert alt genome fasta produced by script to single-line text (first record only)
    awk '/^>/{next} {printf "%s", \$0} END{printf "\\n"}' prep_output_alt_genome.fa > alt_seq.txt
    """
}

process MITOEDIT {
    errorStrategy 'ignore'  // Skip failed MITOEDIT runs

    publishDir "${params.outdir}/mitoedit_results", mode: 'copy', overwrite: true

    input:
    tuple val(row_data), val(position), val(base), val(test_base)
    path mtdna_seq
    path bystander_file

    output:
    tuple val(row_data), val(test_base), path("talen_output_${position}_${test_base}.txt"), optional: true

    script:
    def bystander = bystander_file.name != "NO_FILE" ? "--bystander_file ${bystander_file}" : ""
    """
    set +e
    mitoedit \
        --mtdna_seq_path ${mtdna_seq} \
        ${position} \
        ${test_base} \
        --output_prefix tmp_output \
        --min_spacer ${params.min_spacer} \
        --max_spacer ${params.max_spacer} \
        --array_min ${params.array_min} \
        --array_max ${params.array_max} \
        --filter ${params.mte_filter} \
        --cut_pos ${params.cut_pos} \
        ${bystander} 2>/dev/null

    # Rename output to include position and base for tracking
    if [ -f tmp_output/talen_output.txt ]; then
        mv tmp_output/talen_output.txt talen_output_${position}_${test_base}.txt
        rm -r tmp_output
    fi

    # Omit errors since most will fail
    exit 0
    """
}

process COMBINE_RESULTS {
    publishDir params.outdir, mode: 'copy'

    input:
    path input_tab
    path talen_files

    output:
    path "${params.output_prefix}_combined.tab"

    script:
    """
    #!/usr/bin/env python3
    import pandas as pd
    import glob
    import os
    import sys

    # Read original input (trimmed table produced by PREP_MAF_ALT)
    input_df = pd.read_csv("${input_tab}", sep='\\t')

    # Collect all MITOEDIT talen_output results staged into the work dir
    result_files = sorted(glob.glob('talen_output_*.txt'))

    results = []
    for talen_file in result_files:
        # Extract position and base from filename
        basename = os.path.basename(talen_file)
        # Expected format: talen_output_{position}_{base}.txt
        parts = basename.replace('.txt', '').split('_')
        # guard for unexpected filenames
        try {
            position = int(parts[2])
            test_base = parts[3]
        } catch (Exception) {
            continue
        }

        # Read MITOEDIT talen_output (assuming tab-separated)
        try {
            talen_df = pd.read_csv(talen_file, sep='\\t')
        } catch (Exception) {
            continue
        }

        # Find corresponding row in input (match on mtDNA_pos or Pos column)
        if 'mtDNA_pos' in input_df.columns {
            matches = input_df[input_df['mtDNA_pos'] == position]
        } else {
            matches = input_df[input_df['Pos'] == position]
        }

        if matches.empty {
            continue
        }

        input_row = matches.iloc[0].to_dict()

        # For each row in talen_output, duplicate the input row data
        for idx, talen_row in talen_df.iterrows() {
            combined_row = input_row.copy()
            combined_row['Tested_Base'] = test_base
            # Add all columns from talen_output
            for col, val in talen_row.items() {
                combined_row[col] = val
            }
            results.append(combined_row)
    }

    # Combine all results
    if results {
        combined_df = pd.DataFrame(results)
        # Reorder columns: input columns, Tested_Base, then talen columns
        input_cols = list(input_df.columns)
        talen_cols = [c for c in combined_df.columns if c not in input_cols and c != 'Tested_Base']
        final_cols = input_cols + ['Tested_Base'] + talen_cols
        combined_df = combined_df[final_cols]
        combined_df.to_csv("${params.output_prefix}_combined.tab", sep='\\t', index=False)
    } else {
        # Create empty output if no results
        pd.DataFrame().to_csv("${params.output_prefix}_combined.tab", sep='\\t', index=False)
    """
}

workflow {
    // Prepare trimmed table and alt genome + single-line sequences
    prep = PREP_MAF_ALT( Channel.fromPath(params.input_tab, checkIfExists: true),
                         Channel.fromPath(params.mtdna_seq_path, checkIfExists: true) )

    // Use the trimmed table from PREP for downstream splitting
    rows_ch = prep.out.collectFile().map { files ->
        // locate the trimmed TSV within the emitted files
        def trimmed = files.find { it.name == 'prep_output_trimmed.tsv' }
        return trimmed
    }.flatMap { trimmed_path ->
        Channel.fromPath(trimmed_path)
               .splitCsv(header: true, sep: '\t')
               .flatMap { row ->
                    def pos = row['mtDNA_pos'] ?: row['Pos']
                    def alt_allele = row['Ref. Allele'] ?: row['Alt_Allele'] ?: row['Ref_Allele']
                    def bases = ['A', 'C', 'T', 'G']
                    bases.findAll { it != alt_allele }.collect { test_base ->
                        [row, pos, alt_allele, test_base]
                    }
               }
    }

    // Create file channels for reference and alt single-line seqs
    ref_seq_ch = prep.out.collectFile().map { files -> files.find { it.name == 'ref_seq.txt' } }
    alt_seq_ch = prep.out.collectFile().map { files -> files.find { it.name == 'alt_seq.txt' } }

    bystander_ch = file(params.bystander_file)

    // Run MITOEDIT on reference genome
    mitoedit_ref = MITOEDIT(rows_ch, ref_seq_ch.first(), bystander_ch)

    // Run MITOEDIT on alt genome
    mitoedit_alt = MITOEDIT(rows_ch, alt_seq_ch.first(), bystander_ch)

    // Merge MITOEDIT outputs and pass to COMBINE_RESULTS
    all_talen_files = Channel.merge(
        mitoedit_ref.map { it[2] },
        mitoedit_alt.map { it[2] }
    ).collect()

    COMBINE_RESULTS(
        prep.out.collectFile().map { files -> files.find { it.name == 'prep_output_trimmed.tsv' } }.first(),
        all_talen_files
    )
}