#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

log.info """\
    NEXTFLOW - DSL2 - LDSC - P I P E L I N E
    ==================================================
    run as                  : ${workflow.commandLine}
    run location            : ${launchDir}
    started at              : ${workflow.start}
    python exe              : ${params.my_python}
    container               : ${workflow.containerEngine}:${workflow.container}

    Tools
    ==================================================
    Python                  : ${params.my_python}
    PLINK                   : ${params.my_plink}

    Chromosomes
    ==================================================
    chromosome_list         : ${params.chromosome_list}

    Inputs
    ==================================================
    summary stats manifest  : ${params['input_sumstats_manifest_csv']}

    Parameters
    ==================================================
    plink clump p1          : ${params.clump_options.clump_p1}
    plink clump p2          : ${params.clump_options.clump_p2}
    plink clump r2          : ${params.clump_options.clump_r2}
    plink clump kb          : ${params.clump_options.clump_kb}

    Outputs
    ==================================================
""".stripIndent()

include { BIOFILTER_POSITIONS } from './biofilter_wrapper.nf'

workflow {
    input_sumstats_plink = PLINK_CLUMP_SETUP()
    output_clumps_loci = PLINK_CLUMP(input_sumstats_plink)
}

workflow PLINK_CLUMP_SETUP {
    main:
        sumstats_lines = new File(params['input_sumstats_manifest_csv']).readLines()
        info_row_lists = []
        sumstats_lines[1..-1].each {
            line ->
            /* groovylint-disable-next-line GStringExpressionWithinString */
            lineG = line.replace('${launchDir}', "${launchDir}")
            line_parts = lineG.toString().trim().replace('[', '').replace(']', '').split(',') as List
            info_row_lists.add(line_parts)
        }

        sumstats_channel = Channel.fromList(info_row_lists)

        sumstats_by_chr = sumstats_channel.combine(Channel.fromList(params['chromosome_list']))

        plink_suffixes = [
            '--bfile' : ['.bed', '.bim', '.fam'],
            '--pfile': ['.pgen', '.pvar', '.psam']
        ]

        sumstats_chr_plinkset = sumstats_by_chr.map {
            analysis, sumstats_file, plink_pre, plink_post, plink_flag, chr ->
            new Tuple(analysis, chr, sumstats_file, plink_flag, new Tuple(
                *plink_suffixes[plink_flag].collect { ext -> "${plink_pre}${chr.toString() == '23' ? 'X' : chr}${plink_post}" + ext })
            )
        }

    emit:
        sumstats_chr_plinkset
}

workflow PLINK_CLUMP_SETUP_STITCHED {
    take:
        group_pheno_sumstats
    main:
        group_info_lines = new File(params['group_info_manifest_csv']).readLines()
        
        sumstats_group_prefixes = [:]
        sumstats_group_suffixes = [:]
        sumstats_group_plink_flags = [:]

        group_info_lines[1..-1].each {
            line ->
            /* groovylint-disable-next-line GStringExpressionWithinString */
            lineG = line.replace('${launchDir}', "${launchDir}")
            line_parts = lineG.toString().trim().replace('[', '').replace(']', '').split(',') as List

            sumstats_group_prefixes[line_parts[0]] = line_parts[1]
            sumstats_group_suffixes[line_parts[0]] = line_parts[2]
            sumstats_group_plink_flags[line_parts[0]] = line_parts[3]
        }

        sumstats_by_chr = group_pheno_sumstats.combine(Channel.fromList(params['chromosome_list']))
        sumstats_by_chr_with_plink_info = sumstats_by_chr.map {
            group, pheno, sumstats_file, chr ->
            new Tuple("${group}_${pheno}", sumstats_file, 
                sumstats_group_prefixes[group], sumstats_group_suffixes[group], sumstats_group_plink_flags[group],
                chr
            )
        }.filter {
            it[2] != null
        }

        plink_suffixes = [
            '--bfile' : ['.bed', '.bim', '.fam'],
            '--pfile': ['.pgen', '.pvar', '.psam']
        ]

        sumstats_chr_plinkset = sumstats_by_chr_with_plink_info.map {
            analysis, sumstats_file, plink_pre, plink_post, plink_flag, chr ->
            new Tuple(analysis, chr, sumstats_file, plink_flag, new Tuple(
                *plink_suffixes[plink_flag].collect { ext -> "${plink_pre}${chr.toString() == '23' ? 'X' : chr}${plink_post}" + ext })
            )
        }

    emit:
        sumstats_chr_plinkset
}

workflow PLINK_CLUMP {
    take:
        // Tuple with (analysis, chr, sumstats_file, (plink_trio))
        sumstats_chr_plinkset
    main:
        // Rename the Channel from the setup sub-workflow
        make_clump_chr_input_in = sumstats_chr_plinkset

        // Name all the scripts we'll need
        input_cleaning_script = "${moduleDir}/scripts/make_clump_input_by_chr.py"
        merging_script = "${moduleDir}/scripts/merge_clump_results.py"
        plotting_script = "${moduleDir}/scripts/make_loci_plot.py"

        // Call the process to make chr-separated input for clumping
        (clump_chr_input, min_p_files) = make_clump_chr_input(
            make_clump_chr_input_in,
            params['clump_options'],
            params['input_colnames'],
            input_cleaning_script
        )

        // Parse the min p value files
        // This will return (analysis, chr, min p-value)
        clump_chr_min_p_channel = parse_analysis_chr_pval(min_p_files)

        // Filter on the minimum p-value
        // This is because Plink will throw an error
        // if there are no significant clumps
        // Filter (analysis, chr, min p-value) and map to just (analysis, chr) for join
        clump_calls_to_run = clump_chr_min_p_channel.filter {
            analysis, chr, pval -> pval < params.clump_options.clump_p1
            }.map { analysis, chr, pval -> new Tuple(analysis, chr) }

        // Take each minimum p-value row and merge them to one table
        // Write this output file for our records
        min_p_files_collect = min_p_files.map { analysis, chr, file -> file }.collect()
        min_p_table = make_min_p_table(min_p_files_collect, params.clump_options.clump_p1)

        // Join clump input with the filtered clump calls
        // Then call plink clump for correct set of processes
        // clump_chr_input is Tuple of:
        // (analysis, chr, sumstats_chr, extract_file, plink_flag, plink_set)
        clump_filtered_input = clump_chr_input.join(clump_calls_to_run, by: [0, 1])
        // clump_filtered_input.view()
        (clump_output, logs) = call_plink_clump(
            clump_filtered_input,
            params['clump_options']
        )

        // Group clump results and merge the outputs
        merge_input = clump_output.groupTuple(by: 0)
        (analysis_merged_clumps, analysis_merged_loci) = merge_clump_results(
            merge_input,
            merging_script
        )

        // Collect the clump and locus files into lists
        clump_files_collect = analysis_merged_clumps.map { analysis, merge_clumps -> merge_clumps }.collect()
        loci_files_collect = analysis_merged_loci.map { analysis, merge_clumps -> merge_clumps }.collect()

        analysis_sumstats_files_collect = sumstats_chr_plinkset.map { 
            analysis, chr, sumstats, plink_flag, plink_trio -> new Tuple(analysis, sumstats) 
            }.groupTuple(by: 1)
        
        final_clumps_with_sumstats = analysis_merged_clumps.join(analysis_sumstats_files_collect, by: 0)
        final_loci_with_sumstats = analysis_merged_loci.join(analysis_sumstats_files_collect, by: 0)

        // Collect the actual summary stats input used
        analysis_sumstats_files_collect = clump_filtered_input.map {
            analysis, chr, sumstats_file, extract_file, plink_flag, plink_set -> \
            new Tuple(analysis, sumstats_file)
            }.groupTuple(by: 0)

        // Make plots and summary tables
        if (params.annotate) {
            // Call biofilter annotation sub-workflow
            clump_lead_snp_positions = make_lead_snp_biofilter_input(clump_files_collect)
            biofilter_input = clump_lead_snp_positions.map { filename -> new Tuple('clump_lead_snps', filename) }
            annotations = BIOFILTER_POSITIONS(biofilter_input)

            // Compile summary table
            (clump_summary_file, loci_summary_file) = compile_results_with_annot(clump_files_collect, loci_files_collect, annotations)

            // Make plots using the final locus boundaries
            locus_plots = make_locus_plots_with_annot(final_loci_with_sumstats, annotations, plotting_script)

            // Also make plots using the clump groups
            clump_plots = make_clump_plots_with_annot(final_clumps_with_sumstats, annotations, plotting_script)
        }
        else {
            // Make summary table
            (clump_summary_file, loci_summary_file) = compile_results_no_annot(clump_files_collect, loci_files_collect)

            // Make plots using the final locus boundaries
            locus_plots = make_locus_plots_no_annot(final_loci_with_sumstats, plotting_script)

            // Also make plots using the clump groups
            clump_plots = make_clump_plots_no_annot(final_clumps_with_sumstats, plotting_script)
        }
}

process make_clump_chr_input {
    publishDir "${launchDir}/Clump/Input/"
    machineType 'n2-standard-4'

    input:
        tuple val(analysis), val(chr), path(sumstats), val(plink_flag), path(plink_set)
        val clump_options_map
        val colnames_map
        path make_clump_input_script
    output:
        tuple val(analysis), val(chr), path("${analysis}.${chr}.txt.gz"), path("${analysis}.${chr}.extract.txt"), val(plink_flag), path(plink_set)
        tuple val(analysis), val(chr), path("${analysis}.${chr}.min_p.csv")
    shell:
    """
        ${params.my_python} ${make_clump_input_script} \
            --chr ${chr} \
            --plink-bim ${plink_set[1]} ${plink_flag} \
            --analysis ${analysis} \
            --sumstats ${sumstats} \
            --chr-col "${colnames_map.chr_col}" \
            --pos-col ${colnames_map.pos_col} \
            --id-col  ${colnames_map.id_col} \
            --p-col   ${colnames_map.p_col} \
            --a1-col  ${colnames_map.A1_col} \
            --a2-col  ${colnames_map.A2_col} \
            --p-thresh ${clump_options_map.clump_p2}
    """
// No stub because min p files are required for workflow plan
}

process parse_analysis_chr_pval {
    machineType 'n2-standard-4'

    input:
        tuple val(analysis), val(chr), path(min_p_file) // this is the appropriate path object
    output:
        tuple val(analysis), val(chr), val(min_p)
    script:
        exist_checks = 0
        sleep_interval = 100 // milliseconds
        max_sleep = 10 * 1000 // 10 seconds

        pval_file_string = "${launchDir}/Clump/Input/${analysis}.${chr}.min_p.csv"

        // make sure the file exists
        current_sleep = 0
        while (exist_checks < 3 && current_sleep <= max_sleep) {
            if (file(pval_file_string).exists()) {
                exist_checks += 1
            }
            sleep(sleep_interval)
            current_sleep += sleep_interval
        }

        // make sure the file is not empty
        not_empty_checks = 0
        current_sleep = 0
        while (not_empty_checks < 3 && current_sleep <= max_sleep) {
            if ((new File(pval_file_string)).empty()) {
                not_empty_checks += 1
            }
            sleep(sleep_interval)
            current_sleep += sleep_interval
        }

        // assert that the file both exists and is not empty
        assert file(pval_file_string).exists() && !(new File(pval_file_string)).empty() : 'File System Latency Issues. Please Try Re-Running'

        // read in the p-value information for this summary stats / chr pair
        pval_lines = (new File(pval_file_string)).readLines()
        min_p = Float.parseFloat(pval_lines[2].split(',')[1])

        '''
        echo "done"
        '''
}

process make_min_p_table {
    publishDir "${launchDir}/Summary/"
    machineType 'n2-standard-4'

    input:
        path min_p_files
        val clump_lead_snp_min_p
    output:
        path "all_analysis_chr_min_pvals.csv"
    script:
        """
        #! ${params.my_python}
        import pandas as pd

        input_files = '${min_p_files.join(' ')}'.split()
        output_file = 'all_analysis_chr_min_pvals.csv'

        dfs = []
        for f in input_files:
            dfs.append(pd.read_csv(f, index_col=0, header=None))
        df = pd.concat(dfs, axis=1).transpose()
        df['Run Clump'] = pd.to_numeric(df['MIN_P']) < float(${clump_lead_snp_min_p})
        df = df.sort_values(by=['ANALYSIS', 'CHR'])
        print(df)
        df.to_csv(output_file, index=False)
        """
    stub:
        '''
        touch all_analysis_chr_min_pvals.csv
        '''
}

process call_plink_clump {
    publishDir "${launchDir}/Clump/Output/"
    machineType 'n2-standard-4'

    input:
        tuple val(analysis), val(chr), path(sumstats), path(extract_file), val(plink_flag), path(plink_set)
        val clump_options_map
    output:
        tuple val(analysis), val(chr), path("${analysis}.${chr}.clumps"), path(extract_file)
        tuple val(analysis), val(chr), path("${analysis}.${chr}.log")
    shell:
    """
        ${params.my_plink} \
            ${plink_flag} ${plink_set[0].toString().replace('.bed', '').replace('.pgen', '')} \
            --clump ${sumstats} ${plink_flag == '--pfile' ? '--clump-unphased' : ''} \
            --clump-p1 ${clump_options_map.clump_p1} \
            --clump-p2 ${clump_options_map.clump_p2} \
            --clump-r2 ${clump_options_map.clump_r2} \
            --clump-kb ${clump_options_map.clump_kb} \
            --extract range ${extract_file} \
            --out ${analysis}.${chr}
    """
    stub:
    """
        touch ${analysis}.${chr}.clumps
        touch ${analysis}.${chr}.log
    """
}

process merge_clump_results {
    publishDir "${launchDir}/Clump/Results/"
    machineType 'n2-standard-4'

    input:
        tuple val(analysis), val(chr_list), path(clump_file_list), path(extract_file_list)
        path(merge_script)
    output:
        tuple val(analysis), path("${analysis}.clumps.csv")
        tuple val(analysis), path("${analysis}.loci.csv")
    shell:
    """
        ${params.my_python} ${merge_script} \
         --clumps ${clump_file_list.join(' ')} \
         --extract-files ${extract_file_list.join(' ')} \
         --analysis ${analysis}
    """
    stub:
    """
        touch ${analysis}.clumps.csv
        touch ${analysis}.loci.csv
    """
}

process make_lead_snp_biofilter_input {
    publishDir "${launchDir}/Annotations/"
    machineType 'n2-standard-4'

    input:
        path(merged_clumps)
    output:
        path('clumped_lead_snps_biofilter_input_positions.txt')
    script:
    """
        #! ${params.my_python}
        import pandas as pd

        # List of input files
        clump_files = '${merged_clumps.join(' ')}'.split()

        # Iterate over input files and concatenate
        dfs = []
        for f in clump_files:
            dfs.append(pd.read_csv(f).dropna(subset=['#CHROM', 'Lead_SNP', 'POS'], how='any'))
        df = pd.concat(dfs)

        # Make sure position column is an int for sorting
        df['POS'] = df['POS'].astype(int)
        df = df[['#CHROM', 'Lead_SNP', 'POS']].drop_duplicates().sort_values(by=['#CHROM', 'POS'])

        # Write output
        print(df)
        df.to_csv('clumped_lead_snps_biofilter_input_positions.txt', index=False, header=False, sep=' ')
    """
    stub:
    """
        touch clumped_lead_snps_biofilter_input_positions.txt
    """
}

process compile_results_with_annot {
    publishDir "${launchDir}/Summary/", mode: 'copy', overwrite: true
    machineType 'n2-standard-4'

    input:
        path(merged_clumps)
        path(merged_loci)
        tuple val(biofilter_data_nickname), path(biofilter_annotations)
    output:
        path('all_clumps_with_annot.csv')
        path('all_loci_with_annot.csv')
    script:
    """
        #! ${params.my_python}
        import pandas as pd

        # Read in the biofilter annotations
        annot = pd.read_csv("${biofilter_annotations}", index_col='Var_ID')
        print(annot)

        # List of clump output
        clump_files = '${merged_clumps.join(' ')}'.split()

        # Iterate over clump output and merge
        dfs = []
        for f in clump_files:
            temp_df = pd.read_csv(f, index_col='Lead_SNP')
            dfs.append(temp_df)
        df = pd.concat(dfs)

        # Add RSID and Gene columns
        df[['Lead_SNP_Nearest_Gene', 'Lead_SNP_RSID']] = annot.reindex(df.index)[['Gene', 'RSID']].values
        df.index.name = 'Lead_SNP'
        print(df)

        # Write the output
        df.to_csv('all_clumps_with_annot.csv')

        # Now do loci
        loci_files = '${merged_loci.join(' ')}'.split()

        # Iterate over and concatenate locus files
        dfs = []
        for f in loci_files:
            temp_df = pd.read_csv(f, index_col='Lead_SNP')
            dfs.append(temp_df)
        df = pd.concat(dfs)

        # Add RSID and Gene columns
        df[['Lead_SNP_Nearest_Gene', 'Lead_SNP_RSID']] = annot.reindex(df.index)[['Gene', 'RSID']].values
        df.index.name = 'Lead_SNP'
        print(df)

        # Write the output
        df.to_csv('all_loci_with_annot.csv')
    """
    stub:
        '''
        touch all_clumps_with_annot.csv
        touch all_loci_with_annot.csv
        '''
}

process compile_results_no_annot {
    publishDir "${launchDir}/Summary/", mode: 'copy', overwrite: true
    machineType 'n2-standard-4'

    input:
        path(merged_clumps)
        path(merged_loci)
    output:
        path('all_clumps.csv')
        path('all_loci.csv')
    script:
        """
        #! ${params.my_python}
        import pandas as pd

        # List of clump files
        clump_files = '${merged_clumps.join(' ')}'.split()

        # Iterate over clump files and merge
        dfs = []
        for f in clump_files:
            temp_df = pd.read_csv(f, index_col='Lead_SNP')
            dfs.append(temp_df)
        df = pd.concat(dfs)
        df.index.name = 'Lead_SNP'

        # Write merged clump output
        print(df)
        df.to_csv('all_clumps.csv')

        # Now do loci
        loci_files = '${merged_loci.join(' ')}'.split()

        # Iterate over and concatenate loci files
        dfs = []
        for f in loci_files:
            temp_df = pd.read_csv(f, index_col='Lead_SNP')
            dfs.append(temp_df)
        df = pd.concat(dfs)
        df.index.name = 'Lead_SNP'

        # Write merged output of the loci
        print(df)
        df.to_csv('all_loci.csv')
        """
    stub:
        '''
        touch all_clumps.csv
        touch all_loci.csv
        '''
}

process make_locus_plots_with_annot {
    publishDir "${launchDir}/Plots/", mode: 'copy', overwrite: true
    machineType 'n2-standard-4'

    input:
        tuple val(analysis), path(loci_file), path(sumstats_files)
        tuple val(biofilter_data_nickname), path(biofilter_annotations)
        path plotting_script
    output:
        path("${analysis}.loci.png")
    shell:
        """
        ${params.my_python} ${plotting_script} \
          --analysis ${analysis} \
          --annot ${biofilter_annotations} \
          --sumstats ${sumstats_files.join(' ')} \
          --loci ${loci_file}
        """
    stub:
        """
        touch ${analysis}.loci.png
        """
}

process make_locus_plots_no_annot {
    publishDir "${launchDir}/Plots/", mode: 'copy', overwrite: true
    machineType 'n2-standard-4'

    input:
        tuple val(analysis), path(loci_file), path(sumstats_files)
        path plotting_script
    output:
        path("${analysis}.loci.png")
    shell:
        """
        ${params.my_python} ${plotting_script} \
          --analysis ${analysis} \
          --sumstats ${sumstats_files.join(' ')} \
          --loci ${loci_file}
        """
    stub:
        """
        touch ${analysis}.loci.png
        """
}

process make_clump_plots_with_annot {
    publishDir "${launchDir}/Plots/", mode: 'copy', overwrite: true
    machineType 'n2-standard-4'

    input:
        tuple val(analysis), path(clump_file), path(sumstats_files)
        tuple val(biofilter_data_nickname), path(biofilter_annotations)
        path plotting_script
    output:
        path("${analysis}.clumps.png")
    shell:
        """
        ${params.my_python} ${plotting_script} \
          --analysis ${analysis} \
          --annot ${biofilter_annotations} \
          --sumstats ${sumstats_files.join(' ')} \
          --loci ${clump_file}
        """
    stub:
        """
        touch ${analysis}.clumps.png
        """
}

process make_clump_plots_no_annot {
    publishDir "${launchDir}/Plots/", mode: 'copy', overwrite: true
    machineType 'n2-standard-4'

    input:
        tuple val(analysis), path(clump_file), path(sumstats_files)
        path plotting_script
    output:
        path("${analysis}.clumps.png")
    shell:
        """
        ${params.my_python} ${plotting_script} \
          --analysis ${analysis} \
          --sumstats ${sumstats_files.join(' ')} \
          --loci ${clump_file}
        """
    stub:
        """
        touch ${analysis}.clumps.png
        """
}