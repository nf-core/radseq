//
// This file holds several functions specific to the workflow/radseq.nf in the nf-core/radseq pipeline
//

class WorkflowRadseq {

    //
    // Check and validate parameters
    //
    public static void initialise(params, log) {
        //genomeExistsError(params, log) // function below

        if (!params.method) {
            log.error "type of workflow to execute not specified with e.g. '--method denovo' or via a detectable config file."
            System.exit(1)
        }
        if (params.method == 'reference') {
            if (!params.genome || params.genome == null) {
            log.error "need to specify a genome file with e.g. '--genome fasta' or via a detectable config file."
            System.exit(1)
            }
        }
        if (params.method == 'denovo'){
            if (!params.sequence_type) {
                log.error "need to specify the sequencing method with e.g. '--sequence_type' or via a detectable config file"
                System.exit(1)
            }
            if (!params.minreaddepth_withinindividual || params.minreaddepth_withinindividual == null) {
                log.warn("using default range of values for minReadDepth_withinIndividual")
            }
            if (params.method == 'denovo' && !params.minreaddepth_betweenindividual || params.minreaddepth_betweenindividual == null) {
                log.warn("using default range of values for minReadDepth_BetweenIndividual")
            }   
        }  
    }

    //
    // Get workflow summary for MultiQC
    //
    public static String paramsSummaryMultiqc(workflow, summary) {
        String summary_section = ''
        for (group in summary.keySet()) {
            def group_params = summary.get(group)  // This gets the parameters of that particular group
            if (group_params) {
                summary_section += "    <p style=\"font-size:110%\"><b>$group</b></p>\n"
                summary_section += "    <dl class=\"dl-horizontal\">\n"
                for (param in group_params.keySet()) {
                    summary_section += "        <dt>$param</dt><dd><samp>${group_params.get(param) ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>\n"
                }
                summary_section += "    </dl>\n"
            }
        }

        String yaml_file_text  = "id: '${workflow.manifest.name.replace('/','-')}-summary'\n"
        yaml_file_text        += "description: ' - this information is collected when the pipeline is started.'\n"
        yaml_file_text        += "section_name: '${workflow.manifest.name} Workflow Summary'\n"
        yaml_file_text        += "section_href: 'https://github.com/${workflow.manifest.name}'\n"
        yaml_file_text        += "plot_type: 'html'\n"
        yaml_file_text        += "data: |\n"
        yaml_file_text        += "${summary_section}"
        return yaml_file_text
    }

    //
    // Exit pipeline if incorrect --genome key provided
    //
    private static void genomeExistsError(params, log) {
        if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
            log.error "=============================================================================\n" +
                "  Genome '${params.genome}' not found in any config files provided to the pipeline.\n" +
                "  Currently, the available genome keys are:\n" +
                "  ${params.genomes.keySet().join(", ")}\n" +
                "==================================================================================="
            System.exit(1)
        }
    }
}
