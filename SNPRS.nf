#! /usr/bin/env nextflow
nextflow.enable.dsl=2

// SNPRS Main Script
// Params are read in from command line or from nextflow.config and/or conf/profiles.config

params.runmode = ""

params.composite_read_dir = ""
params.composite_reads = file(params.composite_read_dir)


workflow{
    
    // Runmode: Composite Genome
    if (params.runmode == "composite"){
        composite_data = fetchCompositeData()
    }
}