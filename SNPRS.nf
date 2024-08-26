#! /usr/bin/env nextflow
nextflow.enable.dsl=2

// SNPRS Main Script
// Params are read in from command line or from nextflow.config and/or conf/profiles.config

params.runmode = ""
params.composite_read_dir = ""


workflow{
    
    // Runmode: Composite Genome
    if (params.runmode == "composite"){
        println "Running in composite mode"
    }
}