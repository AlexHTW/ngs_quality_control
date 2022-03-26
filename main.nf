#!/usr/bin/env nextflow

input_path = projectDir + params.fastq_file
model_path = projectDir + params.model_dir
export_path = projectDir + params.export_dir

process callPythonScript {

    input:
    file fastq from input_path
    path model from model_path
    path export from export_path

    output:
    file "$export/predictions.txt" into evaluations

    script:
    """
    python3 $projectDir/main_single_prediction.py $fastq $model $export $params.organism
    """
}

evaluations.subscribe { println "Received: " + it.text}