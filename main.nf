#!/usr/bin/env nextflow

quality_classification_input_path = projectDir + params.quality_classification_fastq_file
quality_classification_model_path = projectDir + params.quality_classification_model_dir
quality_classification_export_path = projectDir + params.quality_classification_export_dir

process qualityEvaluationClassification {

    input:
    file quality_classification_fastq from quality_classification_input_path
    path quality_classification_model from quality_classification_model_path
    path quality_classification_export from quality_classification_export_path

    output:
    file "$quality_classification_export/predictions.txt" into quality_evaluation

    script:
    """
    python3 $projectDir/main_single_prediction.py $quality_classification_fastq $quality_classification_model $quality_classification_export $params.quality_classification_organism
    """
}

quality_evaluation.subscribe { println "Received: " + it.text}