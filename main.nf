#!/usr/bin/env nextflow

input_path = projectDir + params.fastq_file
model_path = projectDir + params.model_dir
export_path = projectDir + params.export_dir


process callPythonScript {

    output:
    stdout result

    """
    python3 $projectDir/main_file_prediction.py $input_path $model_path $export_path $params.organism
    """
}

result.view { it.trim() }