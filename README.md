# Quality Evaluation of Next Generation Sequencing Data

This project aims to provide an automated quality assessment of NGS-data in fastq files using machine learning.  
For feature engineering, model training and evaluation a dataset of 184 fastq files was provided by the Robert Koch-Institut in Berlin, containing sequencing data of three different organisms (E. faecium, E. coli, S. aureus) with their respective manual quality evaluation (good/ugly), totaling 35GB of data.  

The feature and label extraction combines 3 approaches:  
1. The folder structure of the original dataset provides information regarding the organism, technology and evaluation.
2. The raw fastq data is used to calculate features based on per read-position accumulations.
3. Results of the FastQC analysis are used to extract and engineer most of the features, regarding basic statistics and results of the different [FastQC modules](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/).

Different machine learning algorithms have been tested and the most reliable results could be achieved using Random Forests.  

The process for data extraction, feature engineering, data exploration and model training can be followed in the jupyter notebook provided in the **/notebooks** directory.

### Training

The complete pipeline for model training can be executed using **main_training.py**.
Although the original data is not provided in this repository it can be run as is using intermediate result exports for every major step, located in the */exports* directory.  
The final models are saved in the */models* directory where one model is created for each organism type and one additional model for the complete dataset (which for the provided data is equally competent in evaluating a fastq file as the specialized model).

### Evaluation

The evaluation pipeline can be executed for a collection of fastq files using the **main_multiple_prediction.py** script. This will run FastQC, apply the data transformations and return evaluations for all fastq files in the given directory. Arguments for a specific organism (to use that specific model) and a disired export path can also be given.

The evaluation process for a single fastq file is prototyped using Nextflow. The arguments can be specified in *nextflow.config* and the **main_single_prediction.py** script is invoked in the **main.nf** script. It can be executed with *nextflow run main.nf* with in the repository or *nextflow run ngs_quality_control* from outside.
The results are saved in the specified project directories.
