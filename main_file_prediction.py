from cmath import sin
import fastqc_extract
import data_preparation
import feature_engineering
import ml_model
import sys

fastq_file =  sys.argv[1] #'test_data/all/180525-18-3598-UW18720_S3_L001_R1_001.fastq.gz'
model_dir = sys.argv[2] #'models'
exports_dir = sys.argv[3] #'exports_evaluation_single'
organism = sys.argv[4] #complete_data
if organism not in ['Ecoli', 'Efcm', 'Sau']: #replace with dynamic list of available models
    organism = 'complete_data'
 
ngs_reads = fastqc_extract.import_all_reads(fastq_file, exports_dir, force_reimport=False, include_metadata=False, single_file=True)
filenames = ngs_reads.index
ngs_reads = data_preparation.prepare_fastqc_data(ngs_reads)
ngs_reads = feature_engineering.apply_feature_engineering(ngs_reads, exports_dir, force_reimport=True)

pred = ml_model.predict_evaluation(ngs_reads, model_dir+'/model_rf_'+organism+'.pkl')

with open(exports_dir+'/predictions.txt', 'w') as f:
    for id, prediction in enumerate(pred):
        f.write(filenames[id]+': ')
        if prediction == 0:
            print('ugly')
            f.write('ugly\n')
        elif prediction == 1:
            print('good')
            f.write('good\n')
