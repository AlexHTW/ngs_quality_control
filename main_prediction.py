import fastqc_extract
import data_preparation
import feature_engineering
import ml_model

data_dir = 'Documents/ngs_quality_control/test_data/all'
model_dir = 'Documents/ngs_quality_control/exports'
exports_dir = 'Documents/ngs_quality_control/exports_testing'
organism = 'complete_data' #if unknown use 'complete_data'

ngs_reads = fastqc_extract.import_all_reads(data_dir, exports_dir, include_metadata=False)
filenames = ngs_reads.index
ngs_reads = data_preparation.prepare_fastqc_data(ngs_reads)
ngs_reads = feature_engineering.apply_feature_engineering(ngs_reads, exports_dir)

pred = ml_model.predict_evaluation(ngs_reads, model_dir+'/model_rf_'+organism+'.pkl')

#for id, prediction in enumerate(pred):
with open(exports_dir+'/predictions.txt', 'w') as f:
    for id, prediction in enumerate(pred):
        f.write(filenames[id]+': ')
        if prediction == 0:
            f.write('ugly\n')
        elif prediction == 1:
            f.write('good\n')
