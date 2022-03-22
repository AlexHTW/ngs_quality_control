import fastqc_extract
import data_preparation
import feature_engineering
import ml_model

data_dir = 'Documents/ngs_quality_control/classification_data'
model_dir = 'Documents/ngs_quality_control/exported_datasets2'
exports_dir = 'Documents/ngs_quality_control/export_prediction_data'
organism = 'Sau' #if unknown use 'complete_data'

ngs_reads = fastqc_extract.import_all_reads(data_dir, exports_dir, include_metadata=False)
filenames = ngs_reads.index
ngs_reads = data_preparation.prepare_fastqc_data(ngs_reads)
ngs_reads = feature_engineering.apply_feature_engineering(ngs_reads, exports_dir)

pred = ml_model.predict_evaluation(ngs_reads, model_dir+'/model_rf_'+organism+'.pkl')

for id, prediction in enumerate(pred):
    with open(exports_dir+'/predictions.txt', 'w') as f:
        f.write(filenames[id]+': ')
        if prediction == 0:
            f.write('ugly')
        elif prediction == 1:
            f.write('good')
