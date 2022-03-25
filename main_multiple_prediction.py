from src import fastqc_extract, data_preparation, feature_engineering, ml_model
#import sys

data_dir = 'test_data/all' #sys.argv[1]
model_dir = 'models' #sys.argv[2]
exports_dir = 'exports_evaluation_multiple' #sys.argv[3]
organism = 'complete_data' #sys.argv[4]
if organism not in ['Ecoli', 'Efcm', 'Sau']:
    organism = 'complete_data'

ngs_reads = fastqc_extract.import_all_reads(data_dir, exports_dir, include_metadata=False)
filenames = ngs_reads.index
ngs_reads = data_preparation.prepare_fastqc_data(ngs_reads)
ngs_reads = feature_engineering.apply_feature_engineering(ngs_reads, exports_dir)

pred = ml_model.predict_evaluation(ngs_reads, model_dir+'/model_rf_'+organism+'.pkl')

with open(exports_dir+'/predictions.txt', 'w') as f:
    for id, prediction in enumerate(pred):
        f.write(filenames[id]+': ')
        if prediction == 0:
            print(filenames[id]+ ': ugly')
            f.write('ugly\n')
        elif prediction == 1:
            print(filenames[id]+ ': good')
            f.write('good\n')
