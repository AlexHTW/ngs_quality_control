from src import fastqc_extract, data_preparation, feature_engineering, ml_model

training_data_dir = 'training_data'
exports_dir = 'exports'
model_dir = 'models'

ngs_reads = fastqc_extract.import_all_reads(training_data_dir, exports_dir)

ngs_reads = data_preparation.prepare_fastqc_data(ngs_reads)

ngs_reads = feature_engineering.apply_feature_engineering(ngs_reads, exports_dir)

ml_model.train_model(ngs_reads, model_dir, evaluate=True)
ml_model.train_model_per_organism(ngs_reads, model_dir, evaluate=True)