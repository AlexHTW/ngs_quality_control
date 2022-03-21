import fastqc_extract
import data_preparation
import feature_engineering

training_data_dir = 'Documents/ngs_quality_control/training_data_test' #change to original path
exports_dir = 'Documents/ngs_quality_control/exported_datasets2' # change to original path

ngs_reads = fastqc_extract.import_all_reads(training_data_dir, exports_dir)

ngs_reads = data_preparation.prepare_fastqc_data(ngs_reads)

ngs_reads = feature_engineering.apply_feature_engineering(ngs_reads, exports_dir)

#ngs_reads.to_html('Documents/ngs_quality_control/dataset_check.html')