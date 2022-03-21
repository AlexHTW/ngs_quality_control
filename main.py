import fastqc_extract

training_data_dir = 'Documents/ngs_quality_control/training_data_test' #change to original path
exports_dir = 'Documents/ngs_quality_control/exported_datasets2' # change to original path

ngs_reads = fastqc_extract.import_all_reads(training_data_dir, exports_dir)
print(ngs_reads.shape)
print(ngs_reads.columns)
