import os
import pandas as pd
import numpy as np
import pyfastx
from fastqcparser import FastQCParser

def run_fastqc(rootdir):
    for root, dirs, files in os.walk(rootdir):
        for filename in files:
            if filename.endswith("fastq.gz"):
                fastq_path = root+'/'+filename
                # create results folder
                result_path = root + "/fastqc_results"
                if not os.path.exists(result_path):
                    os.system("mkdir " + result_path)

                # check if fastqc results already exist
                result_file_path = root+'/fastqc_results/'+ filename.replace(".fastq.gz", "_fastqc")
                if not os.path.exists(result_file_path):
                    os.system("fastqc -o " + result_path + " --extract " + fastq_path)

def extract_metadata_from_export(outdir):
    metadata = pd.read_json(outdir+'/metadata.json')
    return metadata

def extract_metadata_from_source(rootdir, outdir, export):
    if not os.path.exists(outdir) and export:
        os.system("mkdir " + outdir)

    reads =	{
        "filename": [],
        "organism": [],
        "technology": [],
        "read_number": [],
        "evaluation": []
    }
    for root, dirs, files in os.walk(rootdir):
        if root.endswith("fastqc"):
            extracted_vars = root.split('/')
            reads["organism"].append(extracted_vars[-5])
            reads["evaluation"].append(extracted_vars[-4])
            reads["technology"].append(extracted_vars[-3])
            reads["filename"].append(extracted_vars[-1].replace("_fastqc", ""))
            read_number = extracted_vars[-1].split('_')[-3][1]
            reads["read_number"].append(read_number)
    metadata = pd.DataFrame(reads)
    if export:
        metadata.to_json(outdir+'/metadata.json')
    return metadata

def extract_metadata(rootdir, outdir, force_reimport=False, export=True):
    if os.path.exists(outdir+'/metadata.json') and not force_reimport:
        return extract_metadata_from_export(outdir)
    else:
        return extract_metadata_from_source(rootdir, outdir, export)

def phred_score(quality, offset=33):
    return ord(quality)-offset

def extract_fastq_positions(path, max_len_estim=1000):
    phred_sum = max_len_estim*[0] # sum all qualities per position, initialize with enough space for max_len
    counts = max_len_estim*[0] # count occurences of each position
    base_content = {'G': max_len_estim*[0], #count occurences of each base per position
                    'C': max_len_estim*[0],
                    'A': max_len_estim*[0], 
                    'T': max_len_estim*[0],
                    'N': max_len_estim*[0]}
    # fill lists with values from fastq file
    for name, seq, qual in pyfastx.Fastq(path, build_index=False):
        for i in range(len(seq)):
            phred_sum[i] += phred_score(qual[i])
            counts[i] += 1
            base_content[seq[i]][i] += 1
    
    phred_sum = np.array(phred_sum)
    counts = np.array(counts)
    max_len = np.min(np.where(counts==0)) #find index of first 0 value as maximum length of sequences
    phred_sum = phred_sum[0:max_len]
    counts = counts[0:max_len]
    means = phred_sum/counts
    
    # convert the lists in base_content to np.arrays to perform calculations
    for key in base_content:
        base_content[key] = np.array(base_content[key])
    max_len = np.min(np.where(base_content['G']==0)) #find index of first 0 value as maximum length of sequences
    for key in base_content:
        base_content[key] = base_content[key][0:max_len] #remove unused positions
    #sum all occurences to calculate percentages
    sums = (base_content['G'] + base_content['C'] + base_content['T'] +
            base_content['A'] + base_content['N'])
    ncontent_total = base_content['N'].sum()/sums.sum()
    for key in base_content:
        base_content[key] = base_content[key]/sums
    
    return means, base_content, ncontent_total

def import_position_data_from_source(rootdir, outdir, export=True):
    if not os.path.exists(outdir) and export:
        os.system("mkdir " + outdir)
    
    filenames = []
    all_phred_means = []
    all_base_content = []
    all_n_content = []
    for root, dirs, files in os.walk(rootdir):
        for name in files:
            filepath = root + os.sep + name
            if filepath.endswith(".fastq.gz"):
                phred_means, base_content, n_content = extract_fastq_positions(filepath)
                all_phred_means.append(phred_means)
                all_base_content.append(base_content)
                all_n_content.append(n_content)
                filenames.append(name.replace(".fastq.gz", ""))
                print("finished " + filepath)
    
    fastq_positions = pd.DataFrame(all_base_content, columns=['G', 'C', 'A', 'T', 'N'])
    fastq_positions['filename']= filenames
    fastq_positions['phred_means'] = all_phred_means
    fastq_positions['n_content'] = all_n_content
    if export:
        fastq_positions.to_json(outdir + '/position_data.json')
    return fastq_positions

def import_position_data_from_export(outdir):
    position_data = pd.read_json(outdir+'/position_data.json')
    for category in position_data:
        if type(position_data[category][0]) is list:
            position_data[category] = position_data[category].apply(np.array)
    return position_data

def import_position_data(rootdir, outdir, force_reimport=False, export=True):
    if os.path.exists(outdir+'/position_data.json') and not force_reimport:
        return import_position_data_from_export(outdir)
    else:
        return import_position_data_from_source(rootdir, outdir, export)

def create_single_fastqc_dataframe(fastqc_import, module_list):
    module_result = []
    module_status = []
    
    for module in fastqc_import.modules:
        result_data = pd.DataFrame(fastqc_import.modules[module]['data'])
        result_data.columns = fastqc_import.modules[module]['fieldnames']
        module_status.append(fastqc_import.modules[module]['status'])
        module_result.append(result_data)
    # imported reads don't include the module 'Overrepresented sequences' 
    # if there are none in the read, so we manually add the status list
    # and an empty dataframe /is None better?
    if len(module_status) < 11:
        module_status.insert(9, 'pass')
        module_result.insert(9, pd.DataFrame())
    module_result.append(module_status)
    module_result.append(module_result[0].Value[0].replace(".fastq.gz", ""))
    module_series = pd.Series(data=module_result, index=module_list)
    single_read = module_series.to_frame().T
    return single_read

def import_fastqc_results_from_source(rootdir, outdir, export=True):
    reads = []
    for root, dirs, files in os.walk(rootdir):
        for name in files:
            filepath = root + os.sep + name
            if filepath.endswith("fastqc_data.txt"):
                reads.append(FastQCParser(filepath))

    module_list = list(reads[0].modules.keys())
    # imported reads don't include the module 'Overrepresented sequences' 
    # if there are none in the read, so we manually add the module to the list
    if(len(module_list) < 11):
        module_list.insert(9, 'Overrepresented sequences')
    module_list.append('Module Statuses')
    module_list.append('filename')

    read_list = []
    
    for read in reads:
        single_read = create_single_fastqc_dataframe(read, module_list)
        read_list.append(single_read)
    
    read_results = pd.concat(read_list)
    read_results.reset_index(inplace=True, drop=True)
    if export:
        read_results.to_json(outdir + '/fastqc_results.json')
    return read_results

def import_fastqc_results_from_export(outdir):
    fastq_results = pd.read_json(outdir+'/fastqc_results.json')
    for category in fastq_results:
        if type(fastq_results[category][0]) is dict:
            fastq_results[category] = fastq_results[category].apply(pd.DataFrame)
    return fastq_results

def import_fastqc_results(rootdir, outdir, force_reimport=False, export=True):
    if os.path.exists(outdir+'/fastqc_results.json') and not force_reimport:
        return import_fastqc_results_from_export(outdir)
    else:
        return import_fastqc_results_from_source(rootdir, outdir, export)


def import_all_reads(inputdir, exportdir, force_reimport=False, export=True, 
                     include_metadata=True, include_positiondata=True):
    run_fastqc(inputdir)
    fastqc_results = import_fastqc_results(inputdir, exportdir, force_reimport, export)
    fastqc_results.set_index('filename', inplace=True)
    # merge dataframes
    if include_metadata:
        metadata = extract_metadata(inputdir, exportdir, force_reimport, export)
        metadata.set_index('filename', inplace=True)
        fastqc_results = fastqc_results.join(metadata)
    if include_positiondata:
        position_data = import_position_data(inputdir, exportdir, force_reimport, export)
        position_data.set_index('filename', inplace=True)
        fastqc_results = fastqc_results.join(position_data)
    fastqc_results.to_json(exportdir+'/all_reads.json')
    return fastqc_results