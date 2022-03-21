import pandas as pd
import numpy as np

def create_status_column_name(module_nr):
    return 'module_' + str(module_nr) + '_status'

def prepare_fastqc_data(fastqc_dataset):
    total_sequences = []
    percent_gc = []
    min_sequence_length = []
    max_sequence_length = []

    for i in range(fastqc_dataset.shape[0]):
        total_sequences.append(fastqc_dataset['Basic Statistics'][i]['Value'][3])
        percent_gc.append(fastqc_dataset['Basic Statistics'][i]['Value'][6])
        length_min_max = str(fastqc_dataset['Basic Statistics'][i]['Value'][5]).split('-')
        min_sequence_length.append(length_min_max[0])
        max_sequence_length.append(length_min_max[-1])
    total_sequences = np.asarray(total_sequences, dtype=np.int64)
    percent_gc = np.asarray(percent_gc, dtype=np.int64)
    min_sequence_length = np.asarray(min_sequence_length, dtype=np.int64)
    max_sequence_length = np.asarray(max_sequence_length, dtype=np.int64)
    fastqc_dataset['total_sequences'] = total_sequences
    fastqc_dataset['percent_gc'] = percent_gc
    fastqc_dataset['min_sequence_length'] = min_sequence_length
    fastqc_dataset['max_sequence_length'] = max_sequence_length
    fastqc_dataset.drop(columns='Basic Statistics', inplace=True)

    # convert evaluation to numeric type
    fastqc_dataset['evaluation'].replace({'ugly': 0, 'good': 1}, inplace=True)

    # convert statuses to numeric type
    status_replacements = {'fail':0, 'warn':1, 'pass':2}
    for i in range(fastqc_dataset.shape[0]):
        statuses_list = [status_replacements.get(n, n) for n in fastqc_dataset['Module Statuses'].iat[i]]
        fastqc_dataset['Module Statuses'].iat[i] = statuses_list
    
    # add inner array data from Module Statuses as df columns
    column_names_status = list(map(create_status_column_name, range(len(fastqc_dataset['Module Statuses'][0]))))

    statuses_df = pd.DataFrame(fastqc_dataset['Module Statuses'].to_list(), columns=column_names_status)
    fastqc_dataset = fastqc_dataset.join(statuses_df)
    fastqc_dataset.drop(columns='Module Statuses', inplace=True)

    # remove status of module 0 (Basic Statistics is always pass/2) 
    fastqc_dataset.drop('module_0_status', axis=1, inplace=True)

    return fastqc_dataset