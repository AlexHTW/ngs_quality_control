from operator import mod
import numpy as np
import pandas as pd
from scipy.stats import beta
from scipy.optimize import curve_fit
from scipy.special import gamma
import os

# Module 1 - Per base sequence quality
def poly_fit(x, a0, a1, a2, a3):
    y = a0 + a1*x + a2*x**2 + a3*x**3
    return y

def phred_means_fitting(row):
    y = row['phred_means'][5:] #remove first 5 values, because they disturb the fitting
    x = np.arange(0, len(y), 1)
    a, _ = curve_fit(poly_fit, x, y, maxfev=10000)
    return a # 4 values

def module_1_features(fastqc_dataset):
    fastqc_dataset[['module_1_a0', 'module_1_a1', 'module_1_a2', 'module_1_a3'
            ]] = fastqc_dataset.apply(phred_means_fitting, axis=1, result_type='expand')
    
    fastqc_dataset.drop(columns='phred_means', inplace=True)
    fastqc_dataset.drop(columns='Per base sequence quality', inplace=True)
    return fastqc_dataset

# Module 2 - Per tile sequence quality

def add_weights_tiles(row):
    bases = [int(d) for d in str(row['Base']).split('-')]
    if len(bases) == 1:
        return 1
    elif len(bases) == 2:
        return bases[1]-bases[0]+1

def means_per_tile(row):
    module_2 = row['Per tile sequence quality'].copy()
    module_2["weight"] = module_2.apply(add_weights_tiles, axis=1)
    weighted_avg = module_2.groupby('Tile').apply(lambda x: np.average(x['Mean'], weights=x['weight']))
    
    std_neg = np.std(weighted_avg[weighted_avg < 0])
    std_pos = np.std(weighted_avg[weighted_avg >= 0])
    #mean_neg = np.mean((weighted_avg[weighted_avg < 0]), axis=0)
    #mean_pos = np.mean((weighted_avg[weighted_avg >= 0]), axis=0)
    #sum_negatives = np.sum(weighted_avg[weighted_avg <0])
    #sum_positives = np.sum(weighted_avg[weighted_avg >= 0])
    #a = (sum_positives - np.abs(sum_negatives))/sum_positives
    #mean_negatives = np.mean((weighted_avg[weighted_avg < 0]), axis=0)
    #mean_positives = np.mean((weighted_avg[weighted_avg >= 0]), axis=0)
    #mean = np.mean(weighted_avg)
    return std_neg, std_pos

def module_2_features(fastqc_dataset):
    fastqc_dataset[['module_2_std_neg', 'module_2_std_pos'
            ]] = fastqc_dataset.apply(means_per_tile, axis=1, result_type='expand')
    fastqc_dataset.drop(columns='Per tile sequence quality', inplace=True)
    return fastqc_dataset

# Module 3 - Per sequence quality scores

def beta_fit(x, a, b):
        y = (gamma(a+b)*x**(a-1)*(1-x)**(b-1))/(gamma(a)*gamma(b))
        return y

def per_sequence_quality_scores(row):
    table = row['Per sequence quality scores'].copy()
    x = table['Quality']
    x = x/(x.iloc[-1]+1)
    y = table['Count']

    a_b, _ = curve_fit(beta_fit, x, y, maxfev=10000)

    return a_b #2 values

def module_3_features(fastqc_dataset):
    fastqc_dataset[['module_3_alpha', 'module_3_beta'
            ]] = fastqc_dataset.apply(per_sequence_quality_scores, axis=1, result_type='expand')
    fastqc_dataset.drop(columns='Per sequence quality scores', inplace=True)
    return fastqc_dataset

# Module 4 - Per base sequence content

def base_content_differences(row):
    g_diff = np.absolute(np.ediff1d(row['G']))
    c_diff = np.absolute(np.ediff1d(row['C']))
    t_diff = np.absolute(np.ediff1d(row['T']))
    a_diff = np.absolute(np.ediff1d(row['A']))
    bases_diff = g_diff + c_diff + t_diff + a_diff
    bases_diff_mean = bases_diff.mean()
    bases_diff_std = bases_diff.std()

    return bases_diff_mean, bases_diff_std

def module_4_features(fastqc_dataset):
    fastqc_dataset[['module_4_diff_mean', 'module_4_diff_std'
            ]] = fastqc_dataset.apply(base_content_differences, axis=1, result_type='expand')
    fastqc_dataset.drop(columns=['G', 'C', 'A', 'T'], inplace=True)
    fastqc_dataset.drop(columns='Per base sequence content', inplace=True)
    return fastqc_dataset

# Module 5 - Per sequence GC content (not covered in feature engineering)
def module_5_features(fastqc_dataset):
    return fastqc_dataset.drop(columns='Per sequence GC content')

# Module 6 - Per base N content

def module_6_features(fastqc_dataset):
    fastqc_dataset.rename(columns={"n_content": "module_6_n_content"}, inplace=True)
    fastqc_dataset.drop(columns='N', inplace=True)
    fastqc_dataset.drop(columns='Per base N content', inplace=True)
    return fastqc_dataset

# Module 7 - Sequence length distribution

def sequence_length_distribution(row):
    table = row['Sequence Length Distribution'].copy()
    #impute values for tables with only one value (changes original table in dataset)
    if len(table) < 2:
        table.loc[1] = [table['Length'][0]-1, 0.0]
        table.loc[2] = table.iloc[0]
        table.iloc[0] = [39, 0.0]
        table.reset_index(drop=True, inplace=True) # not necessary, just to be safe

    x = table['Length'].astype(str).str.split("-").str[-1].astype(float).astype(int)
    x = x/(x.iloc[-1]+1) #+1 because highest x_value has to be lower than 1 
    y = table['Count']
    a_b, _ = curve_fit(beta_fit, x, y, maxfev=10000) 
    
    return a_b # 2 values

def module_7_features(fastqc_dataset):
    fastqc_dataset[['module_7_alpha', 'module_7_beta'
                ]] = fastqc_dataset.apply(sequence_length_distribution, axis=1, result_type='expand')
    fastqc_dataset.drop(columns='Sequence Length Distribution', inplace=True)
    return fastqc_dataset

# Module 8 - Sequence duplication levels

def sequence_duplication_levels(row):
    table = row['Sequence Duplication Levels'][0:9].copy() #only first 10 rows because original scale not linear
    x = table.index.to_series().astype(int)
    x = (x+1)/(x.iloc[-1]+2) #x+1 because beta_fit doesn't allow for x=0 values
    y = table['Percentage of total']/101 #percentage value as fraction, 101 to avoid 1.0 values which lead to NaN

    a_b, _ = curve_fit(beta_fit, x, y, maxfev=10000)

    return a_b # 2 values

def module_8_features(fastqc_dataset):
    fastqc_dataset[['module_8_alpha', 'module_8_beta'
            ]] = fastqc_dataset.apply(sequence_duplication_levels, axis=1, result_type='expand')
    fastqc_dataset.drop(columns='Sequence Duplication Levels', inplace=True)
    return fastqc_dataset

# Module 9 - Overrepresented sequences (not covered in feature engineering)

def module_9_features(fastqc_dataset):
    return fastqc_dataset.drop(columns='Overrepresented sequences')

# Module 10 - Adapter content (not covered in feature engineering)

def module_10_features(fastqc_dataset):
    return fastqc_dataset.drop(columns='Adapter Content')

# run all feature engineering functions
def apply_feature_engineering(fastqc_dataset, exportdir, force_reimport=False, export=True):
    if os.path.exists(exportdir+'/finished_features_all.json') and not force_reimport:
        fastqc_dataset = pd.read_json(exportdir+'/finished_features_all.json')
    else:
        fastqc_dataset = module_1_features(fastqc_dataset)
        fastqc_dataset = module_2_features(fastqc_dataset)
        fastqc_dataset = module_3_features(fastqc_dataset)
        fastqc_dataset = module_4_features(fastqc_dataset)
        fastqc_dataset = module_5_features(fastqc_dataset)
        fastqc_dataset = module_6_features(fastqc_dataset)
        fastqc_dataset = module_7_features(fastqc_dataset)
        fastqc_dataset = module_8_features(fastqc_dataset)
        fastqc_dataset = module_9_features(fastqc_dataset)
        fastqc_dataset = module_10_features(fastqc_dataset)
        if export:
            fastqc_dataset.to_json(exportdir+'/finished_features_all.json')
    return fastqc_dataset