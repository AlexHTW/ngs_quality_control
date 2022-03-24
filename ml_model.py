import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import cross_validate
import joblib


def split_by_organism(fastqc_dataset):
    organisms = fastqc_dataset['organism'].unique()
    grouped = fastqc_dataset.groupby(fastqc_dataset.organism)
    data_per_organism = {org: grouped.get_group(org) for org in organisms}
    return data_per_organism

def extract_target(fastqc_dataset):
    target = fastqc_dataset.evaluation.values
    drop_list = ['organism', 'technology', 'read_number', 'evaluation']
    return target, fastqc_dataset.drop(columns=drop_list, errors='ignore')

def train_model(fastqc_dataset, exportdir, evaluate=False, name='complete_data'):
    target, fastqc_dataset = extract_target(fastqc_dataset)
    if evaluate:
        evaluate_model(fastqc_dataset, target, exportdir, name)
    clf_rf = RandomForestClassifier(n_estimators=100)
    clf_rf.fit(fastqc_dataset, target)
    joblib.dump(clf_rf, exportdir+'/model_rf_'+name+'.pkl')

def train_model_per_organism(fastqc_dataset, exportdir, evaluate=False):
    data_per_organism = split_by_organism(fastqc_dataset)
    for org, data in data_per_organism.items():
        train_model(data, exportdir, evaluate, org)

def evaluate_model(fastqc_dataset, target, exportdir, name):
    clf_rf = RandomForestClassifier(n_estimators=100)
    scores = cross_validate(clf_rf, fastqc_dataset, target, cv=10,
                            scoring=('accuracy', 'f1'), return_train_score=True)
    with open(exportdir+'/model_evaluation_'+name+'.txt', 'w') as f:
        f.write("%s Dataset Cross Validation Model Results: \n\nAccuracies: %s \nF1 Scores: %s \nAccuracy: %0.2f +/- %0.2f \nF1 Score: %0.2f +/- %0.2f"
                % (name, scores['test_accuracy'], scores['test_f1'], scores['test_accuracy'].mean(),
                scores['test_accuracy'].std(), scores['test_f1'].mean(), scores['test_f1'].std()))

def predict_evaluation(fastqc_dataset, modelpath):
    if isinstance(fastqc_dataset, pd.Series):
        fastqc_dataset = fastqc_dataset.to_frame().transpose()
    clf = joblib.load(modelpath)
    prediction = clf.predict(fastqc_dataset)
    return prediction