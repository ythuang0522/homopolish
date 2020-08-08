import os
import joblib
import numpy as np
import pandas as pd
from os import listdir
from os.path import isfile, join
from modules import preprocessing as p
from modules.FileManager import FileManager
from sklearn.svm import SVC
from multiprocessing import Pool
from sklearn.model_selection import train_test_split
from sklearn.metrics import confusion_matrix,classification_report, accuracy_score

def get_file_paths_from_directory(directory_path):
    """
    Returns all paths of files given a directory path
    :param directory_path: Path to the directory
    :return: A list of paths of files
    """
    file_paths = [os.path.abspath(join(directory_path, file)) for file in listdir(directory_path)
                  if isfile(join(directory_path, file)) and file[-7:] == 'feather']
    return file_paths

def read_feather(filename):
    'converts a filename to a pandas dataframe'
    return pd.read_feather(filename)

def train_model(train_dir, output_dir, output_prefix, threads):

    output_dir = FileManager.handle_output_directory(output_dir)
    train_file = get_file_paths_from_directory(train_dir)

    with Pool(processes=threads) as pool: 
        df_list = pool.map(read_feather, train_file)
        train = pd.concat(df_list, ignore_index=True)
    train = p.haplotype(train)

    train = p.haplotype(train)
    temp = train[train['label'] != 5]
    
    #drop duplicate for label "no correction"
    dup = train[train['label'] == 5].drop_duplicates(subset=train.columns.difference(['position']))
    uni = train[train['label'] == 5].drop_duplicates(subset=train.columns.difference(['position']), keep=False)
    train = pd.concat((temp, dup, uni)).reset_index(drop=True)
    Y_train = train['label'].values

    print('Label frequency: \n', train.label.value_counts())
    print('Origin Training shape: ', train.shape)

    X_train = p.preprocessing(train)
    print(X_train)
    X_train = pd.DataFrame(X_train.drop(['label'], axis=1))
    X_train, X_val, Y_train, Y_val = train_test_split(X_train, Y_train, test_size = 0.1, random_state=1)
    print('After Train/Test split shape: ', X_train.shape)

    clf = SVC(random_state=0, probability=True, class_weight='balanced')
    clf.fit(X_train, Y_train)
    output = output_dir+'/'+output_prefix + '.pkl'
    joblib.dump(clf, output) 
    Y_pred = clf.predict(X_val)

    print(confusion_matrix(Y_val, Y_pred))
    print('Accuracy score: ', accuracy_score(Y_val, Y_pred))
    print(classification_report(Y_val, Y_pred))