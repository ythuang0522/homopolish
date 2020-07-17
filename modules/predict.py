import feather
import joblib
import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from multiprocessing import Pool
from sklearn.preprocessing import MinMaxScaler
from modules import preprocessing as p 

def predict(dataframe, model, thread, path):    

    df = feather.read_dataframe(dataframe)

    result = df
    df = p.haplotype(df)
    df = p.preprocessing(df)

    size = 1000
    list_of_X = [df.loc[i:i+size-1,:] for i in range(0, len(df), size)]
    
    model = joblib.load(model)
    jobs = (delayed(model.predict_proba)(chunk) for chunk in list_of_X)
    parallel = Parallel(n_jobs=thread)
    result_prob = parallel(jobs)

    prob = []
    for i in result_prob:
        prob.extend(i)

    prob = [tupl for tupl in prob ] #tuple to list
    predict = [ele.argmax() for ele in prob] # get index of max element
    del_prob = [ele[4] for ele in prob]
    # prob = -np.sort([-ele for ele in prob]) # sort element in list
    

    prediction = path + '/result.feather'
    result['del_prob'] = del_prob
    result['predict'] = predict
    result.to_feather(prediction)
    return prediction