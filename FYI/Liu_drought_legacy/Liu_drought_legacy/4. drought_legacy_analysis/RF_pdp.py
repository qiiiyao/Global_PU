#!/usr/bin/env python
# coding: utf-8


import csv
import pandas as pd
from matplotlib import pyplot
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import random as random
import seaborn as sns
import time
import threading
import time
from itertools import product


def train_test_valid(df, train_ratio = 0.6, valid_ratio = 0.2 ):
    
    train, test, validation = np.split(df.sample(frac = 1, 
                                                 random_state = 42), 
                                       [int(train_ratio * len(df)), 
                                        int((train_ratio + valid_ratio) * len(df))])
    
    return train, test, validation
train_data_path_filename = ['RF_drought(OBS-PRE)','RF_drought_type1_l(OBS-PRE)','RF_drought_type2_l(OBS-PRE)',
                            'RF_drought_type3_l(OBS-PRE)']

for i in range(4):
    train_data_path = "D:/RF//"+train_data_path_filename[i]+".csv"
    train_dataset = pd.read_csv(train_data_path)
    train_dataset = pd.DataFrame(train_dataset)
    train_dataset.replace([np.inf, -np.inf], np.nan)
    train_dataset = train_dataset[np.isfinite(train_dataset).all(1)]

    train_dataset_drop = train_dataset[~train_dataset.isin([0])].dropna(axis=0)#.drop('AI',axis = 1).drop('end_month',axis = 1).drop('ndvi_min_l',axis = 1)#.drop('sm_min_l',axis = 1).drop('ndvi_sum_l',axis = 1).drop('MAP',axis = 1).drop('MAT',axis = 1)#.drop('ndvi_sum_l',axis = 1)#.drop('LSD_mean',axis = 1).drop('detrend_WT',axis = 1).drop('detrend_WP',axis = 1)

    import scipy.stats as st

    from scipy.stats import norm

    input_feature_names = ['MAT','MAP','AGB','Biodiversity','Biomes','LGS','Root depth',
                        'Soil clay','Soil sand','GPP anomaly','EOS anomaly','SM anomaly',
                        'Timing of drought','Drought duration','SM loss','SOS_change']   

       for input_features in input_feature_names:
        data = train_dataset_drop[input_features]
        confidence_level = 0.95 
        mean = data.mean()
        std_dev = data.std()
        lower, upper = norm.interval(confidence_level, loc=mean, scale=std_dev)
        train_dataset_drop = train_dataset_drop[ (train_dataset_drop[input_features] >= lower) & (train_dataset_drop[input_features] <= upper) ]

    train_data, test_data, validation_data = train_test_valid(pd.DataFrame(train_dataset_drop))

    # seperate input data from target variable
    x_train, y_train = train_data.iloc[:,:-1], train_data.iloc[:,-1]
    x_test, y_test = test_data.iloc[:,:-1], test_data.iloc[:,-1]
    x_valid, y_valid = validation_data.iloc[:,:-1], validation_data.iloc[:,-1]


    RFmodelpath = "K:/Liu_drought_legacy/RF/model/"+train_data_path_filename[i]+".joblib"
    RFmodel = joblib.load(filename=RFmodelpath)

    feature_importance = RFmodel.feature_importances_

    feature_name = x_valid.columns.tolist()

    feature_importance = pd.DataFrame(feature_importance)


    feature_importance['names'] = feature_name


    feature_importance = feature_importance.sort_values(by=0,ascending=True)
    feature_importance = feature_importance.transpose()

    ############
    target = 'SOS_change'
    for column_interest in feature_name:
        print(column_interest)
        column_of_interest = column_interest
        data = validation_data
        model = RFmodel

        confidence_level = 0.95
        
        data = data.drop(target, axis= 1)
        x_s = data[column_of_interest]
        x_s.dropna(axis = 0)

        tempfile = x_s.round(3).unique()
        print(tempfile.size)
        if (tempfile.size > 100):

            a_min = x_s.min()
            a_max = x_s.max()   
            interval_num = 100
            x_range = np.linspace(a_min,a_max,num=interval_num,endpoint=True)
            x_range_mean = np.zeros(interval_num)
            x_range_upper = np.zeros(interval_num)
            x_range_lower = np.zeros(interval_num)
        else:
            x_range = tempfile
            x_range_mean = np.zeros(tempfile.size)
            x_range_upper = np.zeros(tempfile.size)
            x_range_lower = np.zeros(tempfile.size)
        
        aaa = 0

        for ii in x_range:

            data[column_of_interest]=ii
        
            train_pred = model.predict(data)
            mean = train_pred.mean()
            std_dev = train_pred.std()

            x_range_mean[aaa] = mean
            aaa = aaa+1
            print(aaa)
        

        x_range = pd.DataFrame(x_range)

        x_range_mean = pd.DataFrame(x_range_mean)
        x_range_mean.columns = [column_interest]

        output_mean = pd.concat([x_range,x_range_mean],axis=1)
        path = "K:/Liu_drought_legacy/RF/model/PHP/"+train_data_path_filename[i]+'_'+column_interest+'_feature_mean.csv'
        output_mean.to_csv(path)

    
    print('finish')
    

