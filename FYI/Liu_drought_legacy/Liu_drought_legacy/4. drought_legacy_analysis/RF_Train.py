#!/usr/bin/env python
# coding: utf-8

import csv
import pandas as pd
import numpy as np
import random as random
import seaborn as sns
import time
from sklearn.ensemble import  RandomForestRegressor
from matplotlib import pyplot
import matplotlib.pyplot as plt
import joblib
import scipy.stats as st
from scipy.stats import norm
from sklearn.metrics import mean_squared_error #MSE
from sklearn.metrics import mean_absolute_error #MAE
from sklearn.metrics import r2_score#R 2

train_data_path = r"D:\Work\RF\RF_drought(OBS-PRE).csv"




def train_test_valid(df, train_ratio = 0.6, valid_ratio = 0.2 ):
    
    train, test, validation = np.split(df.sample(frac = 1, 
                                                 random_state = 42), 
                                       [int(train_ratio * len(df)), 
                                        int((train_ratio + valid_ratio) * len(df))])
    
    return train, test, validation


train_dataset = pd.read_csv(train_data_path)
train_dataset = pd.DataFrame(train_dataset)
train_dataset.replace([np.inf, -np.inf], np.nan)
train_dataset = train_dataset[np.isfinite(train_dataset).all(1)]

train_dataset_drop = train_dataset[~train_dataset.isin([0])].dropna(axis=0)#.drop('AI',axis = 1).drop('end_month',axis = 1).drop('ndvi_min_l',axis = 1)#.drop('sm_min_l',axis = 1).drop('ndvi_sum_l',axis = 1).drop('MAP',axis = 1).drop('MAT',axis = 1)#.drop('ndvi_sum_l',axis = 1)#.drop('LSD_mean',axis = 1).drop('detrend_WT',axis = 1).drop('detrend_WP',axis = 1)



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


r2_score_standard = 0.5
estimuters_irrator =[560,570,580,590,600] 
for n_estimators_temp in estimuters_irrator:

    RFmodel = RandomForestRegressor(n_estimators=n_estimators_temp,verbose = 1,n_jobs=-1)

    RFmodel.fit(x_train,y_train)
    # save

    y_test_pred = RFmodel.predict(x_test)
    df_pred_pred = pd.DataFrame(np.array(x_test))
    df_pred_pred['y_test'] = np.array(y_test)
    df_pred_pred['y_test_pred'] = y_test_pred

    y1 = y_test
    y2 = y_test_pred
    if r2_score(y1,y2)>r2_score_standard:
        r2_score_standard =r2_score(y1,y2)
        print(n_estimators_temp)
        print("mean_squared_error",format(mean_squared_error(y1,y2)))
        print("mean_absolute_error",format(mean_absolute_error(y1,y2)))
        print("RMSE",format(np.sqrt(mean_squared_error(y1,y2))))
        print("r2",r2_score(y1,y2))

        importance = RFmodel.feature_importances_
        joblib.dump(RFmodel,"K:/Liu_drought_legacy/RF/model/RF_train_pdp_type2_l(sos).joblib")
        
        

feature_importance = RFmodel.feature_importances_

feature_name = x_valid.columns.tolist()

feature_importance = pd.DataFrame(feature_importance)


feature_importance['names'] = feature_name


feature_importance = feature_importance.sort_values(by=0,ascending=True)
feature_importance = feature_importance.transpose()

plt_featurename=feature_importance.iloc[1].tolist()
feature_importance_values = feature_importance.iloc[0].tolist()

plt.figure(1,figsize=(8, 13))
plt.barh(plt_featurename, feature_importance_values, height=0.8)

feature_importance.to_csv("K:/Liu_drought_legacy/RF/model/RF_variables_feature_importances.csv",index=True,sep=',')

