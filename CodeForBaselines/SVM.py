import warnings
import sys
warnings.filterwarnings("ignore")
import pandas as pd
from sklearn.svm import SVC
from sklearn import metrics
from sklearn.model_selection import train_test_split
from sklearn.feature_extraction.text import CountVectorizer,TfidfVectorizer  
from sklearn.preprocessing import OneHotEncoder,LabelEncoder
from scipy import sparse
import os
import scipy.io as sio



def test2train(save_name):
    y_pred = sio.loadmat(save_name)['output_mode_pred']
    return y_pred
        # sio.savemat(save_name, {'input_h': X_test/10000000,'output_mode':Y_test,'output_mode_pred': y_pred})
        # # Compute the error for iteration t.
           
def cal_acc(a,b):
    n=a.shape[0]
    m=a.shape[1]
    tterr=0
    r_err=0
    for i in range(n):
        cuerr=0
        for j in range(m):
            if a[i][j]!= b[i][j]:
               tterr+=1
               cuerr+=1
        if cuerr>0:
            r_err+=1
            
    return 1-r_err/n, 1-tterr/(n*m)

KK = 1                     # number of users
K= 1


# Load data
measurement = sio.loadmat('./Data/GeneratedDataset/data14_%d' %K)['X_t']
attack = sio.loadmat('.Data/GeneratedDataset/data14_%d' %K)['label_t']
measurement1 = sio.loadmat('./Data/GeneratedDataset/data14_%d' %K)['X_s']
attack1 = sio.loadmat('./Data/GeneratedDataset/data14_1%d' %K)['label_s']

# online
def SVM_predict(train_x,train_y,test_x):
    print("SVM test")
     clf_SVM = SVC(kernel='rbf', C=.1,gamma='auto')
     clf_SVM.fit(X=train_x, y=train_y)
     pre_SVM = clf_SVM.predict(train_x)
    return clg_SVM.predict_proba(test_x)
pred_y=SVM_predict(measurement,attack,measurement1)
print(cal_acc(pred_y,attack1))
