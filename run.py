# Linear
import os
from subprocess import Popen, PIPE
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import random
import math
from sklearn.neighbors import NearestNeighbors
from scipy.stats.stats import pearsonr
from sklearn.linear_model import SGDClassifier
from sklearn import preprocessing
import Queue as Q
import os

n = 0
threshold = 21234/10

def compare(dir):
    paths=[]
    names=[]

    for filename in os.listdir(dir):
        paths.append(dir+filename)
        names.append(filename)
    
    names.sort()
    paths.sort()
    path_num = len(paths)
    exp_num = path_num/2

    #print(names)
    print("Totally {} experiments".format(exp_num))

    for exp_i in range(0,path_num,2): # path_num  path_num
        print paths[exp_i], paths[exp_i+1]
        matrix_file = "/workspace/project1/LINCS/networks/scaled_output"
        proc = Popen(['/workspace/code/3.9/histogram', matrix_file, paths[exp_i], paths[exp_i+1],str(100)], stdout=PIPE)
    	proc.wait(),
	
       
### Main ###
general = '/workspace/code/Dec11/exp/'
compare(general)

        