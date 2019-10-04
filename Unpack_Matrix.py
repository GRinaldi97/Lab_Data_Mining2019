import os
import pandas as pd


data_frame= pd.DataFrame()
path=os.getcwd()
path_tr = path+'/TRANSCRIPTOME_DATA'
list_directories = [path_tr+'/'+ele for ele in os.listdir(path_tr) if os.path.isdir(path_tr+'/'+ele)]
for dir in list_directories:
        path_sample_type = dir
        expr_directories = [path_sample_type+'/'+items for items in os.listdir(path_sample_type) if os.path.isdir(path_sample_type+'/'+items)]
        for expr_dir in expr_directories:
            List = [expr_dir+'/'+items2 for items2 in os.listdir(expr_dir) if os.path.isdir(expr_dir+'/'+items2)]
        for directory in List:
            H = [directory+'/'+roba for roba in os.listdir(directory) if roba.endswith('.gz')==False]
            for elem in H:
                if data_frame.empty == True:
                    data_frame=pd.read_csv(elem, sep='\t')
                else:
                    place_holder= pd.read_csv(elem, sep='\t')
                    data_frame = pd.merge(data_frame, place_holder, on='Gene')
data_frame.to_csv('datas_2.txt', sep='\t', header=True, index=None)


