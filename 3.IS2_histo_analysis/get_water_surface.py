import pandas as pd
import numpy as np

def bin_data(dataset, lat_res, height_res):
    # Calculate number of bins required
    lat_bin_number = round(abs(dataset['latitude'].min() -  dataset['latitude'].max())/lat_res)
    height_bin_number = round(abs(dataset['photon_height'].min() - dataset['photon_height'].max())/height_res)
    # Add bins to dataset
    dataset1 = dataset
    lat_bins = pd.cut(dataset['latitude'], lat_bin_number, labels = np.array(range(lat_bin_number)))
    dataset1['lat_bins'] = lat_bins
    height_bins = pd.cut(dataset['photon_height'], height_bin_number, labels = np.round(np.linspace(dataset1_test['photon_height'].min(), dataset1_test['photon_height'].max(), num=height_bin_number), decimals = 1))
    dataset1['height_bins'] = height_bins
    dataset1 = dataset1.reset_index(drop=True)
    
    return dataset1.values
    
def get_sea_height(binned_data):
    # Create sea height list
    sea_height = []
    # Group data by latitude
    grouped_data = binned_data.groupby(['lat_bins'], group_keys=True)
    data_groups = dict(list(grouped_data))
    # Loop through groups and return average sea height
    for k,v in data_groups.items():
        tmp_df = pd.DataFrame(v.groupby('height_bins').count())
        largest_h_bin = tmp_df['latitude'].argmax()
        largest_h = tmp_df.index[largest_h_bin]
        lat_bin_sea_median = v.loc[v['height_bins']==largest_h, 'photon_height'].median()
        sea_height.append(lat_bin_sea_median)
        del tmp_df
        
    return np.array(sea_height)


'''

binned_data = bin_data(dataset1_test, 10, 0.5)
sea_height = get_sea_height(binned_data)

'''
