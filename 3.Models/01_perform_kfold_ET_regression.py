import regress_sklearn
import stats_tools
import pandas
import numpy
import sklearn.model_selection
import os


# Open the CSV file as a Pandas data frame - the df variable.
# Specify the first column (i.e., independent variable names)
# to be the index values.
df = pandas.read_csv('GBR_depth.csv')
df = df.drop(df[df.s2green == 0].index)

# Get a list of the columns within the df dataframe
cols = list(df.columns)

# Get the indepedent predictor column names
ind_vars = cols[7:9]
#7:9 = good

# Get the dependent response column names
dep_var = cols[4]

print(ind_vars)
print(dep_var)


# Get the predictor variables and dependent variables
# from the dataframe as numpy arrays
x = df[ind_vars].values
y = df[dep_var].values


# Randomly sample the input data to create training and testing (20% sample) datasets
# so we have an independent dataset to test the quality of the relationship
x_train, x_test, y_train, y_test = sklearn.model_selection.train_test_split(x, y, test_size=0.2, random_state=0)

y_train = numpy.expand_dims(y_train, axis=1)
y_test = numpy.expand_dims(y_test, axis=1)

from sklearn.ensemble import ExtraTreesRegressor
skregrs_obj = ExtraTreesRegressor(n_estimators=100)
#from sklearn.neighbors import KNeighborsRegressor
#skregrs_obj = KNeighborsRegressor(n_neighbors=4, algorithm='brute')
#from sklearn.kernel_ridge import KernelRidge
#skregrs_obj = KernelRidge()


metrics, residuals = regress_sklearn.perform_kfold_fit(skregrs_obj, x_train, y_train, n_splits=5, repeats=20, shuffle=False, data_scaler=None)


df_metrics = pandas.DataFrame(data=metrics[0])
# Save the dataframe to a CSV file.
out_csv_file = "IS_GBR_Regres_Metrics_ET_depth_comp_20201030.csv"
df_metrics.to_csv(out_csv_file)

df_residuals = pandas.DataFrame(data=residuals[0])
# Save the dataframe to a CSV file.
out_csv_file = "IS_GBR_Regres_Residuals_ET_depth_20201030.csv"
df_residuals.to_csv(out_csv_file)

