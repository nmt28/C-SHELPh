import pandas
import joblib
from sklearn.ensemble import ExtraTreesRegressor


# Open the CSV file as a Pandas data frame - the df variable.
# Specify the first column (i.e., independent variable names)
# to be the index values.
df = pandas.read_csv('GBR_depth.csv')
df = df.drop(df[df.s2green == 0].index)

# Get a list of the columns within the df dataframe
cols = list(df.columns)

# Get the indepedent predictor column names
ind_vars = cols[7:]
print(ind_vars)

# Get the dependent response column names
dep_var = cols[4]

# Get the predictor variables and dependent variables
# from the dataframe as numpy arrays
x = df[ind_vars].values
y = df[dep_var].values

# Fit regression model:
regrs_mdl = ExtraTreesRegressor(n_estimators=100)
regrs_mdl.fit(x, y)

# Output model file
output_mdl_file = 'IS_Ber_ET_Model.joblib'

# Save regression model to disk:
print('Saving regression model...')
joblib.dump(regrs_mdl, output_mdl_file, ('gzip', 1))
print('Saved regression model.')

# Verify that the model can be read from disk:
print('Attempting to read regression model from disk...')
regrs_test_mdl = joblib.load(output_mdl_file)
regrs_test_mdl.predict(x)
print('Read model and used it predict successfully.')

