import time
import numpy as np
import pandas as pd
import datetime as dt
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn import metrics
from sklearn.model_selection import KFold
from scipy.interpolate import Rbf
from scipy import stats
from neupy import layers, algorithms
from neupy import plots
from neupy import algorithms
from neupy.layers import *
import dill

import tensorflow as tf
from neupy.utils import tensorflow_session

# Import data
data = pd.read_excel('RealAG.xlsx', index_col=0, header=0)
data.columns = ['Temp(alb)','Temp(In)','Tair','Tdew', 'Wind','HR']

start_time = time.time()
data['Month'] = data.index.month
data['MonthCos'] = np.cos(2 * np.pi * data.index.month / 12)
data['MonthSin'] = np.sin(2 * np.pi * data.index.month / 12)

data['DayOfYear'] = data.index.dayofyear
data['DayOfYearCos'] = np.cos(2 * np.pi * data.index.dayofyear / 365)
data['DayOfYearSin'] = np.sin(2 * np.pi * data.index.dayofyear / 365)

data['WeekOfYear'] = data.index.weekofyear
data['WeekOfYearCos'] = np.cos(2 * np.pi * data.index.weekofyear / 52)
data['WeekOfYearSin'] = np.sin(2 * np.pi * data.index.weekofyear / 52)

data['WeekOfHalfYearCos'] = np.cos(2 * np.pi * data.index.weekofyear / 26)
data['WeekOfHalfYearSin'] = np.sin(2 * np.pi * data.index.weekofyear / 26)

# Add filtered data
tmp2 = data.loc[:,['Temp(In)','Tair','Tdew', 'Wind','HR']].rolling(31, center=False, axis=0, min_periods=1).mean()
tmp2.columns = ['TempinMeanFilter_31', 'TairMeanFilter_31', 'TdewMeanFilter_31','TWindMeanFilter_31','HRMeanFilter_31']

tmp3 = data.loc[:,['Temp(In)','Tair','Tdew', 'Wind','HR']].rolling(31, center=False, axis=0, min_periods=1).std()
tmp3.columns = ['TempinStdFilter_31', 'TairStdFilter_31', 'TdewStdFilter_31','TWindStdFilter_31','HRStdFilter_31']

print(f"Data processing time: {time.time() - start_time:.2f}")

data = pd.concat((data, tmp2, tmp3 ** 2), axis=1)

# Drop empty records
data = data.dropna()

#Define data (Temp)

X = data.iloc[:, 1:]
y = data.loc[:, ['Temp(alb)']]

years = data.index.year
yearsTrain, yearsTest = train_test_split(np.unique(years), test_size=0.3, train_size=0.7, random_state=42)
# yearsTrain, yearsTest = train_test_split(np.unique(years), test_size=0.3, train_size=0.7)

XTrain = X.query('@years in @yearsTrain')
yTrain = y.query('@years in @yearsTrain').values.ravel()
XTest = X.query('@years in @yearsTest')
yTest = y.query('@years in @yearsTest').values.ravel()
results = y.query('@years in @yearsTest')


#===============================================================================
# Neural network
#===============================================================================

# Define neural network

training_times = []
prediction_times = []
all_prediction_times_per_sample = []
n_training_runs = 20

for _ in range(n_training_runs):
    # Note: it's important to clean tensorflow graph before starting a new training
    tf.reset_default_graph()

    if hasattr(tensorflow_session, 'cache'):
        del tensorflow_session.cache

    cgnet = algorithms.Momentum(
        network=[
            layers.Input(XTrain.shape[1]),
            layers.Relu(24),
            layers.Linear(1),
        ],
        step=algorithms.step_decay(
            initial_value=0.05,
            reduction_freq=750,
        ),
        loss='mse',
        batch_size=None,
        regularizer=algorithms.l2(0.002),
        shuffle_data=False,
        verbose=False,
        show_epoch=100,
    )

    XScaler = StandardScaler()
    XScaler.fit(XTrain)
    XTrainScaled = XScaler.transform(XTrain)
    XTestScaled = XScaler.transform(XTest)

    yScaler = StandardScaler()
    yScaler.fit(yTrain.reshape(-1, 1))
    yTrainScaled = yScaler.transform(yTrain.reshape(-1, 1)).ravel()
    yTestScaled = yScaler.transform(yTest.reshape(-1, 1)).ravel()

    # Train
    start_time = time.time()
    cgnet.train(XTrainScaled, yTrainScaled, epochs=5000)

    training_time = time.time() - start_time
    training_times.append(training_time)
    print(f"Training took {round(training_time, 3)} seconds")

    start_time = time.time()
    prediction = cgnet.predict(XTestScaled)
    prediction_time = time.time() - start_time
    prediction_times.append(prediction_time)
    print(f"Prediction took {round(prediction_time, 3)} seconds")

    prediction_times_per_sample = []
    for x_test_sample in np.matrix(XTestScaled):
        start_time = time.time()
        prediction = cgnet.predict(x_test_sample)
        prediction_time_per_sample = time.time() - start_time

        prediction_times_per_sample.append(prediction_time_per_sample)
        all_prediction_times_per_sample.append(prediction_time_per_sample)

    avg_prediction_time_per_sample = np.mean(prediction_times_per_sample)
    print(f"Average prediction time per sample: {round(avg_prediction_time_per_sample, 4)} seconds")

print("")
print("General statistic:")
print("  Number of training samples:", XTrainScaled.shape[0])
print("  Number of test samples:", XTrainScaled.shape[0])
print("  Number of features:", XTrainScaled.shape[1])
print("  Number of parameters:", cgnet.network.n_parameters)
print("  Average training time: {:.2f} +/- {:.2f} (2 standard deviations)".format(np.mean(training_times), 2 * np.std(training_times)))
print("  Average prediction time: {:.5f} +/- {:.5f} (2 standard deviations)".format(np.mean(prediction_times), 2 * np.std(prediction_times)))
print("  Average prediction time per sample: {:.5f} +/- {:.5f} (2 standard deviations)".format(
    np.mean(all_prediction_times_per_sample),
    2 * np.std(all_prediction_times_per_sample),
))

yEstTrain = yScaler.inverse_transform(cgnet.predict(XTrainScaled).reshape(-1, 1)).ravel()
mae = np.mean(np.abs(yTrain-yEstTrain))
start_time = time.time()
prediction = cgnet.predict(XTestScaled)
prediction_time = time.time() - start_time
results['ANN'] = yScaler.inverse_transform(prediction.reshape(-1, 1)).ravel()
print(f"Prediction took {round(prediction_time, 5)} seconds")

# Metrics
mse  = np.mean((yTrain-yEstTrain)**2)
mseTes = np.mean((yTest-results['ANN'])**2)
maeTes = np.mean(np.abs(yTest-results['ANN']))
meantrain = np.mean(yTrain)
ssTest = (yTrain-meantrain)**2
r2=(1-(mse/(np.mean(ssTest))))
meantest = np.mean(yTest)
ssTrain = (yTest-meantest)**2
r2Tes=(1-(mseTes/(np.mean(ssTrain))))


# Plot results
print("NN MAE: %f (All), %f (Test) " % (mae, maeTes))
print ("NN MSE: %f (All), %f (Test) " % (mse, mseTes))
print ("NN R2: %f (All), %f (Test) " % (r2, r2Tes))

results.plot()
plt.show()

plt.scatter(yTest,results['ANN'])
plt.plot([7, 27], [7, 27], color='red')
plt.xlabel('True Values')
plt.ylabel('Predictions')

plt.show(block=True)
