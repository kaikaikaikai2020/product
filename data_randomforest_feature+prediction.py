import pandas as pd
import os
import numpy as np

file_name = 'twd_data_processed2.csv'
os.chdir('../data/')
print(os.listdir())
data = pd.read_csv(file_name)
data['random'] = np.random.random(size=len(data))

variable_31 =['px_last_spot', 'usd_swap_ois_6m', 'usd_swap_ois_12m', 'vix_index', 'oil_future', 'em_captial_flow_proxy','tamsci_index', 'spx_index', 'tw_infadj_2y_spd','tw_effec_excerate_jpm', 'tw_cpi_effec_excerate_jpm',\
'tw_ppi_effec_excerate_jpm', 'tw_eco_concensus_change_citi', 'tw_em_eco_surprise_citi', 'tw_inflation_consensus_change_citi','tw_elec_export_yoy', 'dxy_index', 'twd_ind_prod_yoy', 'tw_cpi_yoy','tw_budget_balance', 'tw_10yr_govtnote',\
'ntn31_yield_pd', 'ntn31_return', 'ret_5_ntn31_adj_ccy', 'ret_10_ntn31_adj_ccy', 'ret_20_ntn31_adj_ccy', 'ret_30_ntn31_adj_ccy', 'ret_40_ntn31_adj_ccy', 'ret_5_ntn31_ema', 'ret_10_ntn31_ema', 'ret_20_ntn31_ema'\
'ret_30_ntn31_ema', 'ret_40_ntn31_ema', 'random','target_mix135_ntn31']

variable_61 = ['px_last_spot', 'usd_swap_ois_6m', 'usd_swap_ois_12m', 'vix_index', 'oil_future', 'em_captial_flow_proxy','tamsci_index', 'spx_index', 'tw_infadj_2y_spd','tw_effec_excerate_jpm', 'tw_cpi_effec_excerate_jpm',\
'tw_ppi_effec_excerate_jpm', 'tw_eco_concensus_change_citi', 'tw_em_eco_surprise_citi', 'tw_inflation_consensus_change_citi','tw_elec_export_yoy', 'dxy_index', 'twd_ind_prod_yoy', 'tw_cpi_yoy','tw_budget_balance', 'tw_10yr_govtnote',\
'ntn61_yield_pd', 'ntn61_return', 'ret_5_ntn61_adj_ccy', 'ret_10_ntn61_adj_ccy', 'ret_20_ntn61_adj_ccy', 'ret_30_ntn61_adj_ccy', 'ret_40_ntn61_adj_ccy', 'ret_5_ntn61_ema', 'ret_10_ntn61_ema', 'ret_20_ntn61_ema'\
'ret_30_ntn61_ema', 'ret_40_ntn61_ema', 'random','target_mix135_ntn61']

from sklearn.ensemble import RandomForestRegressor
from matplotlib import pyplot
from sklearn import preprocessing
from sklearn.preprocessing import MinMaxScaler
from sklearn.preprocessing import StandardScaler
values = data[variable_61].dropna().values
values = values.reshape((len(values), len(variable_61)))

scaler = MinMaxScaler(feature_range=(0,1))
scaler = scaler.fit(values)

scaler_std = StandardScaler()
scaler_std = scaler_std.fit(values)

normalize61_std = scaler_std.transform(values)
normalize61_std = scaler.transform(values)
data61 = data[variable_61]
data61 = data61.dropna()
array = normalize61

X61 = array[:, 0:-1]
y61 = array[:,-1]

X61_train, X61_test, y61_train, y61_test = train_test_split(X61, y61, test_size = 0.2, random_state=1, shuffle= False)

model61 = RandomForestRegressor(n_estimators=300, oob_score = True, random_state=1)
model61.fit(X61_train, y61_train)
print(model61.feature_importances_)
names61 = data61.columns.values[0:-1]
ticks61 = [i for i in range(len(names61))]

from sklearn.model_selection import KFold
from sklearn.model_selection import TimeSeriesSplit, cross_val_score
from sklearn.metrics import mean_squared_error
tscv = TimeSeriesSplit(n_splits = 15)
scores = cross_val_score(model61, X61_train, y61_train, cv= tscv, scoring='r2')

print('scores: ', scores)

from sklearn.feature_selection import RFE
from sklearn.ensemble import RandomForestRegressor

rfe61 = RFE(RandomForestRegressor(n_estimators=300, random_state=1), 5)
fit61 = rfe61.fit(X61_train, y61_train)

from sklearn.neural_network import MLPRegressor
nmMod = MLPRegressor(max_iter=250, random_state =1, solver='lbfgs')
nmMod = nmMod.fit(X61_train,y61_train)
y_predict = nmMod.predict(X61_test)

#roll_sharpe = single_return.rolling(180).mean()/single_return.rolling(180).std()*15.5

def cvScore(clf, X,y, sample_weight, scoring = 'neg_log_loss', t1=None, cv=None, cvGen= None, pctEmbargo= None):
	if scoring not in ['neg_log_loss', 'accuracy']:
		raise Exception('wrong scoring method')
	from sklearn.metrics import log_loss, accuracy_score
	from clfSequential import PurgedKFold
	if cvGen is None:
		cvGen = PurgedKFold(n_split=cv, t1=t1, pctEmbargo = pctEmbargo)
	score =[]
	for train, test in cvGen.split(X=X):
		fit = clf.fit(X = X.iloc[train].values)
		if scoring == 'neg_log_loss':
			prob = fit.predict_proba(X.iloc[test:])
			score_ = -log_loss(y.iloc[test], prob, sample_weight=sample_weight.iloc[test].values, labels= clf.classes_)
		else:
			pred = fit.predict(X.iloc[test, :])
			score_ = accuracy_score(y.iloc[test].values, pred, sample_weight= sample_weight.iloc[test].values)
		score.append(score_)
	return np.array(score)

def performTimeSeriesCV(X_train, y_train, number_folds, algorithm, parameters):
	print('parameters: ', parameters)
	print('Size train set: ', X_train.shape)
	k = int(np.floor(X_train.shape[0])/number_folds)
	accuracies = np.zeros(number_folds-1)
	for i in range(2, number_folds-1):
		print('')
		split = float(i-1)/i
		print('splitting the first '+str(i) +'chung at '+str(i-1) +'/'+str(i))
		X = X_train[: (k*i)]
		y = y_train[: (k*i)]
		pritn('size pf train + test: ', X.shape)
		index = int(np.floor(X.shape[0]*split))
		X_trainFolds = X[:index]
		y_trainFolds = y[:index]
		X_testFolds = X[(index+1):]
		y_testFolds = y[(index+1):]
		accuracies[i-2] = performClassification(X_trainFolds, y_trainFolds, X_testFolds, y_testFolds, algorithm, parameters)
		print('accuracy on fold '+str(i)+':', acc[i-2])
	return accuracies.mean()
def performClassification(X_train, y_train, X_test, y_test, method, parameters):
	if method = 'RF':
		return performRFClass(X_train, y_train, X_test, y_test, parameters)
	if method = 'QDA':
		return performQDAClass(X_train, y_train, X_test, y_test, parameters)

def performRFClass(X_train, y_train, X_test, y_test, parameters):
	from sklearn.ensemble RandomForestClassifier
	clf = RandomForestClassifier(n_estimators = 300, n_jobs=-1)
	clf.fit(X_train, y_train)
	accuracy = clf.score(X_test,y_test)
	return accuracy
def performQDAClass(X_train, y_train, X_test, y_test, parameters):
	from sklearn.discriminant_analysis import QuadraticDiscriminantAnalysis
	def replaceTiny(x):
		if(abs(x)<0.0001):
			x=0.0001
	X_train = X_train.apply(replaceTiny)
	X_test = X_test.apply(replaceTiny)
	clf = QDA()
	clf.fit(X_train, y_train)
	accuracy = clf.score(X_test, y_test)
	return accuracy

performTimeSeriesCV(X_train, y_train, 10,'RF',[])

def permutation_importances(rf, X_train, y_trainFolds, metric):
	baseline = metric(rf, X_train, y_train)
	imp=[]
	for col in X_train.columns:
		save = X_train[col].copy()
		X_train[col]=np.random.permutation(X_train[col])
		m = metric(rf, X_train, y_train)
		X_train[col]-save
		imp.apend(baseline -m)
	return np.array(imp)

def dropcol_importances(rf, X_train, y_train):
	rf_ = clone(rf)
	rf_.random_state = 999
	rf_.fit(X_train, y_train)
	baseline = rf_.oob_score_
	imp = []
	for col in X_train.columns:
		X = X_train.drop(col, axis -1 )
		rf_ = clone (rf)
		rf_.random_state =999
		rf_.fit(X, y_train)
		o = rf_.oob_score_
		imp.append(baseline -o)

	imp = np.array(imp)
	I = pd.DataFrame(data={'Feature':X_train.columns,'Importance':imp})
	I = I.set_index('Feature')
	I = I.sort_values('Importance', ascending=True)
	return I 

import os 
file_name = 'twd_ptspd_return.csv'
os.chdir('../data/')
df = pd.read_csv(file_name)
from fbprophet import Prophet
col = ['ds','y']

df_31 = df[col]
m = Prophet()
m.fit(df_31)
future = m.make_future_dataframe (periods=200)
future.tail()
forecast = m.predict(future)


