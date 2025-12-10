############# Import Packages ##############
import argparse
import os
import shap  # SHAP
import pandas as pd
import numpy as np
import pickle  # save object
import matplotlib.pyplot as plt  # plot

from sklearn.ensemble import RandomForestRegressor
from sklearn.ensemble import AdaBoostRegressor
from sklearn import metrics  # measurements for model performance
from sklearn.model_selection import KFold, GridSearchCV  # tuning parameters
from scipy.stats import nbinom 
from sklearn.model_selection import RandomizedSearchCV

############# command line arguments ##############
parser = argparse.ArgumentParser()
parser.add_argument('--inputtraining', help='input training file with subject_id, outcomes, and predictors',
                    type=argparse.FileType('r'), required=True)
parser.add_argument('--inputvalidation', help='input validation file with subject_id, outcomes, and predictors',
                    type=argparse.FileType('r'), required=True)
parser.add_argument('--Xname', help='input file only including one outcome showing the predictors',
                    type=argparse.FileType('r'), required=True)                    
parser.add_argument('--outcome', help='outcome variable name', required=True)
parser.add_argument('--outputpath', help='destination directory', required=True)

parser.add_argument('--scoring', help='scoring function for RF', default='customize', 
                   choices=['neg_mean_squared_error','neg_mean_squared_log_error','neg_median_absolute_error','explained_variance','customize','r2','customize2'])
args = parser.parse_args()

############# Load data #############
Xname = pd.read_csv(args.Xname)
Xname
## training
df_training = pd.read_csv(args.inputtraining)
# df_training = df_training.set_index(df_training.columns[0]) # make 1st column as rownames
#  The set_index() function in pandas is used to set one or more columns as the index of the DataFrame.
df_training = df_training.set_index('src_subject_id')
xnames = Xname['Xname'].tolist()
## validation
df_validation = pd.read_csv(args.inputvalidation)
df_validation = df_validation.set_index('src_subject_id')

## training
Xt = df_training[xnames].copy()
yt = df_training[args.outcome]#.astype('category')
 
## validation
Xv = df_validation[xnames].copy()
yv = df_validation[args.outcome]#.astype('category')
 
# set output path
os.chdir(args.outputpath)
 
print(Xt.columns)
print(Xt.shape)
print(Xv.shape)

############ define a function to calculate mean shap values #############
def calculate_mean_with_sign_check(x):
    positive_counts = np.sum(np.sign(x) == 1, axis=0)
    negative_counts = np.sum(np.sign(x) == -1, axis=0)
    
    if np.any(positive_counts >= 80) or np.any(negative_counts >= 80):
        return np.mean(x, axis=0)
    else:
        return np.zeros(x.shape[1])
    
    
#############define a function to calculate new scoring function #########
from sklearn.metrics import make_scorer, mean_squared_log_error

def rmsle(y_true, y_pred):
    return np.sqrt(mean_squared_log_error(y_true, y_pred))

# Create a custom scorer using make_scorer
rmsle_scorer = make_scorer(rmsle, greater_is_better=False)  # For minimizing RMSLE, set greater_is_better to False

# calculate negative log-likelihood for a negative binomial distribution 
def neg_binomial_log_likelihood(y_true, y_pred, alpha):
    y_pred = np.maximum(y_pred, 1e-6)
    log_likelihood = -nbinom.logpmf(y_true, alpha, alpha / (alpha + y_pred))
    return np.mean(log_likelihood)

# Create a custom scorer using make_scorer
# alpha is estimated from the data alpha=max(1/(var(y_true)/mean(y_true) -1),1-e6)

def estimate_negative_binomial_alpha(y_true):
    # Calculate sample mean and variance
    sample_mean = np.mean(y_true)
    sample_variance = np.var(y_true)
    
    # Calculate sample dispersion
    sample_dispersion = sample_variance / sample_mean - 1
    
    # Handle edge case when mean is close to zero
    if sample_mean < 1e-6:
        return np.inf  # Return infinity for extremely small means
    
    # Calculate alpha, ensuring it's at least 1
    alpha = max(1 / sample_dispersion, 1)
    
    return alpha

# For the negative binomial distribution, the mean (mu) and variance (sigma^2) are related to alpha through mu = alpha * p / (1 - p) and
# sigma^2 = mu + mu^2 / alpha, where p is the probability of success in a Bernoulli trial.
# https://library.virginia.edu/data/articles/getting-started-with-negative-binomial-regression-modeling
def estimate_alpha_method_of_moments(y):
    mean_y = np.mean(y)
    var_y = np.var(y)
    # Ensure denominator is at least a small positive value
    denominator = max(var_y - mean_y, 1e-6)
    
    alpha = mean_y**2 / denominator
    return alpha

# use empirical data to estimate the alpha value - dispersion value
# alpha =  estimate_negative_binomial_alpha(yt)
alpha =  estimate_alpha_method_of_moments(yt)
neg_binomial_scorer = make_scorer(neg_binomial_log_likelihood, alpha=alpha)  # Adjust alpha as needed

# d2 pinball score for quantile regression, if alpha=0.5, then same as d2 absolute error
from sklearn.metrics import d2_pinball_score
d2_pinball_score_075 = make_scorer(d2_pinball_score, alpha=0.75)

if args.scoring == 'customize':
    scoring = neg_binomial_scorer
elif args.scoring == 'customize2':
    scoring = rmsle_scorer
else:
    scoring = args.scoring

############## modeling ########################
### Random Forest ###
# setup tuning parameters
#"scale_pos_weight": [1, 3, 5],
# Define the parameter grid for RandomForestRegressor
# min_samples_split:
# It determines the minimum number of samples required to split an internal node.
# Higher values prevent the tree from splitting nodes that have fewer samples, potentially reducing model complexity.
# min_samples_leaf:
# It sets the minimum number of samples required to be at a leaf node.
# Using this parameter can smooth the model by preventing nodes with very few samples.
# max_samples: This parameter controls the size of the bootstrap samples used when building trees. 
# It can be useful to tune this parameter, especially when dealing with large datasets.
# random_state: Setting a random_state ensures reproducibility of results. While not strictly a hyperparameter, 
# it's good practice to include it for consistency in model evaluation.
# oob_score: Setting oob_score to True allows you to use out-of-bag samples to estimate the R^2 on unseen data. This can be a useful metric for evaluating model performance.

#Increase Regularization: Regularization techniques like increasing the minimum number of samples required to split an internal node (min_samples_split) and the minimum number of samples required to be at a leaf node (min_samples_leaf) can help prevent overfitting. Experiment with higher values for these parameters.

param_grid_rf = {
    "n_estimators": [500,1000],
    "max_depth": [3, 4, 5, 6],
    "min_samples_split": [6, 8, 10],
    "min_samples_leaf": [4, 6, 8, 10],
    "max_features": ['sqrt', 'log2'],
    "criterion": ['squared_error','friedman_mse'],  # Setting criterion as 'squared_error' 'friedman_mse', 'absolute_error', 'poisson'
    "bootstrap": [True], #add bootstrap parameter
    "max_samples":[0.5, 0.6, 0.7, 0.8, 0.9]
    #"min_impurity_decrease": [0.1, 0.15, 0.2]
    #"oob_score": [True] # to use out-of-bag samples to estimate the R^2 on unseen data
    # Add any other Random Forest specific parameters here
}

# Setup a RandomForestRegressor
rf_model = RandomForestRegressor()

# Initialize variables and lists for storing results
# (retain the other variables you've previously used for storing XGBoost results)
# setup iteration number
n = 10

meanabs_table_training = pd.DataFrame({'Predictors': Xt.columns})
meanabs_table_validation = pd.DataFrame({'Predictors': Xv.columns})
yt_pred_table = pd.DataFrame({'yt_true': yt})
yv_pred_table = pd.DataFrame({'yv_true': yv})
perf_training_rmse = []
perf_validation_rmse = []
perf_training_mse = []
perf_validation_mse = []
perf_training_mae = []
perf_validation_mae = []
perf_training_r2 = []
perf_validation_r2 = []
shap_mat_training_list = list()
shap_mat_validation_list = list()

var_imp_training = pd.DataFrame({'Predictors': Xt.columns})

# Experiment with RandomizedSearchCV for hyperparameter tuning
from sklearn.model_selection import RandomizedSearchCV

for i in range(0, n):
    print(i)
    
    # Initialize KFold for cross-validation
    kf = KFold(n_splits=5, shuffle=True, random_state=i) 

    rand_cv = RandomizedSearchCV(estimator=rf_model, param_distributions=param_grid_rf, n_iter=50,
                            n_jobs=-1, cv=kf, scoring=scoring, verbose=1, refit=True)
    
    # Fit the RandomForestRegressor with the randomized search
    _ = rand_cv.fit(Xt, yt)

    # Get the best parameters and best score
    best_params_rand = rand_cv.best_params_
    print(best_params_rand)
    best_score_rand = rand_cv.best_score_

    # Use the best model with the best parameters from RandomizedSearchCV
    best_model = rand_cv.best_estimator_
    
    # Evaluate the best model on the validation set
    validation_score = best_model.score(Xv, yv)
    
    from sklearn.metrics import mean_squared_error, mean_absolute_error
    
    yt_pred_table[i] = best_model.predict(Xt)
    yv_pred_table[i] = best_model.predict(Xv)

    print(best_model.feature_importances_)
    
    var_imp_training[i] = best_model.feature_importances_

    # Score
    print("RMSE of rf for training is: %f" % (np.sqrt(metrics.mean_squared_error(yt, best_model.predict(Xt)))))
    print("RMSE of rf for validation is: %f" % (np.sqrt(metrics.mean_squared_error(yv, best_model.predict(Xv)))))
    print("MSE of rf for training is: %f" % (metrics.mean_squared_error(yt, best_model.predict(Xt))))
    print("MSE of rf for validation is: %f" % (metrics.mean_squared_error(yv, best_model.predict(Xv))))
    print("R2 of rf for training is: %f" % (metrics.r2_score(yt, best_model.predict(Xt))))
    print("R2 of rf for validation is: %f" % (metrics.r2_score(yv, best_model.predict(Xv))))

    perf_training_rmse.append(np.sqrt(metrics.mean_squared_error(yt, best_model.predict(Xt))))
    perf_validation_rmse.append(np.sqrt(metrics.mean_squared_error(yv, best_model.predict(Xv))))
    perf_training_mse.append(metrics.mean_squared_error(yt, best_model.predict(Xt)))
    perf_validation_mse.append(metrics.mean_squared_error(yv, best_model.predict(Xv)))
    perf_training_mae.append(metrics.mean_absolute_error(yt, best_model.predict(Xt)))
    perf_validation_mae.append(metrics.mean_absolute_error(yv, best_model.predict(Xv)))
    perf_training_r2.append(metrics.r2_score(yt, best_model.predict(Xt)))
    perf_validation_r2.append(metrics.r2_score(yv, best_model.predict(Xv)))

    # Create a TreeExplainer with your trained model
    explainer_rf = shap.TreeExplainer(best_model, feature_perturbation="tree_path_dependent", model_output="raw")

    shap_values_rf_training = explainer_rf.shap_values(Xt, tree_limit=-1, approximate=False, check_additivity=True)
    shap_values_rf_training_df = pd.DataFrame(shap_values_rf_training, columns=Xt.columns).set_index(df_training.index)
    shap_mat_training_list.append(shap_values_rf_training_df)
    meanabs_training = shap_values_rf_training_df.abs().mean(axis=0).to_list()
    meanabs_table_training[i] = meanabs_training

    shap_values_rf_validation = explainer_rf.shap_values(Xv, tree_limit=-1, approximate=False, check_additivity=True)
    shap_values_rf_validation_df = pd.DataFrame(shap_values_rf_validation, columns=Xv.columns).set_index(df_validation.index)
    shap_mat_validation_list.append(shap_values_rf_validation_df)
    meanabs_validation = shap_values_rf_validation_df.abs().mean(axis=0).to_list()
    meanabs_table_validation[i] = meanabs_validation

print('The mean RMSE Score of ', n, ' rf runs for training is', sum(perf_training_rmse)/len(perf_training_rmse))
print('The mean RMSE Score of ', n, ' rf runs for validation is', sum(perf_validation_rmse)/len(perf_validation_rmse))
print('The mean MAE Score of ', n, ' rf runs for training is', sum(perf_training_mae)/len(perf_training_mae))
print('The mean MAE Score of ', n, ' rf runs for validation is', sum(perf_validation_mae)/len(perf_validation_mae))

perf_training_rmse.append(sum(perf_training_rmse)/len(perf_training_rmse))
perf_training_rmse = pd.DataFrame({'perf_training_rmse': perf_training_rmse})
perf_training_mse.append(sum(perf_training_mse)/len(perf_training_mse))
perf_training_mse = pd.DataFrame({'perf_training_mse': perf_training_mse})
perf_training_mae.append(sum(perf_training_mae)/len(perf_training_mae))
perf_training_mae = pd.DataFrame({'perf_training_mae': perf_training_mae})
perf_training_r2.append(sum(perf_training_r2)/len(perf_training_r2))
perf_training_r2 = pd.DataFrame({'perf_training_r2': perf_training_r2})
pd.DataFrame.to_csv(perf_training_rmse, 'perf_training_rmse_rf.csv', sep=',', index=True, index_label='runs')
pd.DataFrame.to_csv(perf_training_mse, 'perf_training_mse_rf.csv', sep=',', index=True, index_label='runs')
pd.DataFrame.to_csv(perf_training_mae, 'perf_training_mae_rf.csv', sep=',', index=True, index_label='runs')
pd.DataFrame.to_csv(perf_training_r2, 'perf_training_r2_rf.csv', sep=',', index=True, index_label='runs')

perf_validation_rmse.append(sum(perf_validation_rmse)/len(perf_validation_rmse))
perf_validation_rmse = pd.DataFrame({'perf_validation_rmse': perf_validation_rmse})
perf_validation_mse.append(sum(perf_validation_mse)/len(perf_validation_mse))
perf_validation_mse = pd.DataFrame({'perf_validation_mse': perf_validation_mse})
perf_validation_mae.append(sum(perf_validation_mae)/len(perf_validation_mae))
perf_validation_mae = pd.DataFrame({'perf_validation_mae': perf_validation_mae})
perf_validation_r2.append(sum(perf_validation_r2)/len(perf_validation_r2))
perf_validation_r2 = pd.DataFrame({'perf_validation_r2': perf_validation_r2})
pd.DataFrame.to_csv(perf_validation_rmse, 'perf_validation_rmse_rf.csv', sep=',', index=True, index_label='runs')
pd.DataFrame.to_csv(perf_validation_mse, 'perf_validation_mse_rf.csv', sep=',', index=True, index_label='runs')
pd.DataFrame.to_csv(perf_validation_mae, 'perf_validation_mae_rf.csv', sep=',', index=True, index_label='runs')
pd.DataFrame.to_csv(perf_validation_r2, 'perf_validation_r2_rf.csv', sep=',', index=True, index_label='runs')

pd.DataFrame.to_csv(yt_pred_table, 'yt_pred_rf_table.csv', sep=',', index=True, index_label='src_subject_id')
pd.DataFrame.to_csv(yv_pred_table, 'yv_pred_rf_table.csv', sep=',', index=True, index_label='src_subject_id')

########### training #############
pd.DataFrame.to_csv(meanabs_table_training, 'meanabs_table_rf_training.csv', sep=',', index=False)

pd.DataFrame.to_csv(var_imp_training, 'var_imp_training.csv', sep=',', index=False)

temp = pd.concat(shap_mat_training_list, axis=0)
# temp = temp.abs()
data_temp = temp.groupby(temp.index).mean()
# plot
shap_mat_rf_training100mean = data_temp.to_numpy()
f = plt.figure()
shap.summary_plot(shap_mat_rf_training100mean, Xt, max_display=len(Xt.columns), show=False)
f.savefig("shap_rf_training100mean.png", bbox_inches='tight', dpi=600)

# save
data_temp.insert(0, 'outcome', yt)
pd.DataFrame.to_csv(data_temp, 'shap_mat_rf_training.n.mean.csv', sep=',', index=True)

run = [x for x in list(range(1, n+1)) for _ in range(len(yt))]
temp.insert(0, 'run', run)
pd.DataFrame.to_csv(temp, 'shap_mat_rf_training.n.runs.csv', sep=',', index=True)

############# validation #############
pd.DataFrame.to_csv(meanabs_table_validation, 'meanabs_table_rf_validation.csv', sep=',', index=False)

temp = pd.concat(shap_mat_validation_list, axis=0)
# temp = temp.abs()
data_temp = temp.groupby(temp.index).mean()
# plot
shap_mat_rf_validation100mean = data_temp.to_numpy()
f = plt.figure()
# , plot_type="layered_violin"
shap.summary_plot(shap_mat_rf_validation100mean, Xv, max_display=len(Xv.columns), show=False)
f.savefig("shap_rf_validation100mean.png", bbox_inches='tight', dpi=600)

# save
data_temp.insert(0, 'outcome', yv)
pd.DataFrame.to_csv(data_temp, 'shap_mat_rf_validation.n.mean.csv', sep=',', index=True)

run = [x for x in list(range(1, n+1)) for _ in range(len(yv))]
temp.insert(0, 'run', run)
pd.DataFrame.to_csv(temp, 'shap_mat_rf_validation.n.runs.csv', sep=',', index=True)
