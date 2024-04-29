# Bayesian approach for the Indian Pimas diabetes data set


## Data sets:
### diabetes.csv: original data set obtained from [Kaggle Pima Indian data set](https://www.kaggle.com/datasets/uciml/pima-indians-diabetes-database/data)
### diabetes_imputed.csv: data set after imputation using posterior predictive distributions
### diabetes_imputed_CE.csv: data set after imputation using chained equations

## Python notebooks
### bayesian imputation: performed Bayessian imputation using both posterior predictive distributions and chained equations. Imputed data saved as csv
### EDA: plots and tables for Exploratory Data Analysis
### bayesian logistic regression: fit Bayesian logistic LASSO for variable selection. Then use selected features for Bayesian logistic regression using diabetes_imputed_CE.csv
### sensitivity analysis: bayesian LASSO and logistic regression using the diabetes_imputed.csv 

