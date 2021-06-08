#!/usr/bin/env python3

import sys
from time import time
import logging
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from sklearn.model_selection import train_test_split
#from sklearn.model_selection import RepeatedKFold
from sklearn.model_selection import LeaveOneOut
from sklearn.model_selection import cross_validate
from sklearn.model_selection import cross_val_predict
from sklearn.naive_bayes import BernoulliNB
from sklearn import metrics


parser = argparse.ArgumentParser()
parser.add_argument("-d", "--data", action="store", required=True)
parser.add_argument("-p", "--processors", action="store", default=8)
args = parser.parse_args()

raw = pd.read_csv(args.data,
        dtype={"Positive": np.int32})

bcols = {
        "CT_HPO_FileID": "pid", 
        "all_HPO_clean": "terms", 
        }
tcols = {
        "ResearchID": "pid", 
        "seq_HPO_clean": "terms"
        }
# split raw into background & target dataframes
bdat = raw.loc[raw["ResearchID"].isnull()][["CT_HPO_FileID", "all_HPO_clean"]].rename(columns=bcols)
tdat = raw.loc[raw["ResearchID"].notnull()][["ResearchID", "seq_HPO_clean", "Positive"]].rename(columns=tcols)

# one-hot-encode HPO terms
b_onehot = bdat.join(bdat["terms"].str.get_dummies(sep="_"))
t_onehot = tdat.join(tdat["terms"].str.get_dummies(sep="_"))

b_onehot["outcome"] = 0
t_onehot["outcome"] = 1

# recombine background & target data
df = b_onehot.merge(t_onehot, how="outer", on=None).fillna(0, downcast="infer")

y = df["outcome"].to_numpy()
terms = df.columns[pd.Series(df.columns).str.startswith("HP:")]
X = df[terms].to_numpy()

# cross-validation
loo = LeaveOneOut()

# Bernoulli naive bayes model
clf = BernoulliNB()
scores = cross_validate(clf, X, y, 
        cv=loo,
        scoring=["accuracy"],#,"roc_auc"],
        return_train_score=True,
        n_jobs=args.processors)

print("Average model accuracy across validation sets")
print(scores["test_accuracy"].mean())

# generate class probabilities predictions
predictions = cross_val_predict(clf, X, y, cv=loo, method="predict_proba", n_jobs=args.processors)

# compute AUC & plot ROC
#fpr, tpr, thresholds = metrics.roc_curve(y, predictions[:,1], pos_label=1)
#roc_auc = metrics.auc(fpr, tpr)
#roc_plot = metrics.RocCurveDisplay(fpr=fpr, tpr=tpr, roc_auc=roc_auc)
#roc_plot.plot()
#plt.show()
#plt.clf()

# join probabilities with positive diagnosis indicator
pred_df = pd.DataFrame(predictions, columns=["neg_proba", "pos_proba"]).join(df[["pid", "Positive", "outcome"]])
# calculate log-probabilities and create rank indicator
pred_df["pos_score"] = -np.log(pred_df["pos_proba"])
pred_df["pos_rank"] = pred_df["pos_score"].rank(method="first").astype("int32")

# plot kernel density of log-probabilities
#dist_plot = sns.kdeplot(pred_df["pos_score"], bw_adjust=0.5)
#plt.savefig("pred_log_proba_kde.png")
#plt.clf()

#print(pred_df)
