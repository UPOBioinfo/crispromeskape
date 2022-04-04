
"""
This script is used to derive feature information of multiple randomized random forests, 
to infer the importance of genes as predictors for different crispr subtypes in multiple hospital bug species.
Each species and each crispr sub type respectively is compared to the strains of the same species without crispr information
@Author: Maximilian Sprang, Github: Muedi
"""
# %%
import numpy as np
import pandas as pd
from os import listdir
from os.path import isfile, join
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from numpy import random as rd
from sklearn.metrics import confusion_matrix, roc_auc_score, roc_curve, accuracy_score

# %%
# # data files
filenames = ["data paper/" + f for f in listdir("data paper/")
             if isfile(join("data paper/", f))]

# %%
iters = 20
species = ["_ec",  "_ef", "_ab", "_kp", "_pa", "_sa"]
# loop over species
for s in species:
    # load data and metadata
    file_data = [file_data for file_data in filenames if s in file_data and "gnmatrix" in file_data][0]
    data = pd.read_csv(file_data, compression='infer', sep="\t")
    data = data.set_index("#genes/strains")
    data = data.T  # long format
    file_meta = [
        file_meta for file_meta in filenames if s in file_meta and "metadata" in file_meta][0]
    meta = pd.read_csv(file_meta, compression='infer', sep="\t")
    meta = meta.set_index("ID")

    # drop ambigous
    meta = meta.loc[meta.Ambiguous.ne("ambiguous").values]
    # remove MLST group -
    meta = meta.loc[meta.MLST.ne("-").values]
    # select datapoints that remain
    data = data.loc[meta.index]
    # select all crispr subtypes to iterate over below
    crispr_sub = list(set(meta["CRISPRtype"].dropna()))
    crispr_sub = [subtype for types in crispr_sub for subtype in types.split(",")]
    # join data to use for indexing
    data_joined = data.join(meta)
    data_joined["CRISPRtype"] = data_joined["CRISPRtype"].fillna("No_crispr") 
    # loop over subtypes
    for sub in crispr_sub:

        targets = [sub, "No_crispr"]
        data_matrix = data[data_joined.CRISPRtype.str.contains(sub) | data_joined.CRISPRtype.eq("No_crispr")]
        data_matrix_input = data_joined[data_joined.CRISPRtype.str.contains(sub) | data_joined.CRISPRtype.eq("No_crispr")]
        data_matrix_input["Target"] = ["No_crispr" if x == "No_crispr" else sub for x in data_matrix_input.CRISPRtype]

        # run for X iterations
        all_dem_features = pd.DataFrame(
            [0] * len(data_matrix.columns), index=data_matrix.columns)
        all_dem_features_dropped = pd.DataFrame(
            [0] * len(data_matrix.columns), index=data_matrix.columns)
        for r in range(iters):
            random_state = rd.randint(0, 1000)
            # start with base rf, insert grid/random search at some point.
            tmp_inp = data_matrix
            rf = RandomForestClassifier(random_state=random_state)
            X_train, X_test, Y_train, Y_test = train_test_split(
                tmp_inp,
                data_matrix_input["Target"],
                test_size=0.2,
                random_state=random_state
            )
                      # fit for all genes
            rf.fit(X_train, Y_train)

            # save feature importance
            feature_importances = pd.DataFrame(rf.feature_importances_, index=X_train.columns, columns=[
                'importance']).sort_values('importance', ascending=True)
            # feature_importances.sort_values("importance", ascending=False).iloc[:30].to_csv("acinetobacter/results/feature_importances_top30.allGenes.tsv", sep="\t")
            
            # add selected features to the all_dem_features dataframe
            all_dem_features.loc[feature_importances.sort_values(
                "importance", ascending=False).iloc[:30].index] += 1
            y_predicted_train = rf.predict(X_train)
            y_predicted_test = rf.predict(X_test)
            # metrics
            # confusion_matrix(Y_train, y_predicted_train)
            # confusion_matrix(Y_test, y_predicted_test)
            # accuracy_score(Y_test, y_predicted_test)

        # add foldchange and difference between respective crispr subtype groups and no_crispr
        all_dem_features.columns = ["known_genes_retained"]
        feature_importance_output = all_dem_features
        feature_importance_output["log2fc_gene(mean_{}/mean_{})".format(targets[0], targets[1])] = np.log2((tmp_inp.loc[data_matrix_input["Target"].eq(targets[0])].mean(0) + 1) /
                                                                                                        (tmp_inp.loc[data_matrix_input["Target"].ne(targets[0])].mean(0) + 1))
        feature_importance_output["diff_gene(mean_{}-mean_{})".format(targets[0], targets[1])] = (tmp_inp.loc[data_matrix_input["Target"].eq(targets[0])].mean(0)) - \
            (tmp_inp.loc[data_matrix_input["Target"].ne(targets[0])].mean(0))

        # save feature overview
        feature_importance_output.sort_values(["known_genes_retained"], ascending=False).to_csv(
            "output/{}_selected_features_top30_{}_of_{}_iters.csv".format(s, sub, iters))
        
# %%
