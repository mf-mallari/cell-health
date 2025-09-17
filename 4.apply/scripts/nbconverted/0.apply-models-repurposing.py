#!/usr/bin/env python
# coding: utf-8

# # Apply Cell Health Models to Repurposing Set
# 
# **Gregory Way, 2019**
# 
# The models are trained to predict cell health phenotypes.
# Here, I apply the models to Cell Painting data from the repurposing set.
# 
# I will use these predictions to identify compound perturbation signatures of cell health impact.

# In[1]:


import os
import sys
import numpy as np
import pandas as pd
import scipy.stats
from joblib import load
import umap

from pycytominer import feature_select
from pycytominer.cyto_utils import infer_cp_features

sys.path.append("../3.train")
from scripts.ml_utils import load_train_test, load_models


# In[2]:


np.random.seed(123)


# ## 1) Load Models and Training Data

# In[3]:


consensus = "modz"
output_dir = "data"


# In[4]:


model_dir = os.path.join("..", "3.train", "models")

model_dict, model_coef = load_models(
    model_dir=model_dir,
    consensus=consensus
)


# In[5]:


data_dir = os.path.join("..", "3.train", "data")

x_train_df, x_test_df, y_train_df, y_test_df = load_train_test(
    data_dir=data_dir,
    consensus=consensus,
    drop_metadata=True
)


# ## 2) Load Cell Painting Repurposing Data Files
# 
# These files are available from https://github.com/broadinstitute/lincs-cell-painting

# In[6]:


batch = "2016_04_01_a549_48hr_batch1"
commit_hash = "27a2d7dd74067b5754c2c045e9b1a9cfb0581ae4"

# We have noticed particular technical issues with this platemap
# remove it from downstream consideration
# https://github.com/broadinstitute/lincs-cell-painting/issues/43
filter_platemap = "C-7161-01-LM6-011"


# In[7]:


# Load transformed data for testing (fixes feature mismatch)transformed_file = "data/repurposing_transformed_features_fixed.tsv.gz"if os.path.exists(transformed_file):    print(f"📝 Using transformed data for testing: {transformed_file}")    complete_consensus_df = pd.read_csv(transformed_file, sep='\t')    print(f"  - Shape: {complete_consensus_df.shape}")    print(f"  - Features: {len([c for c in complete_consensus_df.columns if not c.startswith('Metadata_')])}")else:    print("❌ ERROR: Transformed data not found for testing!")    print("Please run: python transform_features.py")    raise FileNotFoundError("Test data not available")


# In[ ]:


# Apply feature selection to the consensus profiles (compatible with transformed data)feature_ops = [    "variance_threshold",    "drop_na_columns",    "blocklist",    "drop_outliers"]consensus_feature_select_df = feature_select(    complete_consensus_df,    operation=feature_ops,    na_cutoff=0)print(consensus_feature_select_df.shape)consensus_feature_select_df.head()# Split metadata and CP Features  # The training data has Image_ features, while LINCS data has Cells_, Nuclei_, Cytoplasm_ features# Get metadata featuresmeta_features = [c for c in complete_consensus_df.columns if c.startswith("Metadata_")]# Get training data features (non-metadata columns)training_features = [c for c in x_test_df.columns if not c.startswith("Metadata_") and c != "level_0"]# Get consensus features after feature selectionconsensus_features = [c for c in consensus_feature_select_df.columns if not c.startswith("Metadata_")]print(f"Training features: {len(training_features)}")print(f"Consensus features after selection: {len(consensus_features)}")# CRITICAL FIX: Find overlapping features between training and consensus dataoverlapping_features = [c for c in training_features if c in consensus_feature_select_df.columns]print(f"Overlapping features: {len(overlapping_features)}")# Use overlapping features if available, otherwise use all consensus features with warningif len(overlapping_features) > 10:  # Reasonable threshold    print(f"✅ Using {len(overlapping_features)} overlapping features")    feature_df = consensus_feature_select_df.loc[:, overlapping_features]else:    print(f"⚠️  WARNING: Only {len(overlapping_features)} overlapping features found.")    print("Using all consensus features - models may not work properly!")    feature_df = consensus_feature_select_df.loc[:, consensus_features]metadata_df = complete_consensus_df.loc[:, meta_features]print(f"Final feature matrix shape: {feature_df.shape}")print("Feature preview:")feature_df.head()


# In[ ]:


# Split metadata and CP Features  
# The training data has Image_ features, while LINCS data has Cells_, Nuclei_, Cytoplasm_ features

# Get metadata features
meta_features = [c for c in complete_consensus_df.columns if c.startswith('Metadata_')]

# Get training data features (non-metadata columns)
training_features = [c for c in x_test_df.columns if not c.startswith('Metadata_') and c != 'level_0']

# Get consensus features after feature selection
consensus_features = [c for c in consensus_feature_select_df.columns if not c.startswith('Metadata_')]

print(f"Training features: {len(training_features)}")
print(f"Consensus features after selection: {len(consensus_features)}")

# CRITICAL FIX: Find overlapping features between training and consensus data
overlapping_features = [c for c in training_features if c in consensus_feature_select_df.columns]
print(f"Overlapping features: {len(overlapping_features)}")

# Use overlapping features if available, otherwise use all consensus features with warning
if len(overlapping_features) > 10:  # Reasonable threshold
    print(f"✅ Using {len(overlapping_features)} overlapping features")
    feature_df = consensus_feature_select_df.loc[:, overlapping_features]
else:
    print(f"⚠️  WARNING: Only {len(overlapping_features)} overlapping features found.")
    print("Using all consensus features - models may not work properly!")
    feature_df = consensus_feature_select_df.loc[:, consensus_features]

metadata_df = complete_consensus_df.loc[:, meta_features]

print(f"Final feature matrix shape: {feature_df.shape}")
print("Feature preview:")
feature_df.head()


# ## 3) Apply all Regression Models to all Repurposing Plates

# In[ ]:


cell_health_features = list(model_dict.keys())

all_scores = {}
for cell_health_feature in cell_health_features:
    model_clf = model_dict[cell_health_feature]
    pred_df = model_clf.predict(feature_df)
    all_scores[cell_health_feature] = pred_df


# ## 4) Output Results

# In[ ]:


# Output scores
all_score_df = pd.DataFrame.from_dict(all_scores)
repurp_predict_df = (
    metadata_df
    .merge(
        all_score_df,
        left_index=True,
        right_index=True
    )
    .query("Metadata_Plate_Map_Name != @filter_platemap")
)

output_real_file = os.path.join(
    output_dir,
    "repurposing_transformed_real_models_{}.tsv.gz".format(consensus)
)
repurp_predict_df.to_csv(output_real_file, sep="\t", index=False, compression="gzip")

print(repurp_predict_df.shape)
repurp_predict_df.head()


# ## 5) Apply UMAP
# 
# ### Part 1: Apply UMAP to Cell Health Transformed Repurposing Hub Features

# In[ ]:


reducer = umap.UMAP(random_state=1234, n_components=2)

predict_embedding_df = pd.DataFrame(
    reducer.fit_transform(repurp_predict_df.loc[:, cell_health_features]),
    columns=["umap_x", "umap_y"]
)

predict_embedding_df = (
    metadata_df
    .merge(
        predict_embedding_df,
        left_index=True,
        right_index=True
    )
    .query("Metadata_Plate_Map_Name != @filter_platemap")
)

output_real_file = os.path.join(
    output_dir,
    "repurposing_umap_transformed_real_models_{}.tsv.gz".format(consensus)
)

predict_embedding_df.to_csv(output_real_file, sep="\t", index=False, compression="gzip")

print(predict_embedding_df.shape)
predict_embedding_df.head()


# ### Part 2: Apply UMAP to All Repurposing Hub Cell Painting Profiles

# In[ ]:


reducer = umap.UMAP(random_state=1234, n_components=2)repurp_embedding_df = pd.DataFrame(    reducer.fit_transform(        consensus_feature_select_df.loc[:, [c for c in consensus_feature_select_df.columns if not c.startswith("Metadata_") and c != "level_0"]]    ),    columns=["umap_x", "umap_y"])repurp_embedding_df = (    metadata_df    .merge(        repurp_embedding_df,        left_index=True,        right_index=True    )    .query("Metadata_Plate_Map_Name != @filter_platemap"))output_real_file = os.path.join(    output_dir,    "repurposing_umap_transformed_cell_painting_{}.tsv.gz".format(consensus))repurp_embedding_df.to_csv(output_real_file, sep="\t", index=False, compression="gzip")print(repurp_embedding_df.shape)repurp_embedding_df.head()


# ## Merge Data Together for Shiny App Exploration

# In[ ]:


# Load MOA file
moa_url = "https://raw.githubusercontent.com/broadinstitute/lincs-cell-painting/"
moa_url = f"{moa_url}/{commit_hash}/metadata/moa/repurposing_info_external_moa_map_resolved.tsv"

moa_df = pd.read_csv(moa_url, sep="\t")

print(moa_df.shape)
moa_df.head(3)


# In[ ]:


core_id = [
    "{}-{}".format(
        x.split("-")[0],
        x.split("-")[1]
    ) if x != "DMSO"
    else x
    for x in repurp_embedding_df.Metadata_broad_sample
]

repurp_embedding_with_pert_df = (
    repurp_embedding_df
    .assign(Metadata_broad_core_id=core_id)
    .sort_index(axis="columns")
    .merge(
        moa_df,
        left_on="Metadata_broad_core_id",
        right_on="broad_id",
        how="left"
    )
)

print(repurp_embedding_with_pert_df.shape)
repurp_embedding_with_pert_df.head()


# In[ ]:


shiny_merge_cols = [
    "Metadata_Plate_Map_Name",
    "Metadata_broad_sample",
    "Metadata_dose_recode",
    "Metadata_mmoles_per_liter",
    "Metadata_pert_well"
]

shiny_df = (
    repurp_embedding_with_pert_df.merge(
        repurp_predict_df,
        left_on=shiny_merge_cols,
        right_on=shiny_merge_cols,
        how="inner"
    )
    .drop(["broad_sample"], axis="columns")
    .query("Metadata_Plate_Map_Name != @filter_platemap")
)

print(shiny_df.shape)
shiny_df.head()


# In[ ]:


shiny_file = os.path.join(
    "repurposing_cellhealth_shiny",
    "data",
    "moa_cell_health_{}.tsv.gz".format(consensus)
)

shiny_df.to_csv(shiny_file, sep='\t', index=False, compression="gzip")


# In[ ]:


# Merge with transformed data compatibility# Use metadata columns for merging since feature structure has changedmetadata_cols = [c for c in complete_consensus_df.columns if c.startswith('Metadata_')]shiny_combined_df = shiny_df.merge(    complete_consensus_df[metadata_cols + ['level_0']],  # Include level_0 for compatibility    on=metadata_cols,    how="inner")


# ## Output Correlation Matrix

# In[ ]:


shiny_features = [c for c in consensus_feature_select_df.columns if not c.startswith("Metadata_") and c != "level_0"]cell_health_features = [x for x in shiny_df if x.startswith("cell_health")]


# In[ ]:


all_results = []
for cell_health_feature in cell_health_features:
    cell_health = shiny_combined_df.loc[:, cell_health_feature]
    for cp_feature in shiny_features:
        feature = shiny_combined_df.loc[:, cp_feature]
        cor_result, pval = scipy.stats.pearsonr(cell_health, feature)
        all_results.append([cell_health_feature, cp_feature, cor_result, pval])


# In[ ]:


# Output correlation matrix for cell health predictions and CellProfiler features
cor_results_df = (
    pd.DataFrame(
        np.array(all_results), columns=["cell_health", "cp_feature", "pearson_cor", "pval"]
    )
    .sort_values(by="pearson_cor", ascending=False)
    .reset_index(drop=True)
)

cor_results_df.pearson_cor = cor_results_df.pearson_cor.astype(float)

cor_results_df = (
    cor_results_df
    .pivot_table(columns=["cell_health"], index=["cp_feature"], values="pearson_cor")
)

print(cor_results_df.shape)
cor_results_df.head(3)

