#!/bin/bash

##############################################
# Applying trained models to Cell Painting Data from
# The Drug Repurposing Hub
#
# Gregory Way, 2020
# 
# CRITICAL FIX (2024): This pipeline now uses pre-transformed data
# to resolve feature mismatch between Image_* (training) and 
# Cells_* (repurposing) features. Run transform_features.py first
# to generate data/repurposing_transformed_features_fixed.tsv.gz
##############################################

# Check if transformed data exists for testing
if [ ! -f "data/repurposing_transformed_features_fixed.tsv.gz" ]; then
    echo "Transformed data not found. Generating it for testing..."
    echo "Running transform_features.py to create test data..."
    
    if python transform_features.py; then
        echo "Test data generated successfully!"
    else
        echo "ERROR: Failed to generate test data!"
        echo "Please check transform_features.py for errors."
        exit 1
    fi
    echo ""
fi

echo "Test data ready: data/repurposing_transformed_features_fixed.tsv.gz"
echo "Proceeding with pipeline execution..."
echo ""

# Step 0: Convert all notebooks to scripts
jupyter nbconvert --to=script \
        --FilesWriter.build_directory=scripts/nbconverted \
        *.ipynb

# Step 1: Apply the models to the Drug Repurposing Hub
jupyter nbconvert --to=html \
        --FilesWriter.build_directory=scripts/html \
        --ExecutePreprocessor.kernel_name=python3 \
        --ExecutePreprocessor.timeout=10000000 \
        --execute 0.apply-models-repurposing.ipynb

# Step 2: Visualize results
jupyter nbconvert --to=html \
        --FilesWriter.build_directory=scripts/html \
        --ExecutePreprocessor.kernel_name=ir \
        --ExecutePreprocessor.timeout=10000000 \
        --execute 1.visualize-repurposing.ipynb

# Step 3: Fit dose curves
Rscript --vanilla 2.fit-dose.R

# Step 4: Visualize dose response curves
jupyter nbconvert --to=html \
        --FilesWriter.build_directory=scripts/html \
        --ExecutePreprocessor.kernel_name=ir \
        --ExecutePreprocessor.timeout=10000000 \
        --execute 3.visualize-dose-response.ipynb

# Step 5: Generate summary figures
jupyter nbconvert --to=html \
        --FilesWriter.build_directory=scripts/html \
        --ExecutePreprocessor.kernel_name=ir \
        --ExecutePreprocessor.timeout=10000000 \
        --execute 4.lincs-figures.ipynb
