#!/bin/bash
#
# Based on Gregory Way 2019
# Predicting Cell Health
# Pipeline to Process Data
#
# 1) Aggregate and Normalize Single Cell Data
# 2) Curate and Normalize Cell Health Assay Data

set -e

# Step 0: Convert all notebooks to scripts
jupyter nbconvert --to=script \
        --FilesWriter.build_directory=scripts/nbconverted \
        *.ipynb

# Step 1: Use pycytominer to generate and process all single cell profiles
python generate_profiles.py

# Step 1.5: Clean and standardize profiles for feature consistency
echo "Step 1.5: Running data cleaning to fix feature consistency..."
python data_cleaning_fixes.py

# Step 2: Process cell health assay data (these are the labels for training models)
jupyter nbconvert --to=html \
        --FilesWriter.build_directory=scripts/html \
        --ExecutePreprocessor.kernel_name=python3 \
        --ExecutePreprocessor.timeout=10000000 \
        --execute 0.process-labels.ipynb

# Step 3: Normalize cell health assay data and generate .gct files for heatmaps
jupyter nbconvert --to=html \
        --FilesWriter.build_directory=scripts/html \
        --ExecutePreprocessor.kernel_name=python3 \
        --ExecutePreprocessor.timeout=10000000 \
        --execute 1.normalize-labels.ipynb

# Step 4: Build consensus signatures for downstream processing
jupyter nbconvert --to=html \
        --FilesWriter.build_directory=scripts/html \
        --ExecutePreprocessor.kernel_name=python3 \
        --ExecutePreprocessor.timeout=10000000 \
        --execute 2.build-consensus-signatures-fix.ipynb

# Step 5: Build a supplementary table describing perturbations
jupyter nbconvert --to=html \
        --FilesWriter.build_directory=scripts/html \
        --ExecutePreprocessor.kernel_name=python3 \
        --ExecutePreprocessor.timeout=10000000 \
        --execute 3.perturbation-table.ipynb

# Step 6: Retrieve an example cell painting image
# Note: This step requires the raw image data!
#jupyter nbconvert --to=html \
#        --FilesWriter.build_directory=scripts/html \
#        --ExecutePreprocessor.kernel_name=python3 \
#        --ExecutePreprocessor.timeout=10000000 \
#        --execute 4.load-example-cell-painting-image.ipynb

# Step 7: Produce a summary of cell counts
jupyter nbconvert --to=html \
        --FilesWriter.build_directory=scripts/html \
        --ExecutePreprocessor.kernel_name=python3 \
        --ExecutePreprocessor.timeout=10000000 \
        --execute 5.cell-count-summary.ipynb

# Step 8: Generate cell health feature heatmaps (local visualization)
jupyter nbconvert --to=html \
        --FilesWriter.build_directory=scripts/html \
        --ExecutePreprocessor.kernel_name=python3 \
        --ExecutePreprocessor.timeout=10000000 \
        --execute 6.visualize-cell-health-heatmap.ipynb

echo "Pipeline complete!"
echo ""
echo "Visualization Options:"
echo "   PDFs: figures/cell_health_feature_heatmap*.pdf"
echo "   Shiny app: cd cell_health_heatmap_shiny && ./launch_app.sh"
echo ""
