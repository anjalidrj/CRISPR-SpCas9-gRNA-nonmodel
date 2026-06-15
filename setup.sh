#!/bin/bash
echo "Setting up CRISPR Guide Design environment..."

# Install mamba
conda install -y -n base -c conda-forge mamba

# Create or update environment
if conda env list | grep -q "guide-design"; then
    echo "Environment exists, updating..."
    mamba env update -f environment.yml --prune
else
    echo "Creating environment..."
    mamba env create -f environment.yml
fi

# Activate and install additional packages
source /opt/conda/etc/profile.d/conda.sh
conda activate guide-design
mamba install -n guide-design -y matplotlib seaborn
python -m ipykernel install --user --name guide-design --display-name "guide-design"

echo ""
echo "✅ Setup complete!"
echo "Run: source /opt/conda/etc/profile.d/conda.sh && conda activate guide-design"