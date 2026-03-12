conda environment activation

conda create -n candidates python=3.10 -y
conda activate candidates
conda install -c bioconda -c conda-forge --file requirements.txt -y
pip install ./SCALOP/

Use case with example V3_mAb_sequences.xlsx file
All input files must be formatted like V3_mAb_sequences.xlsx, where Column 1 is the Sequence number/ID, Column 2 is the heavy chain sequence, and Column 3 is the light chain sequence

python main.py \
  --xlsx V3_mAb_sequences.xlsx \
  --email your_email0@email_provider.com \
  --out results/
