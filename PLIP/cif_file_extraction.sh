cd /NFS/hcao/af3output/wildtype_output
source ~/pyenvs/py_modules/bin/activate

python copy_ref_top_af3_cifs.py \
  --excel ref_af3drugs_all10prot_ranked.xlsx \
  --base . \
  --top 100
