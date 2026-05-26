export MODULEPATH="/software/modules:$MODULEPATH"
#module load cellgen/nichecompass/0.3.0 
source /software/cellgen/team298/ls34/nemo/bin/activate

echo "LOAD"
python3 /nfs/users/nfs_l/ls34/spatial_skin_atlas/reference_mapping_tutorials/niche_level_resources//nemo.py
echo "DONE"
