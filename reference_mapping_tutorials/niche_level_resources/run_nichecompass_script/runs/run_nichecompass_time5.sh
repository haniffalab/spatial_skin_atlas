export MODULEPATH="/software/modules:$MODULEPATH"
module load cellgen/nichecompass/0.3.0                
echo "LOAD"
python3 /nfs/users/nfs_l/ls34/spatial_skin_atlas/reference_mapping_tutorials/niche_level_resources/run_nichecompass_script/runs/nichecompass_refonly_time5.py
echo "DONE"
