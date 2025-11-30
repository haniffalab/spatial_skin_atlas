export MODULEPATH="/software/modules:$MODULEPATH"
module load cellgen/nichecompass/0.3.0                
echo "LOAD"
python3 /nfs/users/nfs_l/ls34/spatial_skin_atlas/reference_mapping_tutorials/nichecompass_ref_mapping_pretrainedmodel.py
echo "DONE"
