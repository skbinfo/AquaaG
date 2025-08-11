#
# MAKER OPTS TEMPLATE
#

#-----Genome (populated by the script)
genome=
organism_type=eukaryotic

#-----Evidence (populated by the script)
#est=
protein=
rmlib=

#-----Repeat Masking
# The following are left blank to force MAKER to use the
# repeat library and masked genome provided by the script.
model_org=
repeat_protein=

#-----Gene Prediction
snaphmm= # Populated by script for Round 2
augustus_species= # Add a species if you have one (e.g., arabidopsis, human)
est2genome=0 
protein2genome=1

#-----Computational Resources
cpus=1 # This will be overwritten by the script

#-----MAKER Behavior Options
max_dna_len=100000
min_contig=1
keep_preds=1
