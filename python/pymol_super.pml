reinitialize
import sys
import os

python 

target_protein=sys.argv[3]
moving_protein=sys.argv[4]

target_protein_name=os.path.basename(target_protein)[:-4]
moving_protein_name=os.path.basename(moving_protein)[:-4]


print target_protein, target_protein_name
print moving_protein, moving_protein_name
python end

cmd.load(target_protein)
cmd.load(moving_protein)
#show_as cartoon, all

cmd.align(moving_protein_name,target_protein_name)

cmd.save(moving_protein, moving_protein_name)
