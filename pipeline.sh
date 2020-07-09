#!/usr/bin/env bash

#Splitting trajectories into single frames
minnie splitpdbs -i sox4 sox18 --pdbs sox4.pdb sox18.pdb

# This part calculates all types of inter-monomer interactions
minnie findbonds -i sox4 -d 'sox4/01_frames/' --itypes all --clean
minnie findbonds -i sox18 -d 'sox18/01_frames/' --itypes all --clean

# Apply critical interaction filter
minnie timefilter -f sox4/02_interfacea_results/*/sox4_merged_*.csv -i sox4 --per 25
minnie timefilter -f sox18/02_interfacea_results/*/sox18_merged_*.csv -i sox18 --per 25

# Calculate common and distinct interactions in two cases
minnie comparecx -i sox4 sox18 --per 25

# aaaaand graphs!
minnie graph -i 'sox4' 'sox18' --per 25 --itypes all -b protein-dna -s specific
minnie graph -i 'sox4' 'sox18' --per 25 --itypes all -b protein-dna -s common


#To visualize minnie results with Pymol!
for projectID in  sox4 sox18
do
  cd "$projectID/04_compare_complex/25_freq_filtered/25_freq_perres/complex_specific"
  sed 's/,/ /g' "$projectID_hbonds_compared_25_perres.csv" | grep "protein-dna" | grep ":" | awk '{printf "show sticks, (resi %s and chain %s) + (resi %s and chain %s) \n", $6,$2,$7,$3}' | sort -u > pymol_hbonds.pml
  sed 's/,/ /g' "$projectID_hbonds_compared_25_perres.csv" | grep "protein-dna" | grep ":" | awk '{printf "distance i. %s and n. %s and chain %s, i. %s and n. %s and chain %s\n", $6,$8,$2,$7,$9,$3}' | sort -u | sed 's/OP1/O1P/g' | sed 's/OP2/O2P/g' >> pymol_hbonds.pml


  sed 's/,/ /g' "$projectID_hydrophobic_compared_25_perres.csv" | grep "protein-dna" | grep ":" | awk '{printf "show sticks, (resi %s and chain %s) + (resi %s and chain %s) \n", $6,$2,$7,$3}' | sort -u > pymol_hydrophobic.pml
  sed 's/,/ /g' "$projectID_hydrophobic_compared_25_perres.csv" | grep "protein-dna" | grep ":" | awk '{printf "distance i. %s and n. %s and chain %s, i. %s and n. %s and chain %s\n", $6,$8,$2,$7,$9,$3}' | sort -u | sed 's/OP1/O1P/g' | sed 's/OP2/O2P/g' >> pymol_hydrophobic.pml


  sed 's/,/ /g' "$projectID_ionic_compared_25_perres.csv" | grep "protein-dna" | grep ":" | awk '{printf "show sticks, (resi %s and chain %s) + (resi %s and chain %s) \n", $6,$2,$7,$3}' | sort -u > pymol_ionic.pml
  sed 's/,/ /g' "$projectID_ionic_compared_25_perres.csv" | grep "protein-dna" | grep ":" | awk '{printf "distance i. %s and n. %s and chain %s, i. %s and n. %s and chain %s\n", $6,$8,$2,$7,$9,$3}' | sort -u | sed 's/OP1/O1P/g' | sed 's/OP2/O2P/g' >> pymol_ionic.pml

  sed 's/,/ /g' "$projectID_ring_stacking_compared_25_perres.csv" | grep "protein-dna" | grep ":" | awk '{printf "show sticks, (resi %s and chain %s) + (resi %s and chain %s) \n", $6,$2,$7,$3}' | sort -u > pymol_ring_stacking.pml
  sed 's/,/ /g' "$projectID_ring_stacking_compared_25_perres.csv" | grep "protein-dna" | grep ":" | awk '{printf "distance i. %s and n. %s and chain %s, i. %s and n. %s and chain %s\n", $6,$8,$2,$7,$9,$3}' | sort -u | sed 's/OP1/O1P/g' | sed 's/OP2/O2P/g' >> pymol_ring_stacking.pml
  cd -
done
