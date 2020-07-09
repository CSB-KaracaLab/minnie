<img src="logo.png" alt="logo" width="200" />


minnie - a structural ensemble analysis package to deduce the fingerprints of binding

developed by Deniz Dogan and Ezgi Karaca

## Motivation
The molecular dynamics (MD) simulations produce extremely valuable information on the functionally important interactions. A prominent way to define these interactions relies on the characterization of inter-molecular contacts that are persistently similar or different across different simulation trajectories. Expanding on this notion, we have developed the python package, minnie (Molecular INteractioN fIngErprints), which:

(i) calculates  inter-/intra-molecular hydrogen bonds, ionic, hydrophobic and ring stacking interactions by using interfacea python package;
(ii) filters out the non-essential interactions via a user-defined observation frequency filter;
(iii) compares the common and distinct interactions across given simulation sets;
(iv) outputs the common and distinct interaction profiles in table format and compares their distribution with box-and-whisker plots.

## Capabilities
minnie compares the simulation sets of any complex type (i.e. protein-protein/protein-DNA/protein-ligand complex). As long as the residue numbering is comparable, minnie finds similar and distinct interactions between two simulations of the same system, between two highly homologous systems, as well as between wild-type and mutant states of the complex under study.

## Clone the repository
```
git clone https://github.com/CSB-KaracaLab/minnie.git
```
or, if you don't use git:
```
wget https://github.com/CSB-KaracaLab/minnie/archive/master.zip
```

## Set up the right environment for minnie
minnie functions under the umbrella of anaconda3 (with Python 3.6). So, before everything, anaconda3 should be installed.
```
cd minnie
conda env create -f environment.yml
conda activate minnie
cd interfacea
python setup.py build
python setup.py install
cd ..
bash execute.sh

```
setup.sh script installs the required dependencies together with the first version of [interfacea](https://zenodo.org/record/3516439#.Xs91eBMzZBw).
setup.sh finally calls execute.sh script, which sets proper aliases and updates your .bashrc file.

## Activate minnie

When the environment is set, preferentially start a new terminal, alternatively:
```
source ~/.bashrc
```

## **Run minnie**

**1) Your trajectory should be saved in an ensemble format while different conformers are separated by ENDMDL. Then `splitpdbs` option of minnie is used to split your trajectories into single frames.**\
\
The resulting frames are saved under projectID1/02_frames & projectID2/02_frames.
```
minnie splitpdbs -cn projectID1 projectID2 -p ensemble1.pdb ensemble2.pdb
```
In our toy model (as provided under example_run), the project id folders (-cn, i.e. complex names) are defined as sox4 & sox18, where the ensembles are fed as sox4.pdb & sox18.pdb.
```
minnie splitpdbs -cn sox4 sox18 -p sox4.pdb sox18.pdb
```
\
**2) When the individual frames are generated in your project folders, their inter-monomer interactions (hydrogen bonds, ionic, hydrophobic and ring stacking interactions) are calculated with the `findbonds` option. The interactions here are calculated with the interfacea python package.**\
\
The results are saved in the csv format under projectID1/03_interfacea_results & projectID2/03_interfacea_results.\
**`-i all`** calculates hbonds, ionic, hydrophobic, ring_stacking interactions.
You also have the option to calculate intra-monomer interactions with **`-intra True`** flag.
```
minnie findbonds -cn projectID1 -p projectID1/02_frames/* -i all

minnie findbonds -cn projectID2 -p projectID2/02_frames/* -i all
```
In the toy model, the project id folders are kept as defined in the previous step, i.e., sox4 & sox18.
```
minnie findbonds -cn sox4 -p sox4/02_frames/* -i all

minnie findbonds -cn sox18 -p sox18/02_frames/* -i all
```
\
**3) The observed interactions are filtered according to a user-defined observation frequency, such that only the interactions that there at least for the x% of the simulation time will be kept. This is achieved with the `timefilter` option.**
```
minnie timefilter -f projectID1/03_interfacea_results/*/projectID1_merged_*.csv -cn projectID1 --per observation_freq

minnie timefilter -f projectID2/03_interfacea_results/*/projectID2_merged_*.csv -cn projectID2 --per observation_freq
```
The toy model interactions, which are observed for at least 25% of the simulation are kept.
```
minnie timefilter -f sox4/03_interfacea_results/*/sox4_merged_*.csv -cn sox4 --per 25

minnie timefilter -f sox18/03_interfacea_results/*/sox18_merged_*.csv -cn sox18 --per 25
```

**4) At this step of minnie, common and distinct interaction profiles among the two filtered cases are calculated with the `compareCX` option.**
```
minnie compareCX -cn projectID1 projectID2 --per observation_freq
```
The 25%-filtered interactions of the toy systems are analyzed.
```
minnie compareCX -cn sox4 sox18 --per 25
```

**5) The common and distinct interaction distribution comparisons are delivered as box-and-whisker plots.**
```
minnie graph -cn 'projectID1' 'projectID2' --per 25 -i all -b interaction_type -s specific

minnie graph -cn 'projectID1' 'projectID2' --per 25 -i all -b interaction_type -s common
```
which translate into the following for our toy model:
```
minnie graph -cn 'sox4' 'sox18' --per 25 -i all -b protein-dna -s specific

minnie graph -cn 'sox4' 'sox18' --per 25 -i all -b protein-dna -s common
```

All of the example toy-model-related commands are provided in pipeline.sh

**6) Aaaand, To visualize minnie results with Pymol**\
\
To be run in projectID/05_compare_complex/\*_freq_filtered/\*_freq_perres/complex_specific_


```
sed 's/,/ /g' *_hbonds_compared_*_perres.csv | grep "protein-dna" | grep ":" | awk '{printf "show sticks, (resi %s and chain %s) + (resi %s and chain %s) \n", $6,$2,$7,$3}' | sort -u > pymol_hbonds.pml
sed 's/,/ /g' *_hbonds_compared_*_perres.csv | grep "protein-dna" | grep ":" | awk '{printf "distance i. %s and n. %s and chain %s, i. %s and n. %s and chain %s\n", $6,$8,$2,$7,$9,$3}' | sort -u | sed 's/OP1/O1P/g' | sed 's/OP2/O2P/g' >> pymol_hbonds.pml
```

 Options to visualize each bond type are added to pipeline.sh !!

**7) Finally, to clean unnecessary folders**

 ```
minnie clean -cn 'projectID1' 'projectID2'
minnie clean -cn 'sox4' 'sox18'
 ```

## Troubleshoot
If you would need to see brief descriptions of the minnie's options

```
minnie --help
minnie <subcommand> --help
```
An example of help command - to retrieve usage information of the `minnie splitpdbs` option

```
$ minnie splitpdbs --help

Usage: minnie splitpdbs
                        -cn, --complexName     <string> <string>
                                                Project ID of your complex(s)

                        -p, --pdbs             [<traj.pdb>] [<traj.pdb>]
                                               Trajectory ID of your complex(s)

Usage example:

 minnie splitpdbs -cn sox4 sox18 -p sox4.pdb sox18.pdb
 minnie splitpdbs -cn sox4  -p sox4.pdb

```

If you would need to remove minnie from your conda setup
```
conda env remove -n  minnie
```
Also remove the minnie-related lines from your .bashrc

## Acknowledgements
We are grateful to Dr. João Rodrigues ([@JoaoRodrigues](https://github.com/JoaoRodrigues)) for making his interfacea package available for minnie! We also thank Tulay Karakulak ([@KarakulakTulay](https://github.com/KarakulakTulay)) for her contribution to the conceptualization of this work.

## Contact
ezgi.karaca@ibg.edu.tr
