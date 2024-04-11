# q_dms_ttr_paper

make sure to download the data used in this project from 10.6084/m9.figshare.25331758
unzip this data and store it in this current directory as data/ <br>

## how to install 

```bash
git clone 
```

## Layout of directory

q_dms_ttr_paper/<br>
outside information not from this project<br>
	resources/ <br>
		csvs/ 
			p5_sequences.csv  			# has all the 5' sequences so can be removed
			sequencing_runs.csv			# has all the sequencing runs required for this paper
			ttr_mutation_dgs_all.csv	# has all dG data from steves ttrs 
			ttr_mutation_dgs_subset.csv # has only dG data from steves ttr used in this study
	data/
		sequencing_runs/ # all runs required for this 

## figure generation 

All final figures were generated with notebooks 
notebooks/final_figures/



## commands to run
All commands should be run from the root directory of the project. All data should be stored in the data/ directory in the root 
directory. If its not code will not work properly.

### designing the constructs for the study 
```python
# step 1
# takes the data from bonilla et al. and selects the 98 that will be used in this study
# data from bonilla et al. is in data/construction_data/bonilla_et_al_pnas_2021
# construct names are reversed to be consistent with this study CAUGG_CCUAAA to CCUAAA_CAUGG
# sets will be outputed to data/constructs_data/sets
python q_dms_ttr_paper/cli.py generate-sets

#step 2 
# takes the sets of TLR mutants selected in step 1 and inserts each into the miniTTR scaffold in replace of 
# the wild-type 11ntR sequence and structure.
# sets will be outputed to data/constructs_data/scaffold_sets
python q_dms_ttr_paper/cli.py generate-scaffold-sets

# step 3 
# takes the sets of step 2 and randomizes their helical regions to increase the diversity of the library.
# sets will be outputed to data/constructs_data/randomized_sets
python q_dms_ttr_paper/cli.py randomize-helices

# step 4 
python q_dms_ttr_paper/cli.py barcode-libraries

# step 5 
python q_dms_ttr_paper/cli.py generate-order

```

## data anaylsis the qDMS-MaPseq data / mutation data / solvent accessibility of structures

```python
# move sequencing data into correct folder
python q_dms_ttr_paper/cli.py get-sequencing-data




# generate initial dataframes with information about each motif
python q_dms_ttr_paper/cli.py process-sequencing-data
# generate information about mutations of each ttr construct 
python q_dms_ttr_paper/cli.py characterize-mutations
# generate the mg 1/2 of all ttr mutants
python q_dms_ttr_paper/cli.py process-mg-1-2

```

----------

## generating plots 




## generating the figures for the paper 



