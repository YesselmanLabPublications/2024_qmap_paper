# q_dms_ttr_paper

Layout of directory

q_dms_ttr_paper/
	# outside information not from this project
	resources/ 
		csvs/ 
			p5_sequences.csv  			# has all the 5' sequences so can be removed
			sequencing_runs.csv			# has all the sequencing runs required for this paper
			ttr_mutation_dgs_all.csv	# has all dG data from steves ttrs 
			ttr_mutation_dgs_subset.csv # has only dG data from steves ttr used in this study
	data/
		sequencing_runs/ # all runs required for this 

```python
# move sequencing data into correct folder
python q_dms_ttr_paper/cli.py get-sequencing-data
# generate initial dataframes with information about each motif
python q_dms_ttr_paper/cli.py process-sequencing-data
# generate information about mutations of each ttr construct 
python q_dms_ttr_paper/cli.py characterize-mutations
# generate the mg 1/2 of all ttr mutants
python q_dms_ttr_paper/cli.py process-mg-1-2

----------
# generate figures


```
