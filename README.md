# 2024_qmap_paper

make sure to download the data used in this project from [10.6084/m9.figshare.25331758](https://figshare.com/articles/dataset/data/25331758) unzip this data and store it in this current directory as data/ <br>

## how to install 

```bash
conda create -n qmap_paper python=3.8
git clone https://github.com/YesselmanLabPublications/2024_qmap_paper.git
cd 2024_qmap_paper
pip install .
# go to figshare and download the data and unzip it in the "2024_qmap_paper" directory
# NOTE if you do not put the data there things will not work properly
```

## figure generation 
All final figures were generated with notebooks `notebooks/`

## commands to run
All commands should be run from the root directory of the project. All data should be stored in the data/ directory in the root 
directory. If its not code will not work properly.

### designing the constructs for the study 
This will regenerate constructs for the study. This is not necessary to run as the data is already generated and stored in the data/constructs_data/ directory. It will overwrite existing constructs with different helical sequences.

#### step 1 

This step involves selecting a subset of data from the study conducted by Bonilla et al. The selected data will be used in the current study.

- The data from Bonilla et al. is located in the directory `data/construction_data/bonilla_et_al_pnas_2021`.
- The construct names in the data are reversed to maintain consistency with the current study. For example, `CAUGG_CCUAAA` is reversed to `CCUAAA_CAUGG`.
- The selected data sets will be outputted to the directory `data/constructs_data/sets`.

To execute this step, run the following command:
```python
# Step 1: Data Selection
python q_dms_ttr_paper/cli.py generate-sets
```

#### step 2

This step involves taking the sets of TLR mutants selected in Step 1 and inserting each mutant into the miniTTR scaffold, replacing the wild-type 11ntR sequence and structure. The resulting sets will be outputted to the directory `data/constructs_data/scaffold_sets`.

To execute this step, run the following command:

```python
# Step 2: Inserting TLR Mutants into miniTTR Scaffold
python q_dms_ttr_paper/cli.py generate-scaffold-sets
```
#### step 3

```python
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



