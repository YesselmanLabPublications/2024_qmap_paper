# 2024_qmap_paper

## Data Download
Make sure to download the data used in this project from [10.6084/m9.figshare.25331758](https://figshare.com/articles/dataset/data/25331758). Unzip this data and store it in the `data/` directory within the current project directory.

## Installation

```bash
# Create a new conda environment
conda create -n qmap_paper python=3.8

# Clone the repository
git clone https://github.com/YesselmanLabPublications/2024_qmap_paper.git

# Navigate to the project directory
cd 2024_qmap_paper

# Install the required packages
pip install .

# Download the data from Figshare and unzip it in the "2024_qmap_paper" directory
# Note: The data must be placed correctly for the code to work properly
wget -O data.zip https://figshare.com/ndownloader/files/44981896
unzip data.zip 

# If you do not have RNAfold installed, you can install it using conda
# Install ViennaRNA using conda
conda install -c bioconda viennarna
```

## Figure generation 
All final figures were generated with notebooks `notebooks/figures`

## Commands to Run
All commands should be run from the root directory of the project. Ensure all data is stored in the `data/` directory in the root directory. If not, the code will not function correctly.

### Designing the Constructs for the Study
This will regenerate constructs for the study. This step is not necessary as the data is already generated and stored in the `data/constructs_data/` directory. Running this will overwrite existing constructs with different helical sequences.

#### Step 1: Data Selection
Select a subset of data from the study conducted by Bonilla et al. to be used in the current study.

- The data from Bonilla et al. is located in `data/construction_data/bonilla_et_al_pnas_2021`.
- Construct names in the data are reversed to maintain consistency with the current study (e.g., `CAUGG_CCUAAA` becomes `CCUAAA_CAUGG`).
- Selected data sets will be outputted to `data/constructs_data/sets`.

To execute this step, run:
```bash
python qmap_paper/cli.py generate-sets
```

#### Step 2: Inserting TLR Mutants into miniTTR Scaffold
Insert each TLR mutant into the miniTTR scaffold, replacing the wild-type 11ntR sequence and structure. The resulting sets will be outputted to `data/constructs_data/scaffold_sets`.

To execute this step, run:
```bash
python qmap_paper/cli.py generate-scaffold-sets
```

#### Step 3: Randomizing Helical Regions
Randomize the helical regions of the sets generated in Step 2 to increase the diversity of the library. The randomized sets will be outputted to `data/constructs_data/randomized_sets`.

To execute this step, run:
```bash
python qmap_paper/cli.py randomize-helices
```

#### Step 4: Barcoding Libraries
Add additional helical barcodes to the end of each construct to increase hamming distance.

To execute this step, run:
```bash
python qmap_paper/cli.py barcode-libraries
```

#### Step 5: Generating Order
Generate the order for the constructs.

To execute this step, run:
```bash
python qmap_paper/cli.py generate-order
```

## Data Analysis: DMS-MaPseq Data / Mutation Data / Solvent Accessibility of Structures

### Move DMS-MaPseq Data into Correct Folder
To move the sequencing data into the correct folder.
Note this will not work if you do not have access to YesselmanLab data. All of this data
is available already on via the figshare link above.

```bash
python qmap_paper/cli.py get-sequencing-data
```

### Process DMS-MaPseq Data
Processes the DMS-MaPseq data to assign mutation fractions to specific motifs.
```bash
python qmap_paper/cli.py process-sequencing-data
```

### Characterize Mutations
Assigns what mutations occurs in what TLR variant
```bash
python qmap_paper/cli.py characterize-mutations
```

### Compute Solvent Accessibility
Compute the solvent accessibility of of the three As in the tetraloop receptor
```bash
python qmap_paper/cli.py compute-sasa
```

### Compute the Mg1/2 of All TTR Mutants
Fits the tetraloop As as a function of mg to compute the mg1/2
```bash
python qmap_paper/cli.py process-mg-1-2
```

## generating plots 




