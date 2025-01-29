# CombCalculator

CombCalculator is an effort for a light and easy-downloaded command-line tool for UniProt ID mapping and combination generation. Through this tool the user can fully map lists of UniProt IDs, fetch the most recent and mapped version of whole UniProt, and generate combinations of UniProt IDs for protein-protein interaction (PPI) studies. For all PPI combinations a highly-curated set of 68 interaction-informative features is calculated, including an interaction prediction and interaction affinity estimation (for more info on how these prediction models were created please check the 'TR_PPI' github<sup>1</sup>).

The tool consists of 3 bash commands that cover different mapping and combination generation needs:

 <code>comb_calculator.py</code> ,  <code>calculate_uniprot.py</code> ,  <code>feature_calculation.py</code> 

For understanding which tool to use according to your input and needs please follow carefully the following flowchart!

<h2>A. CombCalculator Flowchart</h2>

## <h2>B. Table of Contents - Tutorial for Setting Up and Running the Tool</h2>

- [1. Requirements](#1-requirements)
- [2. Setup](#2-setup)
- [3. Running the Tool: Commands & Examples](#3-running-the-tool-commands--examples)
- [4. Output: Mapping & Calculated Features](#4-output-mapping--calculated-features)
- [5. Acknowledgments](#5-acknowledgments)

## 1. Requirements
To run this tool, you'll need the following:

### 1.1 Python Version

- **Python version**: 3.10.12
  - Make sure to use Python version 3.10.12 for compatibility with the tool. You can check your Python version by running:
```bash
python --version
```

If you're running any other version or distribution, there might be compatibility issues. 

### 1.2 Operating System

- **Linux-only**: This tool is designed to run on Linux-based operating systems. Please ensure that you're running a supported Linux distribution. This tool has been developed in:
  - **Ubuntu 22.04.5 LTS**
  - **64-bit OS** 
  - **GNOME version 42.9**
  - **Windowing system: X11**

## 2. Setup
### 2.1 Installing Dependencies
The required libraries are listed in the <code>requirements.txt</code> file located in the root of this repository. To install these dependencies, follow these steps:

1. **Clone the repository:**
```bash
git clone https://github.com/harzav/CombCalculator.git
cd CombCalculator
```

2. **Set up a virtual environment (recommended for managing dependencies):**
If you're using venv, do the following steps:
- **Create a virtual environment:**
```bash
python -m venv venv
```

- **Activate the virtual environment:**
```bash
source venv/bin/activate
```

3. **Install the required libraries, through the virtual environment (preferably):**
Once the virtual environment is activated, install all the required libraries using <code>pip</code> by running:
```bash
pip install -r requirements.txt
```

This will install all the Python dependencies listed in the <code>requirements.txt</code> file.
If you're using a system-wide Python installation without a virtual environment (not recommended), you can install the dependencies directly by running:
```bash
pip install --user -r requirements.txt
```

### 2.2 Downloading the Datasets
This tool requires an additional folder named <code>Datasets</code> to be downloaded from Google Drive and placed in the project's directory. Follow these steps:

1. **Click the following Google Drive link to download the required file:**
[Download File
](https://drive.google.com/drive/folders/1dog4Z5EaBjlzVMcWjLDAXrziYt2yQeL2?usp=drive_link
)
 
2. **Save the downloaded file to the root directory of the project.**
 
If you have downloaded a <code>Datasets</code> folder to your <code>Downloads</code> directory, you need to move it to the root of the project directory. Use the following command:
```bash
mv ~/Downloads/Datasets /path-to-project/CombCalculator
```
 
If your CombCalculator is located in your <code>home</code> directory under <code>CombCalculator</code>, the command would look like this:
```bash
mv ~/Downloads/Datasets ~/CombCalculator
```
**Important!!!**
Make sure the file name and location match the tool's expected setup.
Your directory structure should look like this:
```scss
CombCalculator/
├── requirements.txt
├── Datasets
├── (other files)
```

## 3. Running the Tool: Commands & Examples

The `CombCalculator` project provides three Python scripts to perform the tasks detailed in the `CombCalculator` flowchart. Below are the commands for each script, their example usage, and details for each argument.

### 3.1 `calculate_uniprot.py`

This command is used for mapping either the full UniProt Database, either a random sample of 'n' proteins from UniProt or either a provided list of UniProt IDs. The different example commands for each case are provided below
#### Example Commands
1. **Full UniProt Mapping**
```bash
python3.10 /dir/CombCalculator/calculate_uniprot.py \
--mode fetch_full \
--output_file /dir/mapping_OUTPUT/uniprot_mapped_disorder_TEST.csv \
--cached_disorder_file /dir/CombCalculator/Datasets/uniprot_mapped_disorder.csv
```
2. **UniProt Sample Mapping**
```bash
python3.10 /dir/CombCalculator/calculate_uniprot.py \
--mode fetch_sample \
--sample_size 20 \
--output_file /dir/mapping_OUTPUT/uniprot_mapped_disorder_TEST.csv \
--cached_disorder_file /dir/CombCalculator/Datasets/uniprot_mapped_disorder.csv 
```
3. **List of UIDs for Mapping**
```bash
python3.10 /dir/CombCalculator/calculate_uniprot.py \
--mode map_existing \
--input_file /dir/CombCalculator/example_inputs/uid_list_TEST.csv \
--output_file /dir/mapping_OUTPUT/uniprot_mapped_disorder_TEST.csv \
--cached_disorder_file /dir/CombCalculator/Datasets/uniprot_mapped_disorder.csv
```

#### Positional Arguments

| Position | Parameter       | Values                 | Description                                    |
|----------|-----------------|------------------------|------------------------------------------------|
| 0        | Script directory    | `/dir/CombCalculator/calculate_uniprot.py`     | String value indicating the script directory. Remove `/dir` if it is saved in your home directory|
| 1        | Mode (--mode)   | `fetch_full`, `fetch_sample` , `map_existing` | String value indicating if the mapping is to be done for the full version of UniProt, a sample of it or for a provided list (correspondigly)             |
| 2        | Sample Size (--sample_size)   | Integer, n ∈ ℕ ∩ [2, +∞] | Only for the `fetch_sample` mode. Indicates the desired sample of UniProt for mapping             |
| 3        | Input File (--input_file)     | `/dir/CombCalculator/example_inputs/uid_list_TEST.csv`      | String value indicating the directory of the input dataset containing UIDs in a column named `uid` (see examples) |
| 4        | Output File (--output_file)    | `/dir/mapping_OUTPUT/uniprot_mapped_disorder_TEST.csv`     | String value indicating the directory of the output dataset mapped for `sequence`, `refseq_id`, `disorder`|
| 5        | Disorder File (--cached_disorder_file)   | `/dir/CombCalculator/Datasets/uniprot_mapped_disorder.csv`     | String value indicating the directory of the dataset with the instrinsic disorder values. Please use the one provided in the `CombCalculator/Datasets` directory if you aren't sure about your disorder values|

---

### 3.2 `comb_calculator.py`

The main function. Receives a list of mapped UniProt IDs (according to the mapping detailed in the flowchart, `Section 4. Output: Mapping & Calculated Features` and the `calculate_uniprot.py` function) and then calculates all their possible combinations of PPIs and a set of descriptive features for them (also see `Section 4. Output: Mapping & Calculated Features` for their description)

**Important!!** To use this tool the exact mapping mentioned is necessary. It is best if you use the `calculate_uniprot.py` function in case you have a list of UIDs, to ensure proper mapping.

#### Example Command

```bash
python3.10 /dir/CombCalculator/comb_calculator.py \
--input /dir/CombCalculator/example_inputs/uniprot_mapped_disorder_TEST.csv \
--output /dir/callable_output_TEST/ \
--max_combinations 1000000 \
--threads 4 \
--batch_size 5 \
--additional /dir/CombCalculator/example_inputs/combinationsTEST.csv
```

#### Positional Arguments

| Position | Parameter       | Values                 | Description                                    |
|----------|-----------------|------------------------|------------------------------------------------|
| 0        | Script directory    | `/dir/CombCalculator/comb_calculator.py`     | String value indicating the script directory. Remove `/dir` if it is saved in your home directory|
| 1        | Input Dataset (--input)   | `/dir/CombCalculator/example_inputs/uniprot_mapped_disorder_TEST.csv` | String value indicating the directory of the input dataset containing UIDs in a column named `uid`, and mapped for  `sequence`, `refseq_id`, `disorder` (see examples)          |
| 2        | Output Directory (--output)     | `/dir/callable_output_TEST/`           | String value indicating the output directory, where all the datasets will be saved, containing all the possible PPIs that the input UIDs can create, along with 68 descriptive features (see `Section 4. Output: Mapping & Calculated Features` for their description)|
| 3        | Max. number of combinations (--max_combinations)     | Integer, n ∈ ℕ ∩ [2, +∞] | The maximum number of combinations you want to create. You *should* set a number safely above the maximum possible PPI combinations if you want all of them calculated.  |
| 4        | Threads for parallel processing (--threads)    | Integer, n ∈ ℕ ∩ [1, maximum computational resources]     | The number of threads you want to allocate for processing the PPI features for each combinations batch. Depends on your computational resources. |
| 5        | Batch size of PPIs for parallel processing (--threads)    | Integer, n ∈ ℕ ∩ [1, +∞]     | The number of PPIs per processing batch. Should be similar to `--threads` for efficiency. Depends on your computational resources. |
| 6        | Additional dataset directory (--additional)    | `/dir/CombCalculator/example_inputs/combinationsTEST.csv`     | String value indicating the directory of an additional dataset in which you already have calculated PPIs (they will be excluded from the calculation). You can always use the dataset provided: `example_inputs/combinationsTEST.csv`, which has already calculated PPIs and in which you can add more, newly calculated PPIs (*keeping always the same structure!!*) |

---

### 3.3 `feature_calculation.py`

#### Example Command

```bash
python3.10 /dir/CombCalculator/feature_calculation.py \
/dir/CombCalculator/example_inputs/fc_TEST.csv \
/dir/fc_OUTPUT/OUTPUT.csv \
/dir/fc_OUTPUT/OUTPUT_RAW.csv 
```

#### Positional Arguments

| Position | Parameter       | Values                 | Description                                    |
|----------|-----------------|------------------------|------------------------------------------------|
| 0        | Script directory    | `/dir/CombCalculator/feature_calculation.py`     | String value indicating the script directory. Remove `/dir` if it is saved in your home directory|
| 1        | Input Dataset   | `/dir/CombCalculator/example_inputs/fc_TEST.csv` | String value indicating the directory of the input dataset containing PPIs in two columns named `uidA`, `uidB` |
| 2        | Output File     | `/dir/fc_OUTPUT/OUTPUT.csv`           | String value indicating the directory of the output dataset, with all the possible mapping and features calculated for the PPIs that were provided (see `Section 4. Output: Mapping & Calculated Features` for their description). *Note*: This dataset is scaled using Arithmetic Sample-Wise normalization and imputated using KNN imputation method. If this is not wanted, use `Raw Output File` instead. |
| 3        | Raw Output File | `/dir/fc_OUTPUT/OUTPUT_RAW.csv`           | String value indicating the directory of the output dataset, with without the preprocessing methods applied. |

---

#### Inputs
For a test drive, make sure you use the appropriate inputs that are located in the `CombCalculator/example_inputs/` folder. All the directories in the example commands use a `/dir` prefix in the beginning.  Remove it if `CombCalculator` is  saved in your home directory.
Please be aware of the different input formats for each of the 3 scripts and make sure you keep the same feature names for the mandatory rows. To use `comb_calculator.py` tool the exact mapping mentioned is necessary. It is best if you use the `calculate_uniprot.py` function in case you have a list of UIDs, to ensure proper mapping.Below are some screenshots of input examples for each tool:

1. `calculate_uniprot.py`
2. `comb_calculator.py`
3. `feature_calculation.py`

## 4. Output: Mapping & Calculated Features

The output of each command is either a mapped list of UIDs (in the case of `calculate_uniprot.py`) or a set of 68 calculated features for a protein interacting pair (in the case of `comb_calculator.py` and `feature_calculation.py`)
The features for each case will be explained below.

### 4.1 Mapping Features

| Feature Name    | Values                 | Description                                    |
|-----------------|------------------------|------------------------------------------------|
| `uid`    | String    | The UniProt ID (UID) of the protein|
| `refseq_id`   | String | The NCBI Reference Sequence (RefSeq) database ID of each protein |
| `seq`     | String           | The protein sequence |
| `disorder` | Float, n ∈ ℝ ∩ [0, 1]           | The intrinsic disorder value for the protein. It is predicted using the DisPredict 3.0 Tool <sup>2</sup> for the whole UniProt Database. For new entries in the UniProt Database (beyond Sep. 2024) it is predicted using the IUPred2A Tool  <sup>3</sup>|

## 5. Acknowledgments
Give credit to any contributors, libraries, or other resources that helped in the creation of the tool. You can also mention any external resources or research papers used in your project.

## Citations

1) https://github.com/harzav/TR_PPI_project

2) Kabir, M. W. U., & Hoque, M. T. (2024). DisPredict3.0: Prediction of intrinsically disordered regions/proteins using protein language model. Applied Mathematics and Computation, 472, 128630. https://doi.org/10.1016/j.amc.2024.128630

3) Bálint Mészáros, Gábor Erdős, Zsuzsanna Dosztányi (2018), IUPred2A: context-dependent prediction of protein disorder as a function of redox state and protein binding, Nucleic Acids Research ;46(W1):W329-W337.

