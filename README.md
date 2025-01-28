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

If you're running any other version or distribution, there might be compatibility issues. 

### 1.2 Operating System

- **Linux-only**: This tool is designed to run on Linux-based operating systems. Please ensure that you're running a supported Linux distribution. This tool has been developed in:
  - **Ubuntu 22.04.5 LTS**
  - **64-bit OS** (this tool is not compatible with 32-bit systems)
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
python3.10 /dir/comb_calculator_MASTER/calculate_uniprot.py \
--mode fetch_full \
--output_file /dir/mapping_OUTPUT/uniprot_mapped_disorder_TEST.csv \
--cached_disorder_file /dir/comb_calculator_MASTER/Datasets/uniprot_mapped_disorder.csv
```
2. **UniProt Sample Mapping**
```bash
python3.10 /dir/comb_calculator_MASTER/calculate_uniprot.py \
--mode fetch_sample \
--sample_size 20 \
--output_file /dir/mapping_OUTPUT/uniprot_mapped_disorder_TEST.csv \
--cached_disorder_file /dir/comb_calculator_MASTER/Datasets/uniprot_mapped_disorder.csv 
```
3. **List of UIDs for Mapping**
```bash
python3.10 /dir/comb_calculator_MASTER/calculate_uniprot.py \
--mode map_existing \
--input_file /dir/comb_calculator_MASTER/example_inputs/uid_list_TEST.csv \
--output_file /dir/mapping_OUTPUT/uniprot_mapped_disorder_TEST.csv \
--cached_disorder_file /dir/comb_calculator_MASTER/Datasets/uniprot_mapped_disorder.csv
```

#### Positional Arguments

| Position | Parameter       | Values                 | Description                                    |
|----------|-----------------|------------------------|------------------------------------------------|
| 1        | Mode (--mode)   | `fetch_full`, `fetch_sample` , `map_existing` | String value indicating if the mapping is to be done for the full version of UniProt, a sample of it or for a provided list (correspondigly)             |
| 2        | Sample Size (--sample_size)   | Integers | Only for the `fetch_sample` mode. Indicates the desired sample of UniProt for mapping             |
| 3        | Input File (--input_file)     | `input_file.csv`      | String value indicating the directory of the input dataset containing UIDs in a column named `uid` (see examples) |
| 4        | Output File (--output_file)    | `output_file.csv`     | String value indicating the directory of the output dataset mapped for `sequence`, `refseq_id`, `disorder`|
| 5        | Disorder File (--cached_disorder_file)   | `uniprot_mapped_disorder.csv`     | String value indicating the directory of the dataset with the instrinsic disorder values. Please use the one provided in the examples if you aren't sure about your disorder values|

#### Inputs
Make sure you use the example inputs that are located in the `example_inputs` folder.

---

### 3.2 Script 2: `script2.py`

#### Example Command

```bash
python script2.py example_inputs/Dataset2 param1_value param2_value output_file2.csv
```

#### Positional Arguments

| Position | Parameter       | Values                 | Description                                    |
|----------|-----------------|------------------------|------------------------------------------------|
| 1        | Input Dataset   | `example_inputs/Dataset2` | Path to the input dataset folder.             |
| 2        | Parameter 1     | Custom value           | Specify the first parameter value (e.g., `0.5`). |
| 3        | Parameter 2     | Custom value           | Specify the second parameter value (e.g., `True`). |
| 4        | Output File     | `output_file2.csv`     | Name of the output file where results will be saved. |

#### Inputs
Make sure the required inputs are located in the `example_inputs/Dataset2` folder.

---

### 3.3 Script 3: `script3.py`

#### Example Command

```bash
python script3.py example_inputs/Dataset3 config.json output_file3.csv
```

#### Positional Arguments

| Position | Parameter       | Values                 | Description                                    |
|----------|-----------------|------------------------|------------------------------------------------|
| 1        | Input Dataset   | `example_inputs/Dataset3` | Path to the input dataset folder.             |
| 2        | Config File     | `config.json`          | Path to a configuration JSON file.            |
| 3        | Output File     | `output_file3.csv`     | Name of the output file where results will be saved. |

#### Inputs
Make sure the required inputs are located in the `example_inputs/Dataset3` folder, and provide the configuration file (`config.json`) in the same directory or specify the correct path.
## 4. Output: Mapping & Calculated Features
Describe the format and contents of the output produced by the tool. You might include details about the mapped data, how the features are calculated, and any interpretation of the results.

## 5. Acknowledgments
Give credit to any contributors, libraries, or other resources that helped in the creation of the tool. You can also mention any external resources or research papers used in your project.
