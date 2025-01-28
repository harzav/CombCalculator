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
Provide instructions on how to run the tool, including example commands and their expected results. You can also give sample input and output for clarity.

## 4. Output: Mapping & Calculated Features
Describe the format and contents of the output produced by the tool. You might include details about the mapped data, how the features are calculated, and any interpretation of the results.

## 5. Acknowledgments
Give credit to any contributors, libraries, or other resources that helped in the creation of the tool. You can also mention any external resources or research papers used in your project.
