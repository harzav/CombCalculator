# CombCalculator

# CombCalculator

CombCalculator is an effort for a light and easy-downloaded command-line tool for UniProt ID mapping and combination generation. Through this tool the user can fully map lists of UniProt IDs, fetch the most recent and mapped version of whole UniProt, and generate combinations of UniProt IDs for protein-protein interaction (PPI) studies. For all PPI combinations a highly-curated set of 68 interaction-informative features is calculated, including an interaction prediction and interaction affinity estimation (for more info on how these prediction models were created please check the 'TR_PPI' github<sup>1</sup>).

The tool consists of 3 bash commands that cover different mapping and combination generation needs:

 <code>comb_calculator.py</code> ,  <code>calculate_uniprot.py</code> ,  <code>feature_calculation.py</code> 

For understanding which tool to use according to your input and needs please follow carefully the following flowchart!

<h2>A. CombCalculator Flowchart</h2>

![flowchart](file:///c%3A/Users/harry/Downloads/Comb%20%281%29.png)

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

If you don't have Python 3.10.12 installed, you can download it from the official Python website.

### 1.2 Installing Dependencies
The required libraries are listed in the <code>requirements.txt</code> file located in the root of this repository. To install these dependencies, follow these steps:

1. **Clone the repository:**
   ```bash
   git clone https://github.com/your-username/your-repository.git
   cd your-repository
   
3.  
## 2. Setup
In this section, you will explain how to install and set up the tool. This can include steps such as cloning the repository, installing dependencies, and configuring the environment.

## 3. Running the Tool: Commands & Examples
Provide instructions on how to run the tool, including example commands and their expected results. You can also give sample input and output for clarity.

## 4. Output: Mapping & Calculated Features
Describe the format and contents of the output produced by the tool. You might include details about the mapped data, how the features are calculated, and any interpretation of the results.

## 5. Acknowledgments
Give credit to any contributors, libraries, or other resources that helped in the creation of the tool. You can also mention any external resources or research papers used in your project.
