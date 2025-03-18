# CombCalculator

CombCalculator is an effort for a light and easy-downloaded command-line tool for UniProt ID mapping and combination generation. Through this tool the user can fully map lists of UniProt IDs, fetch the most recent and mapped version of whole UniProt, and generate combinations of UniProt IDs for protein-protein interaction (PPI) studies. For all PPI combinations a highly-curated set of 68 interaction-informative features is calculated, including an interaction prediction and interaction affinity estimation (for more info on how these prediction models were created please check the [TR_PPI GitHub](https://github.com/harzav/TR_PPI_project)<sup>[1]</sup>).

The tool consists of 3 bash commands that cover different mapping and combination generation needs:

 <code>comb_calculator.py</code> ,  <code>calculate_uniprot.py</code> ,  <code>feature_calculation.py</code> 

For understanding which tool to use according to your input and needs please follow carefully the following flowchart!

<h2>A. CombCalculator Flowchart</h2>

<img src="https://github.com/user-attachments/assets/e4d6e5f7-f35f-41d5-8d7a-0b35ccee8665" width="1000">

## <h2>B. Table of Contents - Tutorial for Setting Up and Running the Tool</h2>

- [1. Requirements](#1-requirements)
- [2. Setup](#2-setup)
- [3. Running the Tool: Commands & Examples](#3-running-the-tool-commands--examples)
- [4. Output: Mapping & Calculated Features](#4-output-mapping--calculated-features)
- [5. Acknowledgments](#5-acknowledgments)
- [6. Citations](#6-citations)

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
>[!IMPORTANT] 
> Make sure the file name and location match the tool's expected setup.
> Your directory structure should look like this:
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

>[!WARNING] 
>To use this tool the exact mapping mentioned is necessary. It is best if you use the `calculate_uniprot.py` function in case you have a list of UIDs, to ensure proper mapping.

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

#### Feature Calculation Multiprocessing

For calculating features for a large number of PPIs (e.g. > 10k) it is best to use the **multiprocessing** version of Feature Calculation function named `feature_calculation_multiprocessing.py`. This function achieves the same results with the default function, but processes the PPIs in batches of size *N* set by the user, utilizing a number of cores for parallel processing of the batches, which is manually set based by your systems capabilities.

The example command of calling the script is:

```bash
python3.10 /dir/CombCalculator/feature_calculation_multiprocessing.py \
--input_csv /dir/CombCalculator/example_inputs/fc_TEST.csv \
--batch_size 10 \
--num_threads 4 \
--output_dir /dir/fc__multi_OUTPUT
```

And the Positional arguments are the following:

| Position | Parameter       | Values                 | Description                                    |
|----------|-----------------|------------------------|------------------------------------------------|
| 0        | Script directory    | `/dir/CombCalculator/feature_calculation_multiprocessing.py`     | String value indicating the script directory. Remove `/dir` if it is saved in your home directory|
| 1        | Input Dataset (--input_csv)  | `/dir/CombCalculator/example_inputs/fc_TEST.csv` | String value indicating the directory of the input dataset containing PPIs in two columns named `uidA`, `uidB` |
| 2        | Batch Size (--batch_size) |  Integer, n ∈ ℕ ∩ [10, +∞] | The number of PPIs per batch. Minimum is recommended to be 10 for efficient feature calculation. Maximum is based on your systems capabilities and teh total number of PPIs in your input dataset. |
| 3        | Threads for parallel processing (--num_threads)    | Integer, n ∈ ℕ ∩ [4, maximum computational resources]     | The number of threads you want to allocate for processing the PPI features for each combinations batch. Depends on your computational resources. Minumum is set to be 4- less than 4 is inefficient for parallel feature calculation|
| 4        | Output File (--output_dir)    | `/dir/fc__multi_OUTPUT`           | String value indicating the directory of the output folder. There all the batches and final output CSVs (`/OUTPUT.csv` & `OUTPUT_RAW.csv`) will be saved  |

---

#### Inputs

For a test drive, make sure you use the appropriate inputs that are located in the `CombCalculator/example_inputs/` folder. All the directories in the example commands use a `/dir` prefix in the beginning.  Remove it if `CombCalculator` is  saved in your home directory.
Please be aware of the different input formats for each of the 3 scripts and make sure you keep the same feature names for the mandatory rows. To use `comb_calculator.py` tool the exact mapping mentioned is necessary. It is best if you use the `calculate_uniprot.py` function in case you have a list of UIDs, to ensure proper mapping. 

Below are some screenshots of input examples for each tool:

1. `calculate_uniprot.py`

![calc_uniprot](https://github.com/user-attachments/assets/14eac732-a4ee-4ac4-9012-05a37f7446e4)

2. `comb_calculator.py`

![comb_calculator](https://github.com/user-attachments/assets/32ee6e38-5769-426e-8359-36d62d6cfae9)

3. `feature_calculation.py`

![feature_calc](https://github.com/user-attachments/assets/e716222a-1f16-471c-9636-fd42a5391334)

## 4. Output: Mapping & Calculated Features

The output of each command is either a mapped list of UIDs (in the case of `calculate_uniprot.py`) or a set of 68 calculated features for a protein interacting pair (in the case of `comb_calculator.py` and `feature_calculation.py`). In the case of the PPI datasets, the mapping features are calculated for both proteins in the PPI (with suffixes `_A`, and `_B`) 
The features for each case will be explained below.

### 4.1 Mapping Features

| Feature Name    | Values                 | Description                                    |
|-----------------|------------------------|------------------------------------------------|
| `uid`    | String    | The UniProt ID (UID) of the protein|
| `refseq_id`   | String | The NCBI Reference Sequence (RefSeq) database ID of each protein |
| `seq`     | String           | The protein sequence |
| `disorder` | Float, n ∈ ℝ ∩ [0, 1]           | The Intrinsic Disorder value for the protein. It is predicted using the [DisPredict 3.0 Tool](https://doi.org/10.1016/j.amc.2024.128630) <sup>[2]</sup> for the whole UniProt Database. For new entries in the UniProt Database (beyond Sep. 2024) it is predicted using the [IUPred2A Tool](https://doi.org/10.1093/nar/gky384)<sup>[3]</sup>|

---

### 4.2 Calculated Features

| Feature Role | Number of Columns | Column Names | Description | Values |
|---------|------------------|--------------|-------------|--------|
| **PPI's Disorder Prediction** | 1 | `disorder_check` | Indicates if the Intrinsic Disorder of the pair is above average, in which case features can't be calculated. If one of the two proteins in the PP has an Intrinsic Disorder >0.5, then the pair is considered Disordered (`1`), else not (`0`)   | 0, 1 |
| **GO term similarity** | 3 | `BP_similarity`, `MF_similarity`, `CC_similarity` | The similarity of the proteins in the pair based on the similarity of their respective GO terms. Range is from 0-1, where `0` is non-similarity, `1` is total similarity and  | Float, n ∈ ℝ ∩ [0, 1] , NaN |
| **Existence of a Homologous Interacting Pair in other organisms** | 4 | `Homologous in Mouse`, `Homologous in Drosophila`, `Homologous in Yeast`, `Homologous in Ecoli` | Existence of corresponding homologous PPIs in Mouse, Yeast, Drosophila, E. coli. If a homologous pair exists, it is calculated whether this pair is also an interacting pair, and if it is, `1` is appended in this column, else `0` is appended. `NaN` values are appended to the interacting pairs with no homology between Human and each other examined organism. The PPI Datasets of other organisms are derived from the DIP database, while the mapping of homologous proteins was done via the NCBI- HomoloGene Dataset. | 0, 1, NaN |
| **Existence of the PPI in other Databases** | 4 | `Exists in DIP?`, `Exists in APID?`, `Exists in BIOGRID?`, `Exists in MINT?` | Existence of PPI in APID, DIP, BIOGRID, and MINT databases. `0` if it doesn't exist, `1` if it exists | 0, 1 |
| **Sequence Similarity** | 1 | `Sequence_similarity` | Sequence SImilarity E-value of the sequences of the proteins in each PPI. Lower (i.e., stronger) E-values indicate more “significant” alignments, suggesting a higher probability that the sequences share a common evolutionary origin. Higher (i.e., weaker) E-value indicates that the alignment might be a random event. | Float, n ∈ ℝ ∩ [-∞, +∞] |
| **Domain Interactions** | 1 | `pfam_interaction` | Presence of known domain interactions between the proteins in each PPI pair based on the Pfam Database. `0` if a domain interaction doesn't exist, `1` if it exists, and `NaN` if there is no info on the DB. | 0, 1, NaN |
| **Existence in iRefIndex** | 1 | `irefindex_check` | Presence of the PPI in [iRefIndex](https://doi.org/10.1186/1471-2105-9-405)<sup>[4]</sup> Database. Specifically it checks if the PPI is a unique, binary, human and experimentally verified PPI in iRefIndex. `0` if it isn't, `1` if it is.  | 0, 1 |
| **iRefIndex Method of Identification** | 1 | `irefindex_method` | The experimental method of identification in iRefIndex. For more info on these methods check [iRefIndex MITAB Documentation](https://irefindex.vib.be/wiki/index.php/README_MITAB2.6_for_iRefIndex_19.0#Column_number:_7_.28Method.29)<sup>[5]</sup>. Applies if `1` is appended in `irefindex_check` column. Example: `two hybrid prey pooling approach` | String, NaN |
| **iRefIndex PPI Detection Citation** | 1 | `irefindex_pmid` | The PMID of the paper for the experimental PPI identification. Example: `pubmed:23275563`| String, NaN |
| **Subcellular co-localization** | 1 | `Subcellular Co-localization?` | Subcellular co-localization of the two proteins in the pair in eukaryotic cells.  `0` if a co-localization doesn't exist, `1` if it exists, and `NaN` if there is no info. | 0, 1, NaN |
| **Gene expression profile similarity** | 15 | `0`, `1` , ...., `14` | The Spearman Index of the similarity of the two proteins in terms of their Gene Expression among fifteen NCBI GEO gene expression datasets | Float, n ∈ ℝ ∩ [-1, 1] |
| **Amino acid difference** | 20 | `A%`, ...., `G %` | The absolute difference in the percentage of every amino acid between the PP | Float, n ∈ ℝ ∩ [0, 1] |
| **Molecular weight difference** | 1 | `MW dif` | The absolute difference in molecular weight between PP | Float, n ∈ ℝ ∩ [0, +∞] |
| **Aromaticity index difference** | 1 | `Aromaticity dif` | The absolute difference in aromaticity index difference between PP. | Float, n ∈ ℝ ∩ [0, 1] |
| **Instability index difference** | 1 | `Instability dif` | The absolute difference in instability index between PP | Float, n ∈ ℝ ∩ [0, 100] |
| **Amino acid fraction difference** | 3 | `helix_fraction_dif`, `turn_fraction_dif`, `sheet_fraction_dif` | The difference in fraction of total amino acids that are contained in 3 areas: The fraction of aa in helix, the fraction of aa in turn, the fraction of aa in sheet | Float, n ∈ ℝ ∩ [0, 1] |
| **Molar extinction coefficient difference** | 2 | `cys_reduced_dif`, `cys_residues_dif` | The difference in molar extinction coefficient when: (a) The molar extinction coefficient is calculated assuming cysteines(reduced) and (b) The molar extinction coefficient is calculated assuming cystines residues (Cys-Cys-bond) | Integer, n ∈ ℕ ∩ [0, +∞] |
| **GRAVY (Grand Average of Hydropathy) difference** | 1 | `gravy_dif` | The absolute difference in GRAVY index between PP | Float, n ∈ ℝ ∩ [0, +∞]  |
| **pH charge difference** | 1 | `ph7_charge_dif` | The absolute difference in the protein charge when pH= 7 | Float, n ∈ ℝ ∩ [0, +∞]  |
| **RNA expression profile similarity** | 2 | `GSE227375_spearman`, `GSE228702_spearman` | The Spearman Index of the similarity between each protein in the pair, in terms of their RNA expression profiles (NCBI GEO GSE227375 and GSE228702 datasets) | Float, n ∈ ℝ ∩ [-1, 1] |
| **Classifier Prediction** | 1 | `class_value` | The predicted class of the PPI. It is calculated using an [Evolutionary Optimization Algorithm](https://github.com/harzav/TR_PPI_project)<sup>[1]</sup> , along with SVM and RF classifiers| 0, 1  |
| **Classifier Confidence** | 1 | `classifier_confidence` | The confidence of the Classification Algorithm in its prediction | Float, n ∈ ℝ ∩ [0, 1]  |
| **Interaction Affinity Estimation** | 1 | `regression_value` | The estimation of the interaction's Affinity. It is calculated using an [Evolutionary Optimization Algorithm](https://github.com/harzav/TR_PPI_project)<sup>[1]</sup>, along with SVR, RBF-SVR, RFR, and CNN regressors | Float, n ∈ ℝ ∩ [0, 1]  |

---

>[!NOTE]
>- "PP" refers to the "Protein Pair" in each PPI
>- "NaN" indicates missing or unavailable values.
>- Spearman index measures correlation between values, ranging from -1 to 1.
>- For more info on how the Classifier and the Interaction Affinity Estimation models were created please check the [TR_PPI GitHub](https://github.com/harzav/TR_PPI_project)<sup>[1]</sup>.

## 5. Acknowledgments

This Tool and Project is directly linked to the [TR_PPI Project](https://github.com/harzav/TR_PPI_project)<sup>[1]</sup>, and therefore is considered part of the VIRTUOUS project (https://www.virtuoush2020.com/), funded by the European Union’s Horizon 2020 research and innovation program under the Marie Sklodowska-Curie-RISE (Grant Agreement No 872181). The tool cannot be copied and distributed without the direct approval of the VIRTUOUS consortium members. Contact [InSyBio PC](https://insybio.com/)<sup>[6]</sup> or this GitHub's authors for more info.

<img src="https://github.com/user-attachments/assets/00ee01f1-4ff0-4ae9-ba58-676fe18eff09" width="297" height="87">
<img src="https://github.com/user-attachments/assets/ae4117a9-404d-49c2-8b16-ddcc531d734c" width="297" height="67">

## 6. Citations

1) TR_PPI Project GitHub: https://github.com/harzav/TR_PPI_project

2) Kabir, M. W. U., & Hoque, M. T. (2024). DisPredict3.0: Prediction of intrinsically disordered regions/proteins using protein language model. Applied Mathematics and Computation, 472, 128630. https://doi.org/10.1016/j.amc.2024.128630

3) Mészáros B, Erdos G, Dosztányi Z. (2018). IUPred2A: context-dependent prediction of protein disorder as a function of redox state and protein binding. Nucleic Acids Res. Jul 2;46(W1):W329-W337. https://doi.org/10.1093/nar/gky384. 

4) Razick, S., Magklaras, G., Donaldson, I. M. (2008). iRefIndex: A Consolidated Protein Interaction Database with Provenance. BMC Bioinformatics. 9 (1), 405. https://doi.org/10.1186/1471-2105-9-405.

5) iRefIndex Documentation Link: https://irefindex.vib.be/wiki/index.php/README_MITAB2.6_for_iRefIndex_19.0#Column_number:_7_.28Method.29

6) For contacting InSyBio: https://insybio.com/ or [info@insybio.com](https://info@insybio.com)
