# Eutectic
Calculates phase diagram and eutectic composition of n-component mixtures of liquid crystalline (and other) materials. 

- [Features](#features)
- [Requirements](#requirements)
- [Installation](#installation)
- [Usage](#usage)
  - [eutectic.py](#eutecticpy)
  - [mix_scout.py](#mix_scoutpy)
- [Examples](#examples)
- [Background](#background)

## Features
- **Calculate and display phase diagrams** for 2, 3, or 4-component systems.
- **Find eutectic points** and print their compositions.
- **Plot phase diagrams** with automatic handling of different component numbers.
- **Interactive command-line interface** for inputting data.
- **Automated mixture scouting** with `mix_scout.py` to find optimal mixtures from a dataset.
- **Can estimate the eutectic point for a mixture of _any_ number of components

## Requirements
- Python 3.6 or higher
- Libraries:
  - `numpy`
  - `matplotlib`
  - `pandas` (for `mix_scout.py`)
  - `tqdm` (for `mix_scout.py`)
    
## Installation
Clone the repository and install the required libraries:

```
git clone https://github.com/RichardMandle/eutectic-phase-diagram.git
cd eutectic-phase-diagram
pip install -r requirements.txt
```

Alternatively, install the required libraries manually:

```pip install numpy matplotlib pandas tqdm```

## Usage
This script calculates phase diagrams and eutectic points based on user-provided data. You can either run it from ```eutectic.py```, which will prompt you for input (enter values for melting point, enthalpy and clearing-point for each materal, separated by commas)
## eutectic.py
**Running the Script**
Run the script from the command line:
```python eutectic.py```
**Input Data**
The script will prompt you to input the following data:

* Melting points (in °C): Enter the melting points of each component, separated by commas.
* Enthalpies of fusion (in kJ/mol): Enter the enthalpies of fusion for each component, separated by commas.
* Clearing points (in °C): Enter the clearing points (isotropic transition temperatures) for each component, separated by commas.
* Note: Ensure that all inputs have the same number of components.

## Examples
** Using ```eutectic.py```
* Run the script:
```python eutectic.py```
* Input data when prompted:
```
Please input melting points (°C), separated by commas: 100, 80, 60
Please input enthalpies of fusion (kJ/mol), separated by commas: 30, 25, 20
Please input clearing points (°C), separated by commas: 120, 110, 90
```
* View the generated phase diagram and compoisition:
```
Eutectic composition:
14.286% Component A
42.857% Component B
42.857% Component C
Melting point: 52.5 °C
Clearing point: 96.43 °C
Enantiotropic range: 43.93 °C
```
### Example with Plotting
```
Please input melting points (°C), separated by commas: 100, 50
Please input enthalpies of fusion (kJ/mol), separated by commas: 30, 20
Please input clearing points (°C), separated by commas: 120, 60
```
This then produces a plot:<br>
![image](https://github.com/user-attachments/assets/83f170a8-fa54-4751-9fcd-75693a06032a)<br><br>
and data on the composition is returned to the terminal window:<br>
![image](https://github.com/user-attachments/assets/8590faa4-ef2b-44be-b65f-b24881288aec)<br><br>
<br><br>
By entering additional number of components we can visualise 3-component mixtures:<br>
<img width="866" height="666" alt="image" src="https://github.com/user-attachments/assets/cfd9df02-aa60-415c-9020-019cd4c3214b" />
And again, information is returned to the terminal window:
```
--- Eutectic Diagram Generator ---

Enter melting points (°C), comma-separated: 20, 40, 30
Enter enthalpies of fusion (kJ/mol), comma-separated: 20, 20, 10
Enter clearing points (°C), comma-separated: 100, 50, 30

Eutectic composition:
31.549% Component A
18.991% Component B
49.461% Component C
Melting point: -15.68 °C
Clearing point: 55.88 °C
Enantiotropic range: 71.56 °C
```

We can also plot 4-component mixtures; we can choose between different visualisers:<br>
<img width="808" height="267" alt="image" src="https://github.com/user-attachments/assets/ecb54073-e2fb-441d-8da7-71a52dc41428" /><br>
Here we'll use the Mayavi Scatter (#2); the matplotlib one is slow. This gives:
<img width="523" height="435" alt="image" src="https://github.com/user-attachments/assets/69503da1-a043-4e70-aeca-ff300a5760b7" />
<br> you can adjust opacity etc with the normal Mayavi tools.

There are no visualisation options for higher dimensional mixtures.

### What it does Internally###
The code executes a few functions internally:<br>
* make_phase_diagram(): Calculates the phase diagram based on provided data.
* plot_phase_diagram(): Plots the phase diagram for 2 or 3-component systems.
* print_composition(): Prints the eutectic composition and related data.
<br>
## mix_scout.py
This script automates the search for optimal mixtures from a dataset of pure compounds. It evaluates many combinations of compounds to find mixtures with desired properties, such as low melting points

Command-Line Arguments
Run the script with the following options:

```
python mix_scout.py -i compounds.csv -mp 100 -n 2 -s 10000
-i, --input: Path to the CSV file containing compound data.
-mp, --max_melt: Maximum acceptable melting point for mixtures (default: 100°C).
-n, --max_comps: Maximum number of components in a mixture (default: 2).
-s, --search_size: Number of points to evaluate per mixture (default: 10,000).
```

** Input Data File **
Prepare a CSV file (e.g., compounds.csv) with the following format (no header):

```
Compound A, Melting Point A, Clearing Point A, Enthalpy A
Compound B, Melting Point B, Clearing Point B, Enthalpy B
```

Columns:
* Compound Name
* Melting Point (°C)
* Clearing Point (°C)
* Enthalpy of Fusion (kJ/mol)

Example CSV Content:
```
Compound1, 100, 120, 30
Compound2, 50, 60, 20
Compound3, 80, 110, 25
```
** Output **
The script will evaluate all possible combinations up to the specified maximum number of components and output the mixtures that meet the criteria. It prints:
* The lowest melting point mixture.
* The highest clearing point mixture.
* The mixture with the widest enantiotropic range.

### Example:
First, we'll use data published in this paper (https://pubs.acs.org/doi/10.1021/jacs.4c16555). We store this as a .csv file:

```
RM2931,105.9,106,38.1
RM2930,102.5,82,25.3
RM2929,128.2,145.4,25.6
CJG106,108,98.7,36
CJG104,158.7,133.5,38.8
CJG97,122.4,95.4,30.9
CJG98,149.2,130.6,35.5
```

We feed this .csv file into mix_scout.py. We'll also specify some other options, namely ```-mp 80``` to set the maximum melting point to 80 °C, ```-s 1000000``` to increase the number of compositons evaluated to a million, and ```-n 5``` to look at up to 5-component mixtures:<br>```python mix_scout.py -i jacs_cpd_data.csv -mp 80 -s 10000000 -n 5```. 

In the terminal, we see the output:<br>
<img width="628" height="590" alt="image" src="https://github.com/user-attachments/assets/611516b6-67f6-4370-a7d9-b8ebfec6446e" />
There are many possible mixtures, but here we just get the lowest melting point, the highest clearing point and the widest range one. The data can be saved to CSV with the ```-export``` flag, just pass a filename to save as.

