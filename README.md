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
  - 
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

Example Input
```
Please input melting points (°C), separated by commas: 100, 50
Please input enthalpies of fusion (kJ/mol), separated by commas: 30, 20
Please input clearing points (°C), separated by commas: 120, 60
```
This then produces a plot:<br>
![image](https://github.com/user-attachments/assets/83f170a8-fa54-4751-9fcd-75693a06032a)<br><br>
and data on the composition is returned to the terminal window:<br>
![image](https://github.com/user-attachments/assets/8590faa4-ef2b-44be-b65f-b24881288aec)<br><br>

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

** Functions **
* initialize(): Parses command-line arguments.
* load_compound_data(): Loads compound data from a CSV file.
* generate_combinations(): Generates all possible combinations of compounds.
* evaluate_combinations(): Evaluates each combination to find suitable mixtures.

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
etc
