# Eutectic
Calculates phase diagram and eutectic composition of n-component mixtures of liquid crystalline materials. 

Takes melting points (degrees C), enthalpies (kJ mol^-1) and some other phase transition temperature (e.g. N-Iso, SmA-N_F). Generates a large number of possible concentrations and evaluates the melting point and clearing point of each _via_ the method described by E.P.Raynes _et al_ in DOI: 10.1039/C39740000098

Plotting: with 2-components gives you get a simple scatter diagram; with 3-components you get a triangular ternary colourmap plot; with 4-components you get a sort of 3D pyramidal scatter plot with colour and alpha mapping. With 5+ components ye get nowt.

# Issues
Reliable plotting limited to 2- or 3- component systems; visualising 4 component systems is sub-optimal, 5-components and above is wild.

# Future Plans
Ability to plot a 'slice' through 3+ component data by fixing one or more concentrations; reduces the dimensionality. 

# Usage
Give it some data:
_2-components_
Here we'll use data for RM734 and RM734-F (10.1002/chem.201702742) and a few other favorites to find the eutectic composition and a enantiotropic NF material:

~~~
import Eutectic 
Melts = [139.8, 165.24] # melting points, in degrees C
Enthalp = [29.89,50.06] # enthalpy of fusion, in kJ mol-1
Isos = [132.7,139.6]    # Nf-N temperature in degrees C

mix_concs,mix_melt,mix_transition = Eutectic.DoPhaseDiagram(Melts,Nf2N,Enthalp) # calculate phase diagram

~~~
![image](https://github.com/RichardMandle/Eutectic/assets/101199234/21a105c3-ef73-4047-9860-0fdfcc4905d6)
~~~
Eutectic.PrintComposition(mix_concs,mix_melt,mix_transition)
~~~
~~~
eutectic composition:
0.741 mol% A
0.259 mol% B
Melting point: [126.03] °C
Clearing point: 134.49 °C
enantiotropic range = 8.46 °C
~~~

_3-components_
~~~
import Eutectic 
Melts = [139.8, 139.0, 165.24] # melting points, in degrees C
Enthalp = [29.89,34.79,50.06]  # enthalpy of fusion, in kJ mol-1
Isos = [132.7,86.5,139.6]      # Isotropisation temperature in degrees C (set to -273.15 for non LC)

mix_concs,mix_melt,mix_transition = Eutectic.DoPhaseDiagram(Melts,Nf2N,Enthalp) # calculate phase diagram

~~~
![image](https://github.com/RichardMandle/Eutectic/assets/101199234/db939071-01c3-4b63-82e8-6a5c7f34e943)
~~~
Eutectic.PrintComposition(Concs,Melts,Clear)                   # print eutectic composition
~~~
~~~
eutectic composition:
0.464 mol% A
0.418 mol% B
0.119 mol% C
Melting point: 106.35 °C
Clearing point: 114.23 °C
enantiotropic range = 7.88 °C
~~~

_4-components_
~~~
import Eutectic 
Melts = [139.8, 62.2, 37.2, 99.7]
Enthalp = [29.89,24.12,20.6, 27.5]
Nf2N = [132.7,78.5,15.4, 68]

mix_concs,mix_melt,mix_transition = Eutectic.DoPhaseDiagram(Melts,Nf2N,Enthalp)
~~~
![image](https://github.com/RichardMandle/Eutectic/assets/101199234/7ba94d00-5331-4471-bdda-6b1f773d794b)
~~~
Eutectic.PrintComposition(mix_concs,mix_melt,mix_transition)
~~~
~~~
eutectic composition:
0.028 mol% A
0.278 mol% B
0.614 mol% C
0.08 mol% D
Melting point: [19.35] °C
Clearing point: 40.41 °C
enantiotropic range = 21.06 °C
~~~

_>5-components_

~~~
import Eutectic 
Melts =   [139.8, 139.0, 155.2, 164.7, 138.8, 143.3, 143.2, 141.8, 165.24]
Enthalp = [29.89, 34.79, 33.12, 28.79, 38.75, 36.05, 40.61, 43.66, 50.06]
Nf2N = [132.7, 86.5, 42.6, 66.1, 32.7, 91, 77.6, 117.1, 139.56 ]

mix_concs,mix_melt,mix_transition = Eutectic.DoPhaseDiagram(Melts,Nf2N,Enthalp)
~~~
~~~
Eutectic.PrintComposition(mix_concs,mix_melt,mix_transition)

Not plotting; 9 components
eutectic composition:
0.184 mol% A
0.161 mol% B
0.121 mol% C
0.107 mol% D
0.09 mol% E
0.129 mol% F
0.092 mol% G
0.088 mol% H
0.028 mol% I
Melting point: [76.23] °C
Clearing point: 86.64 °C
enantiotropic range = 10.4 °C
~~~

# Data Export
We can write the data to .csv for later retrieval like this:
~~~
Eutectic.SaveData(mix_concs,mix_melt,mix_transition,'MyEutecticData.csv')
~~~
