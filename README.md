# Eutectic
Calculates phase diagram and eutectic composition of n-component mixtures of liquid crystalline materials.

# Usage
Give it some data:
~~~
import Eutectic 
Melts = [139.8, 139.0, 165.24] # melting points, in degrees C
Enthalp = [29.89,34.79,50.06]  # enthalpy of fusion, in kJ mol-1
Isos = [132.7,86.5,139.6]      # Isotropisation temperature in degrees C (set to -273.15 for non LC)

Concs,Melt,Clear = Eutectic.DoPhaseDiagram(Melts,Isos,Enthalp) # calculate phase diagram

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

