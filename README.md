Methanol synthesis model developed by Kristian Tiiro from University of Oulu.
Funded by Business Finland via project ‘Bio-CCU - Creating sustainable value of the bio-based CO2’ (2352/31/2022).

Model description

Model portrays a methanol synthesis plant including:
 - Lurgi-type reactor with Cu/ZnO/Al2O3 catalyst
 - Knockout drum for separating crude methanol stream from reactor outlet
 - Recycle stream

The model can be simulated from main.m
The simulation is able to model the steady state behaviour of the plant 
at different catalyst activation values, and age the catalyst dynamically.

Mass and energy balances are considered, and the most important outputs 
from the model are:
 - Product stream flow
 - Yields, conversions and selectivities
 - Specific electricity consumption
 - Specific net heat consumption (assuming all heat streams are 100 % utilizable)

The literature and assumptions used for the model are displayed extensively 
in the comments of the code.
