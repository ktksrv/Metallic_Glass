# Bulk Metallic Glass Generation Scripts
For the three Bulk Metallic Glasses tested with the workflow (monolithic (Cu-Cu), Tertiary (Pd-Ni-P), Multi-Component (Vitreloy-105), the LAMMPS script for generating 32,800 atoms BMG is provided with the potential used ( Neuro-evolution and Embedded-Atom ) and the initial structure.
## Method
1. A 32,800-atom FCC structure with a vacuum at the top and bottom is taken as the initial structure (to impose an artificial barostat with zero target pressure). For multi-component and tertiary glasses, the different atom types are randomly distributed in the FCC lattice to satisfy the proportion requirement of the glass.
2. The initial structure is heated above its melting point, equilibrated, and then cooled fast enough ( between 10<sup>10</sup>-10<sup>13</sup> K/s cooling rate) to produce a metallic glass. 

