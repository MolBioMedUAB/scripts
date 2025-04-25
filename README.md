# scripts
Collection of useful scripts for the kind of calculations performed in the group


## get_RESP.py

Script for calculating RESP charges of a molecular system. The system can be optimised by the program. As an output, a list of the atoms (with the input order) and the corresponding RESP charge of each atom is delivered.
Source: https://wires.onlinelibrary.wiley.com/doi/full/10.1002/wcms.1457

### TODO
- [ ] Generate mol2 with GAFF atom types and frcmod with GAFF ff using antechamber/parmchk2 and parse RESP charges.
- [ ] Test multiconformer RESP charge fitting
- [ ] Add automatic designation of equivalent Hs
