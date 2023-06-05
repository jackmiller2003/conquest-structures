# conquest-structures
Structures from CONQUEST letter to JCTC. The file structure for each of the three zip files is explained below. In addition, the exact structures used during benchmarking are provided in the tables.

## File Structure

**_GMTKN55_and_bayesian_**
  - (Molecule)
    - crest: The crest structures
    - CONQUEST-CC-new: The CONQUEST structures found via cc-pVDZ
    - CONQUEST-6-31G-new: The CONQUEST structures found via 6-31G*
    - CONUQEST-STO-old: The CONQUEST structures found via STO-3G with an old selection mechanism
    - setConformerSearch: The ground-truth set of conformers
  - (Molecule)
    - ...
    - ...
  - ...

**_ligands_**
  - (Molecule)
    - crest: The crest structures
    - Molecule-(i): The ith run of CONQUEST with the molecule
  - ...

**_rna_**
  - (Molecule)
    - crest: The crest structures
    - Molecule-ngpus4-(i): The ith run of CONQUEST with the molecule with 4 GPUs
  - ...

## Included Molecules in GMTKN55 Dataset

| Molecule Name       | Included            | Reasoning |
| :---                | :---:               | ---: |
| H_ttt               | :heavy_check_mark:  | |
| P_TT                | :heavy_check_mark:  | |
| B_T                 | :heavy_check_mark:  | |
| C1                  | :heavy_check_mark:  | |
| G1                  | :heavy_check_mark:  | |
| 1                   | :heavy_check_mark:  | |
| 99                  | :heavy_check_mark:  | |
| N4H6_1              | :heavy_check_mark:  | |
| H4P2O7              | :heavy_check_mark:  | |
| H2S2O7_1            | :heavy_check_mark:  | |
| SI5H12_1            | :heavy_check_mark:  | |
| ALA_xai             | :heavy_check_mark:  | |
| ARG_xak             | :heavy_check_mark:  | |
| ASN_xab             | :heavy_check_mark:  | |
| ASP_xau             | :heavy_check_mark:  | |
| GLN_xai             | :heavy_check_mark:  | |
| GLU_xal             | :heavy_check_mark:  | |
| GLY_xac             | :heavy_check_mark:  | |
| HIS_xau             | :heavy_check_mark:  | |
| ILE_xaj             | :heavy_check_mark:  | |
| LEU_xad             | :heavy_check_mark:  | |
| LYS_xan             | :heavy_check_mark:  | |
| PHE_xaw             | :heavy_check_mark:  | |
| PRO_xac             | :heavy_check_mark:  | |
| THR_xab             | :heavy_check_mark:  | |
| TRP_xaf             | :heavy_check_mark:  | |
| TYR_xar             | :heavy_check_mark:  | |
| VAL_xak             | :heavy_check_mark:  | |
| THR_xab             | :heavy_check_mark:  | |
| TRP_xaf             | :heavy_check_mark:  | |
| TYR_xar             | :heavy_check_mark:  | |
| VAL_xak             | :heavy_check_mark:  | |
| tryptophan          | :heavy_check_mark:  | |
| cysteine            | :heavy_check_mark:  | |
| serine              | :heavy_check_mark:  | |
| aspartic-acid       | :heavy_check_mark:  | |
| MET_xbm             | :heavy_check_mark:  | |
| SER_xak             | :heavy_check_mark:  | |
| B11                 | :heavy_check_mark:  | |
| GLY_b               | :heavy_check_mark:  | |
| CYS_xag             | :heavy_check_mark:  | |
| N3H5_1              | :x:                 | CONQUEST could not identify enough Torsion angles |
| SER_aL              | :x:                 | Beyond 40 atom cap |
| pc22                | :x:                 | Single ring |
| ant                 | :x:                 | No torison angles |
| P7H7_1              | :x:                 | Single ring |
| SI6H12_1            | :x:                 | Single ring |
| octane2             | :x:                 | Too many torsion angles for finite difference optimisation with compute capacity |
| 2p                  | :x:                 | Beyond 40 atom cap |
| F14f                | :x:                 | Too many torsion angles for finite difference optimisation with compute capacity |
| S8_1                | :x:                 | Single ring |
| N3P3H12_1           | :x:                 | Single ring |
| antdimer            | :x:                 | No torsion angles |
| octane1             | :x:                 | No torsion angles |
| octane2             | :x:                 | Too many torsion angles for finite difference optimisation with compute capacity |
| F22f                | :x:                 | Too many torsion angles for finite difference optimisation with compute capacity |

## Included Structures from OMEGA Validation Set

| Molecule Name       | Included            |
| :---                | :---:               |
| PEZ                 | :heavy_check_mark:  |
| ANP                 | :heavy_check_mark:  |
| L4G                 | :heavy_check_mark:  |

## Random RNA Sequences

| Molecule Name       | Included            |
| :---                | :---:               |
| GUA                 | :heavy_check_mark:  |
| UCC                 | :heavy_check_mark:  |
| GUC                 | :heavy_check_mark:  |
| GAC                 | :heavy_check_mark:  |
