# conquest-structures

Structures from CONQUEST letter to JCTC. The file structure for each of the three zip files is explained below. In addition, the exact structures used during benchmarking are provided in the tables.

## File Structure

**_GMTKN55_and_bayesian_**

- (Molecule)
  - crest: The crest structures
  - genetic: The CONQUEST structures found via cc-pVDZ
  - setConformerSearch: The ground-truth set of conformers
- (Molecule)
  - ...
  - ...
- ...

**_ligand_benchmark**

- (Molecule)
  - crest: The crest structures
  - Molecule-(i): The ith run of CONQUEST with the molecule
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
| SER_aL              | :heavy_check_mark:  | |
| H2S2O7_1            | :x:                 | CREST failed to run |
| 2p                  | :x:                 | Reached 4 hour in-node limit on the Gadi supercomputer |
| N3H5_1              | :x:                 | CONQUEST could not identify enough Torsion angles |
| pc22                | :x:                 | Single ring |
| ant                 | :x:                 | No torison angles |
| P7H7_1              | :x:                 | Single ring |
| SI6H12_1            | :x:                 | Single ring |
| octane2             | :x:                 | Too many torsion angles for finite difference optimisation with compute capacity |
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
| DG                  | :heavy_check_mark:  |
| FC3                 | :heavy_check_mark:  |
| I1P                 | :heavy_check_mark:  |

## Included Molecules in 37Conf8 Dataset

| Molecule Name       | Included            | Reasoning |
| :---                | :---:               | ---: |
| cembrene            | :heavy_check_mark:  | |
| cimetidine          | :heavy_check_mark:  | |
| dapiprazole         | :heavy_check_mark:  | |
| eperisone           | :heavy_check_mark:  | |
| fotemustin          | :heavy_check_mark:  | |
| gaba                | :heavy_check_mark:  | |
| glufosinate         | :heavy_check_mark:  | |
| indometacin         | :heavy_check_mark:  | |
| isradipine          | :heavy_check_mark:  | |
| knoevenagel         | :heavy_check_mark:  | |
| liu                 | :heavy_check_mark:  | |
| mannosulfan         | :heavy_check_mark:  | |
| n_iprit             | :heavy_check_mark:  | |
| oseltamivir         | :heavy_check_mark:  | |
| pantoprazole        | :heavy_check_mark:  | |
| propranolol         | :heavy_check_mark:  | |
| prostaglandin       | :heavy_check_mark:  | |
| pyridoxalphosphate  | :heavy_check_mark:  | |
| quinuclidinyl       | :heavy_check_mark:  | |
| rh34                | :heavy_check_mark:  | |
| rolipram            | :heavy_check_mark:  | |
| rosmaridiphenol     | :heavy_check_mark:  | |
| shiepox_depox       | :heavy_check_mark:  | |
| takemoto            | :heavy_check_mark:  | |
| triacetin           | :heavy_check_mark:  | |
| uh232               | :heavy_check_mark:  | |
| VX                  | :heavy_check_mark:  | |
| WAY100635           | :heavy_check_mark:  | |
| xylitolpentanit     | :heavy_check_mark:  | |
| atropine            | :x:                 | RDKIT failed to embed molecule |
| artemizinin         | :x:                 | RDKIT failed to embed molecule |
| 18crown6            | :x:                 | Macrocycle        |
| bispidin            | :x:                 | Macrocycle        |
| cryptand            | :x:                 | Macrocycle        |
| gemcitabine         | :x:                 | RDKIT failed to embed molecule |
| rambo               | :x:                 | RDKIT failed to embed molecule |
| tsogoeva            | :x:                 | RDKIT failed to embed molecule |

## CONF196

| Molecule Name       | Included            | Reasoning |
| :---                | :---:               | ---: |
| COHVAW              | :heavy_check_mark:  | |
| FGG                 | :heavy_check_mark:  | |
| GFA                 | :heavy_check_mark:  | |
| GGf01               | :heavy_check_mark:  | |
| GS464               | :heavy_check_mark:  | |
| GS557               | :heavy_check_mark:  | |
| SANGLI              | :heavy_check_mark:  | |
| WG                  | :heavy_check_mark:  | |
| WGG                 | :heavy_check_mark:  | |
| YIVNOG              | :heavy_check_mark:  | |
| POXTRD              | :x:                 | Macrocycle        |
| CAMVES              | :x:                 | Macrocycle        |
| CHPSAR              | :x:                 | Macrocycle        |
| 18crown6            | :x:                 | Macrocycle        |
