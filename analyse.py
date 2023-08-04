from analysis import *
import copy
import pprint

experiments = []

lmolecules = [
    "H_ttt",
    "P_TT",
    "B_T",
    "C1",
    "G1",
    "1",
    "99",
    "N4H6_1",
    "H4P2O7",
    # "H2S2O7_1", # crest could not generate a conformer
    "SI5H12_1",
    "ALA_xai",
    "ARG_xak",
    "ASN_xab",
    "ASP_xau",
    "GLN_xai",
    "GLU_xal",
    "GLY_xac",
    "HIS_xau",
    "ILE_xaj",
    "LEU_xad",
    "LYS_xan",
    "PHE_xaw",
    "PRO_xac",
    "THR_xab",
    "TRP_xaf",
    "TYR_xar",
    "VAL_xak",
    "MET_xbm",
    "SER_xak",
    "B11",
    "GLY_b",
    "CYS_xag",
    "tryptophan",
    "aspartic-acid",
    "cysteine",
    "serine",
    "SER_aL",
]


for name in lmolecules:
    exp = Experiment(
        f"/g/data/kx58/jm0124/gResearch/gResults/rdkit_random/{name}",
        name,
        baselines=[
            "setConformerSearch",
            "genetic",
            "crest",
        ],
    )
    experiments.append(exp)

expCollection = ExperimentCollection(experiments, name="Ensemble")


for nConfs in [1, 3, 5]:
    print(f"### {nConfs} ###")

    copied_collection = copy.deepcopy(expCollection)

    copied_collection.buildRelevantStats(
        nConfs=nConfs,
        baselines=["setConformerSearch", "genetic", "crest"],
    )
    pprint.pprint(copied_collection.getTopMetrics(top=nConfs, cutoff=0.00190439916))


nConfs = 1

copied_collection = copy.deepcopy(expCollection)

copied_collection.buildRelevantStats(
    nConfs=nConfs,
    baselines=["setConformerSearch", "genetic", "crest"],
)

copied_collection.graphEnergyBarPlot(
    "graphs"
)