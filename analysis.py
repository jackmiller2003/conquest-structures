import copy
import os
import pprint as pp

import matplotlib.font_manager
import matplotlib.pyplot as plt
import numpy as np


class ExperimentCollection:
    def __init__(self, experiments, name="", directory="graphs"):
        self.experiments = experiments
        self.name = name
        self.directory = directory

    def buildRelevantStats(self, nConfs, baselines):
        self.statsDict = {}
        self.baselines = baselines
        self.nConfs = nConfs

        for experiment in self.experiments:
            self.statsDict[experiment.moleculeName] = {
                "lowestEnergy": {},
                "allEnergies": {},
                "nCands": {},
            }

            experiment.loadResults(nConfs)

            for name, energyList in experiment.energyResults:
                if energyList == []:
                    continue
                self.statsDict[experiment.moleculeName]["lowestEnergy"][name] = min(
                    energyList
                )

            mostMin = min(
                self.statsDict[experiment.moleculeName]["lowestEnergy"].values()
            )

            for name, energyList in experiment.energyResults:
                self.statsDict[experiment.moleculeName]["allEnergies"][name] = [
                    energy - mostMin for energy in energyList
                ]

    def graphEnergyBarPlot(self, path="/scratch/kx58/jm0124/gResearch/graphs"):
        baselines = self.baselines

        dictOfEnergiesAndBaselines = {}

        for baseline in baselines:
            dictOfEnergiesAndBaselines[baseline] = []

        for statData in self.statsDict.values():
            for baseline, data in statData["allEnergies"].items():
                dictOfEnergiesAndBaselines[baseline].extend(data)

        fig, ax = plt.subplots(dpi=300, figsize=(7.5, 5))
        fig.tight_layout()

        bins = [0.00190439916]
        for i in range(1, 3):
            bins.append(round(bins[i - 1] * 5, 5))

        # Create a list of colors for the colormap
        colors = []
        for i, baseline in enumerate(self.baselines):
            if baseline == "setConformerSearch":
                colors.append("black")
            else:
                colors.append(plt.cm.tab20c(4 * i))

        # print(bins)
        width = (1 / (len(baselines))) - 0.05

        x = np.linspace(-5, len(bins) + 10, len(bins) + 16)

        print("\n")
        pp.pprint(dictOfEnergiesAndBaselines)
        print("\n")

        for i, bl in enumerate(self.baselines):
            energies = dictOfEnergiesAndBaselines[bl]

            prev_binned = -10
            label = bl
            if bl == "setConformerSearch":
                label = "Reference (GMTKN55 + BOSS)"
            if bl == "genetic":
                label = "CONQUEST (SRS-RI-MP2/cc-pVDZ)"
            if bl == "crest":
                label = "CREST (GFN2-xTB)"

            c = colors[i]

            for j, binned in enumerate(bins):
                numberInBin = len(
                    [
                        energy
                        for energy in energies
                        if energy < binned and energy >= prev_binned
                    ]
                )
                rects1 = ax.bar(
                    x[j + 5] + (i + 0.5) * (width),
                    int(numberInBin),
                    width,
                    color=c,
                    alpha=0.8,
                    label=label,
                    edgecolor="black",
                )
                prev_binned = binned

        font = {
            "family": "Times New Roman",
            "size": 18,
        }

        matplotlib.rc("font", **font)

        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["bottom"].set_color("black")
        ax.spines["left"].set_color("black")

        ax.set_xticks(x[6:-10])

        # Create a new list of bins in kJ/mol
        nbins = [(r"" + str(round(b * 2625.5))) for b in bins]
        ax.set_xticklabels(nbins, fontsize=18, fontdict=font)

        ax.set_yticklabels(ax.get_yticks(), fontsize=18, fontdict=font)

        ax.set_xlabel("Energy (kJ/mol)", fontsize=18, fontdict=font)
        ax.set_ylabel("Number of conformers", fontsize=18, fontdict=font)

        plt.setp(ax.get_xticklabels(), rotation=30, ha="right")

        # plt.title(f'Bar Chart of Conformers Taking the Top {self.nConfs} Lowest Energy Conformers')
        handles, labels = plt.gca().get_legend_handles_labels()
        by_label = dict(zip(labels, handles))
        plt.legend(by_label.values(), by_label.keys())

        plt.savefig(
            f"{path}/energyDistributionBar-{self.name}-{self.nConfs}.pdf",
            dpi=300,
            bbox_inches="tight",
        )

    def getTopMetrics(self, top, cutoff):
        topMetrics = {}

        for molecule, statData in self.statsDict.items():
            topMetrics[molecule] = {}

            oneBaselineWithoutEnough = False
            for baseline in self.baselines:
                if len(statData["allEnergies"][baseline]) < top:
                    oneBaselineWithoutEnough = True
                    break

            # If one baseline does not have enough structures then skip this molecule
            if oneBaselineWithoutEnough:
                continue

            bestWorstEnergySet = float("inf")

            for baseline, energies in statData["allEnergies"].items():
                sortedEnergies = list(sorted(energies))
                worstEnergy = sortedEnergies[top - 1]
                if worstEnergy < bestWorstEnergySet:
                    bestWorstEnergySet = worstEnergy

            for baseline, energies in statData["allEnergies"].items():
                sortedEnergies = list(sorted(energies))
                worstEnergy = sortedEnergies[top - 1]

                if worstEnergy < cutoff + bestWorstEnergySet:
                    topMetrics[molecule][baseline] = True
                else:
                    topMetrics[molecule][baseline] = False

        successRates = {}
        for baseline in self.baselines:
            numerOfEntries = 0
            successRates[baseline] = 0
            for molecule, metrics in topMetrics.items():
                # print(molecule, metrics)
                if baseline in metrics.keys():
                    numerOfEntries += 1
                    if metrics[baseline] == True:
                        successRates[baseline] += 1

                    if baseline == "setConformerSearch" and metrics[baseline] == False:
                        print(f"{molecule} is not best for setConformerSearch")
                        print(
                            f'Its energies are {self.statsDict[molecule]["allEnergies"][baseline]}'
                        )
                        print(
                            f'Whereas for genetic... {self.statsDict[molecule]["allEnergies"]["genetic"]}'
                        )

            print(
                f"Baseline {baseline} has {successRates[baseline]} successes out of {numerOfEntries} molecules"
            )

            successRates[baseline] = successRates[baseline] / numerOfEntries

        pp.pprint(successRates)


class Experiment:
    """
    Defines an experiment involving one genetic search.
    """

    def __init__(self, directory, moleculeName, baselines):
        self.directory = directory
        self.baselines = baselines
        self.moleculeName = moleculeName
        self.info = {}
        self.geneticRoundInfo = {}

    def getEnergies(self, name):
        structureNames = os.listdir(f"{self.directory}/{name}")
        energy_order_pairs = []

        for structureName in structureNames:
            if structureName.endswith(".sdf"):
                order = int(structureName.split("_")[0])
                num = float(structureName.split("-")[-1][:-4])
                energy_order_pairs.append((-num, order))

        if name == "setConformerSearch" or name == "genetic":
            sorted_energy_order_pairs = sorted(energy_order_pairs, key=lambda x: x[0])
        else:
            sorted_energy_order_pairs = sorted(energy_order_pairs, key=lambda x: x[1])

        return [pair[0] for pair in sorted_energy_order_pairs]

    def loadResults(self, numConfs):
        self.energyResults = []
        for name in self.baselines:
            try:
                energies = self.getEnergies(name)
            except Exception as e:
                print(f"Failed to get energies for {name} with {self.moleculeName}")
                print(f"Error: {e}")
                continue
            self.energyResults.append((name, energies))

        truncatedResults = []

        for name, energies in self.energyResults:
            newEnergies = copy.deepcopy(energies)

            truncatedResults.append((name, newEnergies[:numConfs]))

        self.energyResults = truncatedResults
