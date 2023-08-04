# --- Internal Imports --- #
# from src.molecule import *

# --- External Imports --- #
import scipy.stats
import os
import json
import pprint as pp
import matplotlib.pyplot as plt
import numpy as np
import statistics as stats
from collections import OrderedDict
from scipy.optimize import curve_fit
import matplotlib.patheffects as pe
from scipy.optimize import minimize
import matplotlib.font_manager
import copy

class ExperimentCollection():
    """
    Defines a collection of experiments.
    """

    def __init__(self, experiments, name='', directory='/scratch/kx58/jm0124/gResearch/graphs'):
        self.experiments = experiments
        self.name = name
        self.directory = directory

        # for i in range(0, len(directories)):
        #     self.experiments.append(Experiment(directories[i], names[i]))
        
    def getAllResults(self):
        for experiment in self.experiments:
            experiment.loadResults()
    
    def getBiggestDifference(self, baseline1, baseline2):
        """
        Having built relevant statistics, finds the biggest difference between two baselines.
        """

        biggestDiff = 0
        biggestDiffMolecule = ''

        for molecule, statData in self.statsDict.items():
            diff = abs(statData['lowestEnergy'][baseline1] - statData['lowestEnergy'][baseline2])
            if diff > biggestDiff:
                biggestDiff = diff
                biggestDiffMolecule = molecule
    
        return biggestDiff, biggestDiffMolecule

    def buildRelevantStats(self, nConfs, baselines):
        """
        Builds relevant stats for later graphing
        """

        self.statsDict = {}
        self.baselines = baselines
        self.nConfs = nConfs

        for l, experiment in enumerate(self.experiments):
            self.statsDict[experiment.moleculeName] = {
                'lowestEnergy':{},
                'allEnergies':{},
                'nCands':{},
            }

            experiment.loadResults(nConfs)

            for (name, energyList) in experiment.energyResults:
                if energyList == []: continue
                self.statsDict[experiment.moleculeName]['lowestEnergy'][name] = min(energyList)
            
            mostMin = min(self.statsDict[experiment.moleculeName]['lowestEnergy'].values())

            for (name, energyList) in experiment.energyResults:
                self.statsDict[experiment.moleculeName]['nCands'][name] = len(energyList)
                self.statsDict[experiment.moleculeName]['allEnergies'][name] = [energy - mostMin for energy in energyList]

                
    def graphEnergyBarPlot(self, path='/scratch/kx58/jm0124/gResearch/graphs'):

        baselines = self.baselines

        dictOfEnergiesAndBaselines = {}

        for baseline in baselines:
            dictOfEnergiesAndBaselines[baseline] = []

        for molecule, statData in self.statsDict.items():
            for baseline, data in statData['allEnergies'].items():
                dictOfEnergiesAndBaselines[baseline].extend(data)
            

        fig, ax = plt.subplots(dpi=300, figsize=(7.5, 5))
        fig.tight_layout()

        bins = [0.00190439916]
        for i in range(1, 3):
            bins.append(round(bins[i-1]*5,5))

        # Create a list of colors for the colormap
        colors = []
        for i, baseline in enumerate(self.baselines):
            if baseline == 'setConformerSearch': colors.append('black')
            else: colors.append(plt.cm.tab20c(4*i))
        
        # print(bins)
        width = (1/(len(baselines))) - 0.05

        x = np.linspace(-5,len(bins)+10, len(bins)+16)

        print("dict of energies and baselines")

        pp.pprint(dictOfEnergiesAndBaselines)


        for i, bl in enumerate(self.baselines):
            energies = dictOfEnergiesAndBaselines[bl]

            prev_binned = -10
            label = bl
            if bl == 'setConformerSearch': label = 'Reference'
            if bl == 'genetic': label = 'CONQUEST (CC-pVDZ)'
            if bl == 'genetic-6-31G': label = 'CONQUEST (6-31G*)'
            if bl == 'crest-seeded': label = 'CREST SEEDED (GFN-2)'
            if bl == 'crest': label = 'CREST (GFN-2)'

            """
            if bl == 'genetic-CC-new': c = plt.cm.tab20c(4*2)
            elif bl == 'crest': c = plt.cm.tab20c(4*3)
            else: 
            """

            c = colors[i]

            for j, binned in enumerate(bins):
                numberInBin = len([energy for energy in energies if energy < binned and energy >= prev_binned])
                rects1 = ax.bar(x[j+5] + (i+0.5)*(width), int(numberInBin), width, color=c, alpha=0.8, label=label, edgecolor='black')
                prev_binned = binned
        
        font = {'family' : 'Times New Roman',
        # 'weight' : 'bold',
        'size'   : 18}

        matplotlib.rc('font', **font)
        # Set the axis labels adn everything to hav ethe same font as well
        

        # Axis formatting.
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_color('black')
        ax.spines['left'].set_color('black')

        # Set the xticks to be the bin edges
        ax.set_xticks(x[6:-10])
        # Create a new list of bins in kJ/mol 
        nbins = [(r'' + str(round(b*2625.5))) for b in bins]
        ax.set_xticklabels(nbins, fontsize=18, fontdict=font)
        # Set font of yticks and round them to integers
        ax.set_yticklabels(ax.get_yticks(), fontsize=18, fontdict=font)
        
        ax.set_xlabel('Energy (kJ/mol)', fontsize=18, fontdict=font)
        ax.set_ylabel('Number of conformers', fontsize=18, fontdict=font)

        plt.setp(ax.get_xticklabels(), rotation=30, ha='right')

        # plt.title(f'Bar Chart of Conformers Taking the Top {self.nConfs} Lowest Energy Conformers')
        handles, labels = plt.gca().get_legend_handles_labels()
        by_label = dict(zip(labels, handles))
        plt.legend(by_label.values(), by_label.keys())

        plt.savefig(f'{path}/energyDistributionBar-{self.name}-{self.nConfs}.pdf', dpi=300, bbox_inches='tight')
        plt.show()

    def getTopMetrics(self,top=1,cutoff=1.5e-3):
        """
        Returns a dictionary of the top metrics for each molecule.
        """

        # Iterate through all molecules

        topMetrics = {}

        for molecule, statData in self.statsDict.items():
            # Iterate through all baselines
            topMetrics[molecule] = {}
            # If any baseline does not have enough conformers, then do not consider it in the top metrics
            
            oneBaselineWithoutEnough = False
            for baseline in self.baselines:
                if len(statData['allEnergies'][baseline]) < top:
                    oneBaselineWithoutEnough = True
                    break
            
            if oneBaselineWithoutEnough: continue

            # Find the best worst performer in the energy list
            bestWorstEnergySet = 100

            for baseline, energies in statData['allEnergies'].items():
                sortedEnergies = list(sorted(energies))
                worstEnergy = sortedEnergies[top-1]
                if worstEnergy < bestWorstEnergySet:
                    bestWorstEnergySet = worstEnergy
            
            # print(f'Best worst energy set for {molecule} is {bestWorstEnergySet} with cutoff {cutoff} and top {top}')

            for baseline, energies in statData['allEnergies'].items():
                sortedEnergies = list(sorted(energies))
                worstEnergy = sortedEnergies[top-1]
                
                if worstEnergy < cutoff+bestWorstEnergySet:
                    topMetrics[molecule][baseline] = True
                else:
                    topMetrics[molecule][baseline] = False
        
        # For each baseline get the number of successes and divide by the total number of molecules
        successRates = {}
        for baseline in self.baselines:
            numerOfEntries =  0
            successRates[baseline] = 0
            for molecule, metrics in topMetrics.items():
                # print(molecule, metrics)
                if baseline in metrics.keys():
                    numerOfEntries += 1
                    if metrics[baseline] == True:
                        successRates[baseline] += 1
                    
                    if baseline == 'setConformerSearch' and metrics[baseline] == False: 
                        print(f'{molecule} is not best for setConformerSearch')
                        print(f'Its energies are {self.statsDict[molecule]["allEnergies"][baseline]}')
                        print(f'whereas for genetic... {self.statsDict[molecule]["allEnergies"]["genetic"]}')
            
            print(f'Baseline {baseline} has {successRates[baseline]} successes out of {numerOfEntries} molecules')
            successRates[baseline] = successRates[baseline]/numerOfEntries
           
        
        print(successRates)

    def graphEnergyDistribution(self, path='/scratch/kx58/jm0124/gResearch/graphs', distribution='gaussian', top=1e4, left=1e-4, right=1e-1):
        """
        Graphs energy distribution for each molecule.
        """
        

        baselines = self.baselines

        dictOfEnergiesAndBaselines = {}

        for baseline in baselines:
            dictOfEnergiesAndBaselines[baseline] = []

        for molecule, statData in self.statsDict.items():
            for baseline, data in statData['allEnergies'].items():
                dictOfEnergiesAndBaselines[baseline].extend(data)
        
        # Bar chart of energy distribution

        # mpl.font_manager._rebuild()
        fig, ax = plt.subplots(dpi=350, figsize=(9, 5.5))
        fig.tight_layout()
 
        # Create a list of colors for the colormap
        colors = []
        for i, baseline in enumerate(self.baselines):
            if baseline == 'setConformerSearch': colors.append('black')
            else: colors.append(plt.cm.tab20c(4*i))

        forCurveFitting = []

        distParams = {}
        try: 
            setEnergies = dictOfEnergiesAndBaselines['setConformerSearch']
            setEnergiesNew = []

            for energy in setEnergies:
                if energy > 0:
                    setEnergiesNew.append(energy)
                else: setEnergiesNew.append(1e-6)
            
            setEnergies = setEnergiesNew
        except:
            # print('No set conformer search energies found.')
            pass

        for i, bl in enumerate(self.baselines):
            energies = dictOfEnergiesAndBaselines[bl]
            # print(f'Baseline: {bl} has {energies} energies')

            initParams = [1.1, 0.05]
            energiesNew = []
            for energy in energies:
                if energy > 0:
                    energiesNew.append(energy)
                else: energiesNew.append(1e-6)
            energies = energiesNew
            # print(energies)

            
            if bl == 'setConformerSearch':
                for testDist in ['gaussian', 'gamma', 'chi2']:
                    if testDist == 'gaussian': 
                        results = minimize(self.gaussian, initParams, energies, method='Nelder-Mead')
                        nll = self.gaussian(results.x, energies)
                    elif testDist == 'gamma': 
                        results = scipy.stats.gamma.fit(energies, floc=0)
                        nll = self.gamma(results, energies)
                    elif testDist == 'chi2': 
                        results = scipy.stats.chi2.fit(energies, floc=0)
                        nll = self.chi2(results, energies)
                
                    # print(f'NLL for {testDist} is {nll}')
                    
            gammaResults = scipy.stats.gamma.fit(energies, floc=0)
            try:
                nll = self.gamma(gammaResults, setEnergies)
            except:
                nll = 0
            # print(f'Gamma NLL for {bl} is {nll}')
                    

            if distribution == 'gaussian': results = minimize(self.gaussian, initParams, energies, method='Nelder-Mead')
            elif distribution == 'gamma': results = scipy.stats.gamma.fit(energies, floc=0)
            elif distribution == 'chi2': results = scipy.stats.chi2.fit(energies, floc=0)

            if distribution == 'gaussian': distParams[bl] = results.x
            elif distribution == 'gamma': distParams[bl] = results
            elif distribution == 'chi2': distParams[bl] = results

            fx = np.logspace(-8, -1.5, 300)
            if distribution == 'gaussian': y1 = self.guassFunc(fx, distParams[bl][0], distParams[bl][1])
            elif distribution == 'gamma': y1 = self.gammaFunc(fx, distParams[bl][0], distParams[bl][1], distParams[bl][2])
            elif distribution == 'chi2': y1 = self.chi2Func(fx, distParams[bl][0])
            # print(f'Dist params for {bl}: {distParams[bl]}')
            # print(y1[:10])
            label = bl
            if bl == 'setConformerSearch': label = 'reference dataset'
            if bl == 'genetic': label = 'genetic-cc-pVDZ'
            if bl == 'crest': label = 'CREST'
            ax.plot(fx,y1, label=label, color=colors[i], linewidth=2)

        ax.set_xscale('log')
        # ax.set_yscale('log')
        ax.set_ylim(0, top)
        ax.set_xlim(left, right)
        ax.set_xlabel('Energy (Ha)')
        ax.set_ylabel('Probability Density')

        plt.setp(ax.get_xticklabels(), rotation=45, ha='right')

        plt.title(f'Distribution of Conformers Taking the Top {self.nConfs} Lowest Energy Conformers')

        plt.legend()
        plt.savefig(f'{path}/energyDistribution-{self.name}-{self.nConfs}-{left}-{right}.pdf', dpi=300,  bbox_inches='tight')
        plt.show()

    # --- Different distributions --- #

    def guassFunc(self,sample_data,mean,sd):
        return scipy.stats.norm.pdf(sample_data, loc=mean, scale=sd)

    def gaussian(self, params, sample_data):
        mean = params[0]   
        sd = params[1]

        # Calculate negative log likelihood
        nll = -np.sum(scipy.stats.norm.logpdf(sample_data, loc=mean, scale=sd))

        return nll

    def gamma(self, params, sample_data):
        shape = params[0]   
        scale = params[1]

        # Calculate negative log likelihood
        nll = -np.sum(scipy.stats.gamma.logpdf(sample_data, shape, scale))

        return nll

    def gammaFunc(self,sample_data,shape,loc,scale):
        return scipy.stats.gamma.pdf(sample_data, a=shape, loc=loc, scale=scale)

    def chi2(self, params, sample_data):
        df = params[0]   

        # Calculate negative log likelihood
        nll = -np.sum(scipy.stats.chi2.logpdf(sample_data, df))

        return nll
    
    def chi2Func(self,sample_data,df):
        return scipy.stats.chi2.pdf(sample_data, df)


    def graphAllResults(self, nConfs=3):
        """
        Graphs all results of the experiments
        using a similar method as presented in that class
        but presenting them all on a large plot with subplots
        """

        # mpl.font_manager._rebuild()
        fig, ax = plt.subplots(dpi=250, figsize=(12, 5))
        fig.tight_layout()

        avgDict = {
            'lowestEnergy':[],
            'averageEnergy':[],
            'std':[],
            'nCands':[],
            'timing':[]
        }

        allEnergiesAndBaseline = {}

        for l, experiment in enumerate(self.experiments):

            experiment.loadResults(nConfs)
            lowestEnergy, nCands, averageEnergy, std, names = [], [], [], [], []
            timing = experiment.timing
            timing = [t.second + 60*t.minute + 3600*t.hour for t in timing]

            for (name, energyList) in experiment.energyResults:
                if energyList == []: continue
                lowestEnergy.append(min(energyList))
            
            mostMin = min(lowestEnergy)

            avgDict['lowestEnergy'].append(lowestEnergy)

            for (name, energyList) in experiment.energyResults:
                # if energyList == []: continue
                if not name in allEnergiesAndBaseline.keys():
                    allEnergiesAndBaseline[name] = []

                allEnergiesAndBaseline[name].extend([energy - mostMin for energy in energyList])

                nCands.append(len(energyList))
                averageEnergy.append((sum(energyList) / len(energyList)) - mostMin)
                
                if len(energyList) > 1:
                    std.append(stats.stdev(energyList)/(float(3)))
                else: std.append(0)

                names.append(name)

            if l == 0: oNames = names

            avgDict['averageEnergy'].append(averageEnergy)
            avgDict['std'].append(std)
            avgDict['nCands'].append(nCands)
            avgDict['timing'].append(timing)

            if names != oNames:
                raise Exception('Name order changes')

        # print(allEnergiesAndBaseline)

        for item, infoList in avgDict.items():
            avgDict[item] = [sum(x) / len(x) for x in zip(*infoList)]

        # print(f'Average dict is: {avgDict}')

        self.baselines = names
        averageEnergy = avgDict['averageEnergy']
        std = avgDict['std']
        nCands = avgDict['nCands']
        lowestEnergy = avgDict['lowestEnergy']
        timing = avgDict['timing']
        self.baselines, averageEnergy, std, nCands, lowestEnergy, timing = zip(*sorted(zip(self.baselines, averageEnergy, std, nCands, lowestEnergy, timing), key=lambda x: x[1]))


        tlabels = []

        for i, baseline in enumerate(self.baselines):
            minutes = int(timing[i]/60)
            # print(nCands[i])
            if baseline == 'setConformerSearch': baseline = 'provided'
            # else: labelName = nCands[i]
            tlabels.append(baseline + f" ({nCands[i]}, {minutes})")
        
        # print(tlabels)

        bars = ax.bar(
            x=np.arange(len(averageEnergy)),
            height=averageEnergy,
            tick_label=tlabels,
            yerr=std,
            alpha=0.6,
            capsize=5
        )
        
        fbar = bars[0].get_height()

        for i, bar in enumerate(bars):
            break
            # Write at bottom of every bar the value of the bar
            ax.text(
                bar.get_x() + bar.get_width() / 2,
                bar.get_height()  + fbar*0.02,
                f" {len(self.energyResults[i][1])}",
                ha='left',
                color='black',
                fontsize=15
            )
            # Change label of bar


        # Axis formatting.
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.spines['bottom'].set_color('#DDDDDD')
        ax.tick_params(bottom=False, left=False)
        ax.set_axisbelow(True)
        ax.yaxis.grid(True, color='#EEEEEE')
        # Set yaxis to be log
        ax.set_yscale('log')
        ax.xaxis.grid(False)

        # ax2.set_yticks(ax.get_yticks())
        # ax.set_yscale('log')
        # ax.set_ylim(ax.get_ylim())
        # ax.set_axisbelow(True)
        # print(ax.get_ylim())
        x = np.arange(len(averageEnergy))
        # print(x)
        y = [energy-mostMin + ax.get_ylim()[0]/5 for energy in lowestEnergy]
        # Set bottom of y lim to be the negative of the current bottom
        # ax.set_ylim(bottom=-ax.get_ylim()[0])

        # print(y)
        ax.scatter(x, y, color='red', marker='o', zorder=10)

        # Add text annotations to the top of the bars.
        bar_color = bars[0].get_facecolor()
        fbar = 0

        # Add labels and a title.
        ax.set_xlabel('Method of Computation', labelpad=15, color='#333333')
        ax.set_ylabel('Difference to Minimum', labelpad=15, color='#333333')
        ax.set_title(f'Comparison of Conformer Search Methods for {self.name}', pad=15, color='#333333',
                    weight='bold')

        if not os.path.exists(f'{self.directory}/{self.name}'):
            os.mkdir(f'{self.directory}/{self.name}')

        plt.savefig(f'{self.directory}/{self.name}/average-results.pdf', dpi=250)

        plt.show()

        # mpl.font_manager._rebuild()
        fig, ax = plt.subplots(dpi=250, figsize=(12, 5))
        fig.tight_layout()

        # Create a bar chart using allEnergiesAndBaseline by binning each energy
        # and coloring according to the baseline using a colormap

        # Create a list increasing expoentially from 1e-4 to 1e1
        bins = [4.5e-4]
        for i in range(1, 14):
            bins.append(round(bins[i-1]*1.7,5))

        # Create a list of colors for the colormap
        colors = []
        for i in range(len(self.baselines)):
            colors.append(plt.cm.tab20(i))
        
        # print(bins)
        width = 0.2

        x = np.arange(len(bins))

        forCurveFitting = []

        for i, bl in enumerate(self.baselines):
            energies = allEnergiesAndBaseline[bl]
            # print(energies)
            # Create a histogram of the energies with grouped bins
            prev_binned = 0
            numbers = []
            for j, binned in enumerate(bins):
                numberInBin = len([energy for energy in energies if energy < binned and energy >= prev_binned])
                numbers.append(numberInBin)
                rects1 = ax.bar(x[j] + i*(width), numberInBin, width, color=colors[i], alpha=0.4)
                prev_binned = binned
            
            forCurveFitting.append((bl,numbers))
            # n, bins, patches = ax.hist(energies, bins=bins, histtype='bar', color=colors[i], label=bl, alpha=0.6)

            # Set the color of the bins to be the color of the baseline
            # for i in range(len(patches)):
            #     patches[i].set_facecolor(colors[self.baselines.index(bl)])

        # Axis formatting.
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.spines['bottom'].set_color('#DDDDDD')

        # Set the xticks to be the bin edges
        
        ax.set_xticks(x)
        ax.set_xticklabels(bins)
        # ax.ticklabel_format(axis='x', style='sci')
        
        # ax.set_xlim(left=-0.5, right=len(bins)-0.5)
        plt.setp(ax.get_xticklabels(), rotation=45, ha='right')
        ax.set_ylim(bottom=0)
        # Set xlim to be the number of baselines
        ax.set_xlim(left=-2, right=len(bins)+5)

        ax.tick_params(bottom=False, left=False)
        ax.set_axisbelow(True)

        # Add labels and a title.
        ax.set_xlabel('Energy (Ha)', labelpad=15, color='#333333')
        ax.set_ylabel('Number of Conformers', labelpad=15, color='#333333')
        ax.set_title(f'Comparison of Conformer Search Methods for {self.name}', pad=15, color='#333333')

        # ax2 = ax.twinx()

        # print(forCurveFitting)

        for (bl, numbers) in forCurveFitting:
            # print(f'fitting {bl}, {numbers}')
            # print(f'fitting {x}, {numbers}')
            # Initial nonzero number from numbers
            nonzero = [i for i,n in enumerate(numbers) if n != 0][0]
            # Last nonzero number from numbers
            lastNonzero = [i for i,n in enumerate(numbers) if n != 0][-1]
            # if lastNonzero - nonzero < 3:
            #     continue
            xtoFit = [-5, -4, -3, -2, -1]
            xtoFit.extend(x)
            ytoFit = [0,0,0,0,0]
            ytoFit.extend(numbers)
            ytoFit.extend([0,0,0,0,0])
            xtoFit.extend([i for i in range(max(xtoFit), max(xtoFit)+5)])
            if bl == 'conformator': 
                print(f'{xtoFit} for {ytoFit} with {bl}')
                popt, pcov = curve_fit(self.gaussForFitting, xtoFit, ytoFit, p0=[5, 14, 1])
            else:
                popt, pcov = curve_fit(self.gaussForFitting, xtoFit, ytoFit)
            
            x1 = np.linspace(min(xtoFit), max(xtoFit), 1000)
            y = self.gaussForFitting(x1, *popt)
            if bl == 'conformator': 
                print(popt)
            ax.plot(x1, y, label=bl, color=colors[self.baselines.index(bl)], alpha=1, path_effects=[pe.Stroke(linewidth=3)])

        # Add a legend
        ax.legend()

        plt.savefig(f'{self.directory}/{self.name}/histogram-results.pdf', dpi=250)

        handles, labels = plt.gca().get_legend_handles_labels()
        by_label = OrderedDict(zip(labels, handles))
        plt.legend(by_label.values(), by_label.keys())

        plt.show()
        


class Experiment():
    """
    Defines an experiment involving one genetic search.
    """

    def __init__(self, directory, moleculeName, baselines=['rdkit', 'crest', 'balloon', 'genetic', 'genetic-rdkit', 'genetic-rdkit-nobc']):
        self.directory = directory
        self.baselines = baselines
        self.moleculeName = moleculeName
        self.info = {}
        self.geneticRoundInfo = {}
    
    def graphResults(self):
        """
        Graphs the results of the experiment.

        Partly taken from https://www.pythoncharts.com/matplotlib/beautiful-bar-charts-matplotlib/
        """
        # mpl.font_manager._rebuild()
        fig, ax = plt.subplots(dpi=250, figsize=(12, 5))
        
        #Include tight layout
        fig.tight_layout()

        lowestEnergy = []
        averageEnergy = []
        std = []
        nCands = []

        for (name, energyList) in self.energyResults:
            lowestEnergy.append(min(energyList))

        mostMin = min(lowestEnergy)

        for (name, energyList) in self.energyResults:
            nCands.append(len(energyList))
            averageEnergy.append((sum(energyList) / len(energyList)) - mostMin)
            if len(energyList) > 1:
                std.append(stats.stdev(energyList)/(float(3)))
            else: std.append(0)


        # Sort baselines and average energy
        self.baselines, averageEnergy, std, nCands, lowestEnergy, self.timing = zip(*sorted(zip(self.baselines, averageEnergy, std, nCands, lowestEnergy, self.timing), key=lambda x: x[1]))

        tlabels = []
        for i, baseline in enumerate(self.baselines):
            minutes = self.timing[i].minute
            if minutes < 10:
                minutes = f"0{minutes}"
            # print(nCands[i])
            if baseline == 'setConformerSearch': baseline = 'provided'
            # else: labelName = nCands[i]
            tlabels.append(baseline + f" ({nCands[i]}, {self.timing[i].hour}:{minutes})")

        

        bars = ax.bar(
            x=np.arange(len(averageEnergy)),
            height=averageEnergy,
            tick_label=tlabels,
            yerr=std,
            alpha=0.6,
            capsize=5
        )
        
        fbar = bars[0].get_height()

        for i, bar in enumerate(bars):
            break
            # Write at bottom of every bar the value of the bar
            ax.text(
                bar.get_x() + bar.get_width() / 2,
                bar.get_height()  + fbar*0.02,
                f" {len(self.energyResults[i][1])}",
                ha='left',
                color='black',
                fontsize=15
            )
            # Change label of bar


        # Axis formatting.
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.spines['bottom'].set_color('#DDDDDD')
        ax.tick_params(bottom=False, left=False)
        ax.set_axisbelow(True)
        ax.yaxis.grid(True, color='#EEEEEE')
        # Set yaxis to be log
        ax.set_yscale('log')
        ax.xaxis.grid(False)

        # ax2.set_yticks(ax.get_yticks())
        # ax.set_yscale('log')
        # ax.set_ylim(ax.get_ylim())
        # ax.set_axisbelow(True)
        print(ax.get_ylim())
        x = np.arange(len(averageEnergy))
        print(x)
        y = [energy-mostMin + ax.get_ylim()[0]/5 for energy in lowestEnergy]
        # Set bottom of y lim to be the negative of the current bottom
        # ax.set_ylim(bottom=-ax.get_ylim()[0])

        print(y)
        ax.scatter(x, y, color='red', marker='o', zorder=10)

        # Add text annotations to the top of the bars.
        bar_color = bars[0].get_facecolor()
        fbar = 0

        # Add labels and a title.
        ax.set_xlabel('Method of Computation', labelpad=15, color='#333333')
        ax.set_ylabel('Difference to Minimum', labelpad=15, color='#333333')
        ax.set_title(f'Comparison of Conformer Search Methods for {self.moleculeName}', pad=15, color='#333333',
                    weight='bold')
        # Save figure

        # Set ax2 to have same y axis as ax

        # If doesn't exist create a graph directory
        if not os.path.exists(f'{self.directory}/graphs'):
            os.makedirs(f'{self.directory}/graphs')

        print(f"Save dir: {self.directory}/graphs/energy_comparison.pdf")
        plt.savefig(f'{self.directory}/graphs/energy_comparison.pdf', dpi=300, bbox_inches='tight')

    def graphGeneticProgress(self, geneticDir, name):
        directory = f"{geneticDir}/{name}/rounds"

        roundDirs = os.listdir(directory)

        self.geneticRoundInfo[name] = []

        for roundDir in roundDirs:
            with open(f'{directory}/{roundDir}/round_info.json', 'r') as f:
                roundData = json.load(f)
            
            self.geneticRoundInfo[name].append(roundData)
        
        # mpl.font_manager._rebuild()
        fig, ax = plt.subplots(dpi=250, figsize=(10, 5))

        # Line plot with marker over the average energy for each round data

        # Include tight layout
        fig.tight_layout()

        # Axis formatting.
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.spines['bottom'].set_color('#DDDDDD')
        ax.tick_params(bottom=False, left=False)
        ax.set_axisbelow(True)

        ax.yaxis.grid(True, color='#EEEEEE')
        ax.xaxis.grid(False)

        # Add labels and a title.
        ax.set_xlabel('Round', labelpad=15, color='#333333')
        ax.set_ylabel('Average Energy', labelpad=15, color='#333333')
        ax.set_title(f'Genetic Algorithm Progress for {name}', pad=15, color='#333333',
                    weight='bold')

        # Plot the data for each round
        plt.scatter(
            [i for i in range(len(self.geneticRoundInfo[name]))],
            [self.geneticRoundInfo[name][i]['energy_stats']['avg'] for i in range(len(self.geneticRoundInfo[name]))],
            color='blue',
            marker='o'
        )

        plt.errorbar(
            [i for i in range(len(self.geneticRoundInfo[name]))],
            [self.geneticRoundInfo[name][i]['energy_stats']['avg'] for i in range(len(self.geneticRoundInfo[name]))],
            yerr=[self.geneticRoundInfo[name][i]['energy_stats']['std'] for i in range(len(self.geneticRoundInfo[name]))],
            color='black',
            capsize=5
        )

        plt.show()
    
    def openInfo(self, name):
        if name.startswith('genetic'):
            nname = 'genetic'
        else:
            nname = name
        # with open(f'{self.directory}/{name}/{nname}_info.json', 'r') as f:
        #     self.info[name] = json.load(f)

    def getEnergies(self, name):
        structureNames = os.listdir(f'{self.directory}/{name}')
        energy_order_pairs = []

        for structureName in structureNames:
            if structureName.endswith('.sdf'):
                order = int(structureName.split('_')[0])
                num = float(structureName.split('-')[-1][:-4])
                energy_order_pairs.append((-num, order))
        
        if name == "setConformerSearch" or name == "genetic":
            sorted_energy_order_pairs = sorted(energy_order_pairs, key=lambda x: abs(x[0]), reverse=True)
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

        for (name, energies) in self.energyResults:
            newEnergies = copy.deepcopy(energies)
            
            truncatedResults.append((name, newEnergies[:numConfs]))
        
        self.energyResults = truncatedResults
            