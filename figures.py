import matplotlib.pyplot as plt
from pandas import *
from scipy import stats
import numpy as np
import math
import random

def gen_CIs(count_list, num_counts, num_permutations):
    
    lower_CI_list = []
    upper_CI_list = []
    lower_cutoff = int(num_permutations*0.025)
    upper_cutoff = int(num_permutations*0.975)
    
    for ct in count_list:
        
        sample_counts = [np.random.binomial(num_counts,ct/float(num_counts)) for i in range(num_permutations)]
        sample_counts.sort()
        
        lower_CI_list.append(abs(sample_counts[lower_cutoff]-ct))
        upper_CI_list.append(abs(sample_counts[upper_cutoff]-ct))
    
    return [list(lower_CI_list),list(upper_CI_list)]

random.seed(38737)

data_avida = read_csv("../data_avida.csv")
data_largegenomes_analysis = read_csv("../data_largegenomes_analysis.csv")
data_stochastic = read_csv("../data_stochastic.csv")

data_avida_norevertlethal = data_avida[data_avida["Revert"] == "None"]
data_avida_varlength_slip = data_avida_norevertlethal[data_avida_norevertlethal["Treatment"] == "varlength_slip"]
data_avida_fixedlength = data_avida_norevertlethal[data_avida_norevertlethal["Treatment"] == "fixedlength"]

data_avida_varlength_slip_001 = data_avida_varlength_slip[data_avida_varlength_slip["MutationRate"] == 0.01]
data_avida_varlength_slip_01 = data_avida_varlength_slip[data_avida_varlength_slip["MutationRate"] == 0.1]
data_avida_fixedlength_001 = data_avida_fixedlength[data_avida_fixedlength["MutationRate"] == 0.01]
data_avida_fixedlength_01 = data_avida_fixedlength[data_avida_fixedlength["MutationRate"] == 0.1]

data_avida_varlength_slip_001_extinct = data_avida_varlength_slip_001[data_avida_varlength_slip_001["Extinction"] == 1]
data_avida_varlength_slip_01_extinct = data_avida_varlength_slip_01[data_avida_varlength_slip_01["Extinction"] == 1]
data_avida_fixedlength_001_extinct = data_avida_fixedlength_001[data_avida_fixedlength_001["Extinction"] == 1]
data_avida_fixedlength_01_extinct = data_avida_fixedlength_01[data_avida_fixedlength_01["Extinction"] == 1]

data_avida_varlength_slip_001_survived = data_avida_varlength_slip_001[data_avida_varlength_slip_001["Extinction"] == 0]
data_avida_varlength_slip_01_survived = data_avida_varlength_slip_01[data_avida_varlength_slip_01["Extinction"] == 0]
data_avida_fixedlength_001_survived = data_avida_fixedlength_001[data_avida_fixedlength_001["Extinction"] == 0]
data_avida_fixedlength_01_survived = data_avida_fixedlength_01[data_avida_fixedlength_01["Extinction"] == 0]

data_avida_varlength_slip_001_nonzerofitness = data_avida_varlength_slip_001[data_avida_varlength_slip_001["FinalFitness"] != 0]
data_avida_varlength_slip_001_extinct_nonzerofitness = data_avida_varlength_slip_001_extinct[data_avida_varlength_slip_001_extinct["FinalFitness"] != 0]
data_avida_varlength_slip_001_survived_nonzerofitness = data_avida_varlength_slip_001_survived[data_avida_varlength_slip_001_survived["FinalFitness"] != 0]

data_avida_revertlethal = data_avida[data_avida["Revert"] == "Lethal"]
data_avida_revertlethal_varlength_slip = data_avida_revertlethal[data_avida_revertlethal["Treatment"] == "varlength_slip"]
data_avida_revertlethal_fixedlength = data_avida_revertlethal[data_avida_revertlethal["Treatment"] == "fixedlength"]
data_avida_revertlethal_varlength_slip_001 = data_avida_revertlethal_varlength_slip[data_avida_revertlethal_varlength_slip["MutationRate"] == 0.01]
data_avida_revertlethal_fixedlength_001 = data_avida_revertlethal_fixedlength[data_avida_revertlethal_fixedlength["MutationRate"] == 0.01]

data_avida_revertlethal_varlength_slip_001_extinct = data_avida_revertlethal_varlength_slip_001[data_avida_revertlethal_varlength_slip_001["Extinction"] == 1]
data_avida_revertlethal_fixedlength_001_extinct = data_avida_revertlethal_fixedlength_001[data_avida_revertlethal_fixedlength_001["Extinction"] == 1]
data_avida_revertlethal_varlength_slip_001_survived = data_avida_revertlethal_varlength_slip_001[data_avida_revertlethal_varlength_slip_001["Extinction"] == 0]
data_avida_revertlethal_fixedlength_001_survived = data_avida_revertlethal_fixedlength_001[data_avida_revertlethal_fixedlength_001["Extinction"] == 0]
data_avida_revertlethal_varlength_slip_001_survived_nonzerofitness = data_avida_revertlethal_varlength_slip_001_survived[data_avida_revertlethal_varlength_slip_001_survived["FinalFitness"] != 0]
data_avida_revertlethal_varlength_slip_001_extinct_nonzerofitness = data_avida_revertlethal_varlength_slip_001_extinct[data_avida_revertlethal_varlength_slip_001_extinct["FinalFitness"] != 0]
data_avida_revertlethal_varlength_slip_001_nonzerofitness = data_avida_revertlethal_varlength_slip_001[data_avida_revertlethal_varlength_slip_001["FinalFitness"] != 0]

data_avida_revertdel = data_avida[data_avida["Revert"] == "Del"]
data_avida_revertdel_varlength_slip = data_avida_revertdel[data_avida_revertdel["Treatment"] == "varlength_slip"]
data_avida_revertdel_fixedlength = data_avida_revertdel[data_avida_revertdel["Treatment"] == "fixedlength"]
data_avida_revertdel_varlength_slip_001 = data_avida_revertdel_varlength_slip[data_avida_revertdel_varlength_slip["MutationRate"] == 0.01]
data_avida_revertdel_fixedlength_001 = data_avida_revertdel_fixedlength[data_avida_revertdel_fixedlength["MutationRate"] == 0.01]

data_avida_revertdel_varlength_slip_001_extinct = data_avida_revertdel_varlength_slip_001[data_avida_revertdel_varlength_slip_001["Extinction"] == 1]
data_avida_revertdel_fixedlength_001_extinct = data_avida_revertdel_fixedlength_001[data_avida_revertdel_fixedlength_001["Extinction"] == 1]
data_avida_revertdel_varlength_slip_001_survived = data_avida_revertdel_varlength_slip_001[data_avida_revertdel_varlength_slip_001["Extinction"] == 0]
data_avida_revertdel_fixedlength_001_survived = data_avida_revertdel_fixedlength_001[data_avida_revertdel_fixedlength_001["Extinction"] == 0]
data_avida_revertdel_varlength_slip_001_survived_nonzerofitness = data_avida_revertdel_varlength_slip_001_survived[data_avida_revertdel_varlength_slip_001_survived["FinalFitness"] != 0]
data_avida_revertdel_varlength_slip_001_extinct_nonzerofitness = data_avida_revertdel_varlength_slip_001_extinct[data_avida_revertdel_varlength_slip_001_extinct["FinalFitness"] != 0]
data_avida_revertdel_varlength_slip_001_nonzerofitness = data_avida_revertdel_varlength_slip_001[data_avida_revertdel_varlength_slip_001["FinalFitness"] != 0]

data_stochastic_revertnone = data_stochastic[data_stochastic["RevertStatus"] == "revertnone"]
data_stochastic_revertlethal = data_stochastic[data_stochastic["RevertStatus"] == "revertlethal"]
data_stochastic_revertnone_extinct = data_stochastic_revertnone[data_stochastic_revertnone["Extinction"] == 1]
data_stochastic_revertlethal_extinct = data_stochastic_revertlethal[data_stochastic_revertlethal["Extinction"] == 1]
data_stochastic_revertnone_survived = data_stochastic_revertnone[data_stochastic_revertnone["Extinction"] == 0]
data_stochastic_revertlethal_survived = data_stochastic_revertlethal[data_stochastic_revertlethal["Extinction"] == 0]
data_stochastic_revertlethal_allviable = data_stochastic[data_stochastic["PercentNonViableAnalyses"] == 0]
data_stochastic_revertlethal_someviable = data_stochastic[data_stochastic["PercentNonViableAnalyses"] != 0]

lowmut_range = [5,6,7,8,10,15,20]
lowmut_range_revertlethal = [5,6,7,8,10,15,20]
highmut_range = [10,12,15,16,17,20,25]

def figure1():
    
    varlength_slip_001_extinct_counts = [len(data_avida_varlength_slip_001_extinct[data_avida_varlength_slip_001_extinct["PopulationSize"] == N]) for N in lowmut_range]
    fixedlength_001_extinct_counts = [len(data_avida_fixedlength_001_extinct[data_avida_fixedlength_001_extinct["PopulationSize"] == N]) for N in lowmut_range]
    
    varlength_slip_01_extinct_counts = [len(data_avida_varlength_slip_01_extinct[data_avida_varlength_slip_01_extinct["PopulationSize"] == N]) for N in highmut_range]
    fixedlength_01_extinct_counts = [len(data_avida_fixedlength_01_extinct[data_avida_fixedlength_01_extinct["PopulationSize"] == N]) for N in highmut_range]

    fig,(a,b) = plt.subplots(1,2,figsize=(7.5,2.5))
    a.errorbar(lowmut_range,varlength_slip_001_extinct_counts,yerr=gen_CIs(varlength_slip_001_extinct_counts,100,10000),fmt="-ko",capsize=2.5)
    a.errorbar(lowmut_range,fixedlength_001_extinct_counts,yerr=gen_CIs(fixedlength_001_extinct_counts,100,10000),fmt="--ko",capsize=2.5)
    a.errorbar(highmut_range,varlength_slip_01_extinct_counts,yerr=gen_CIs(varlength_slip_01_extinct_counts,100,10000),fmt="-ko",markerfacecolor="w",capsize=2.5)
    a.errorbar(highmut_range,fixedlength_01_extinct_counts,yerr=gen_CIs(fixedlength_01_extinct_counts,100,10000),fmt="--ko",markerfacecolor="w",capsize=2.5)

    b.errorbar([5,6,7,8,10],
           [np.mean(data_avida_varlength_slip_001_extinct[data_avida_varlength_slip_001_extinct["PopulationSize"] == N]["TimeToExtinction"]) for N in [5,6,7,8,10]],
           yerr=[2.0*np.std(data_avida_varlength_slip_001_extinct[data_avida_varlength_slip_001_extinct["PopulationSize"] == N]["TimeToExtinction"])/math.sqrt(len(data_avida_varlength_slip_001_extinct[data_avida_varlength_slip_001_extinct["PopulationSize"] == N]["TimeToExtinction"])) for N in [5,6,7,8,10]],
           fmt="-ko")
    b.errorbar([10,12,15,16,17],
           [np.mean(data_avida_varlength_slip_01_extinct[data_avida_varlength_slip_01_extinct["PopulationSize"] == N]["TimeToExtinction"]) for N in [10,12,15,16,17]],
           yerr=[2.0*np.std(data_avida_varlength_slip_01_extinct[data_avida_varlength_slip_01_extinct["PopulationSize"] == N]["TimeToExtinction"])/math.sqrt(len(data_avida_varlength_slip_01_extinct[data_avida_varlength_slip_01_extinct["PopulationSize"] == N]["TimeToExtinction"])) for N in [10,12,15,16,17]],
           fmt="-ko",markerfacecolor="w")
    b.errorbar([5,6],
           [np.mean(data_avida_fixedlength_001_extinct[data_avida_fixedlength_001_extinct["PopulationSize"] == N]["TimeToExtinction"]) for N in [5,6]],
           yerr=[2.0*np.std(data_avida_fixedlength_001_extinct[data_avida_fixedlength_001_extinct["PopulationSize"] == N]["TimeToExtinction"])/math.sqrt(len(data_avida_fixedlength_001_extinct[data_avida_fixedlength_001_extinct["PopulationSize"] == N]["TimeToExtinction"])) for N in [5,6]],
           fmt="--ko")
    b.errorbar([10,12,15,16,17],
           [np.mean(data_avida_fixedlength_01_extinct[data_avida_fixedlength_01_extinct["PopulationSize"] == N]["TimeToExtinction"]) for N in [10,12,15,16,17]],
           yerr=[2.0*np.std(data_avida_fixedlength_01_extinct[data_avida_fixedlength_01_extinct["PopulationSize"] == N]["TimeToExtinction"])/math.sqrt(len(data_avida_fixedlength_01_extinct[data_avida_fixedlength_01_extinct["PopulationSize"] == N]["TimeToExtinction"])) for N in [10,12,15,16,17]],
           fmt="--ko",markerfacecolor="w")
    
    a.set_xlim(0,30)
    a.set_xticks(range(5,30,5))
    a.set_xticklabels([str(i) for i in range(5,30,5)],fontsize=8)
    a.set_yticks(range(0,120,20))
    a.set_yticklabels([str(i) for i in range(0,120,20)],fontsize=8)
    a.set_xlabel("Population Size",fontsize=10)
    a.set_ylabel("Extinct Populations",fontsize=10)
    b.set_xticks(range(5,25,5))
    b.set_yticks(range(0,100000,20000))
    b.set_yticklabels([str(i) for i in range(0,10,2)],fontsize=8)
    b.set_xlabel("Population Size",fontsize=10)
    b.set_ylabel("Generations to Extinction ($x10^{4}$)",fontsize=10)
    b.set_xticklabels([str(i) for i in range(5,25,5)],fontsize=8)
    plt.text(0,75000,"A",fontsize=12)
    plt.text(19,75000,"B",fontsize=12)
    plt.savefig("figure1.eps",bbox_inches = "tight", dpi=600)
    
    plt.show()

def figure2():
    
    print("N=7, U=0.15 ",stats.mannwhitneyu(list(data_avida_varlength_slip_001_extinct[data_avida_varlength_slip_001_extinct["PopulationSize"] == 7]["FinalGenomeSize"]),
                                          list(data_avida_varlength_slip_001_survived[data_avida_varlength_slip_001_survived["PopulationSize"] == 7]["FinalGenomeSize"])))
    print("N=8, U=0.15 ",stats.mannwhitneyu(list(data_avida_varlength_slip_001_extinct[data_avida_varlength_slip_001_extinct["PopulationSize"] == 8]["FinalGenomeSize"]),
                                          list(data_avida_varlength_slip_001_survived[data_avida_varlength_slip_001_survived["PopulationSize"] == 8]["FinalGenomeSize"])))
    print("N=10, U=0.15 ",stats.mannwhitneyu(list(data_avida_varlength_slip_001_extinct[data_avida_varlength_slip_001_extinct["PopulationSize"] == 10]["FinalGenomeSize"]),
                                          list(data_avida_varlength_slip_001_survived[data_avida_varlength_slip_001_survived["PopulationSize"] == 10]["FinalGenomeSize"])))
    
    print("N=15, U=1.5 ",stats.mannwhitneyu(list(data_avida_varlength_slip_01_extinct[data_avida_varlength_slip_01_extinct["PopulationSize"] == 15]["FinalGenomeSize"]),
                                          list(data_avida_varlength_slip_01_survived[data_avida_varlength_slip_01_survived["PopulationSize"] == 15]["FinalGenomeSize"])))
    print("N=16, U=1.5 ",stats.mannwhitneyu(list(data_avida_varlength_slip_01_extinct[data_avida_varlength_slip_01_extinct["PopulationSize"] == 16]["FinalGenomeSize"]),
                                          list(data_avida_varlength_slip_01_survived[data_avida_varlength_slip_01_survived["PopulationSize"] == 16]["FinalGenomeSize"])))
    print("N=17, U=1.5 ",stats.mannwhitneyu(list(data_avida_varlength_slip_01_extinct[data_avida_varlength_slip_01_extinct["PopulationSize"] == 17]["FinalGenomeSize"]),
                                          list(data_avida_varlength_slip_01_survived[data_avida_varlength_slip_01_survived["PopulationSize"] == 17]["FinalGenomeSize"])))
    
    fig,(a,b) = plt.subplots(1,2,figsize=(7.5,2.5))

    extinct_lowmut = a.boxplot([list(data_avida_varlength_slip_001_extinct[data_avida_varlength_slip_001_extinct["PopulationSize"] == N]["FinalGenomeSize"]) for N in [5,6,7,8,10]],
           positions=[4.75,5.75,6.75,7.75,9.75],patch_artist=True)
    survive_lowmut = a.boxplot([list(data_avida_varlength_slip_001_survived[data_avida_varlength_slip_001_survived["PopulationSize"] == N]["FinalGenomeSize"]) for N in [7,8,10,15,20]],
           positions=[7.25,8.25,10.25,15.25,20.25],patch_artist=True)
    
    a.plot([6.75,7.25],[800,800],"k-")
    a.plot([6.75,6.75],[775,800],"k-")
    a.plot([7.25,7.25],[775,800],"k-")
    
    a.plot([7.75,8.25],[800,800],"k-")
    a.plot([7.75,7.75],[775,800],"k-")
    a.plot([8.25,8.25],[775,800],"k-")
    
    a.plot([9.75,10.25],[800,800],"k-")
    a.plot([9.75,9.75],[775,800],"k-")
    a.plot([10.25,10.25],[775,800],"k-")
    
    for patch, color in zip(extinct_lowmut['boxes'], ["white" for j in range(5)]):
        patch.set_facecolor(color)
    for patch, color in zip(survive_lowmut['boxes'], ["gray" for j in range(5)]):
        patch.set_facecolor(color)
        
    extinct_highmut = b.boxplot([list(data_avida_varlength_slip_01_extinct[data_avida_varlength_slip_01_extinct["PopulationSize"] == N]["FinalGenomeSize"]) for N in [10,12,15,16,17]],
           positions=[9.75,11.75,14.75,15.75,16.75],patch_artist=True)
    survive_highmut = b.boxplot([list(data_avida_varlength_slip_01_survived[data_avida_varlength_slip_01_survived["PopulationSize"] == N]["FinalGenomeSize"]) for N in [15,16,17,20,25]],
           positions=[15.25,16.25,17.25,20.25,25.25],patch_artist=True)
    
    for patch, color in zip(extinct_highmut['boxes'], ["white" for j in range(5)]):
        patch.set_facecolor(color)
    for patch, color in zip(survive_highmut['boxes'], ["gray" for j in range(5)]):
        patch.set_facecolor(color)
        
    b.plot([14.75,15.25],[800,800],"k-")
    b.plot([14.75,14.75],[775,800],"k-")
    b.plot([15.25,15.25],[775,800],"k-")

    b.plot([15.75,16.25],[800,800],"k-")
    b.plot([15.75,15.75],[775,800],"k-")
    b.plot([16.25,16.25],[775,800],"k-")
    
    b.plot([16.75,17.25],[800,800],"k-")
    b.plot([16.75,16.75],[775,800],"k-")
    b.plot([17.25,17.25],[775,800],"k-")

    a.set_xlim(4,21)
    a.set_ylabel("Genome Size",fontsize=10)
    a.set_xticks([5,6,7,8,10,15,20])
    a.set_xticklabels([5,6,7,8,10,15,20],fontsize=8)
    a.set_xlabel("Population Size",fontsize=10)
    a.set_ylim(-40,1010)
    a.set_yticks(range(0,1010,200))
    a.set_yticklabels([str(i) for i in range(0,1010,200)],fontsize=8)
    b.set_xlim(9,26)
    b.set_ylim(-40,1010)
    b.set_xticks([10,12,15,16,17,20,25])
    b.set_xticklabels([10,12,15,16,17,20,25],fontsize=8)
    b.set_xlabel("Population Size",fontsize=10)
    b.set_yticks(range(0,1010,200))
    b.set_yticklabels([str(i) for i in range(0,1010,200)],fontsize=8)
    
    plt.text(-11,900,"A",fontsize=12)
    plt.text(9.5,900,"B",fontsize=12)
    plt.text(-5.85,810,"**")
    plt.text(-7.85,810,"**")
    plt.text(-8.85,810,"**")
    plt.text(14.75,810,"*")
    plt.text(15.55,820,"N.S.",fontsize=6)
    plt.text(16.60,820,"N.S.",fontsize=6)
    plt.savefig("figure2.eps",bbox_inches = "tight", dpi=600)
    plt.show()
    
def figure3():
    
    print("All",stats.spearmanr(data_avida_varlength_slip_001_nonzerofitness["FinalGenomeSize"],data_avida_varlength_slip_001_nonzerofitness["FinalMuLethal"]))
    print("Extinct",stats.spearmanr(data_avida_varlength_slip_001_extinct_nonzerofitness["FinalGenomeSize"],data_avida_varlength_slip_001_extinct_nonzerofitness["FinalMuLethal"]))
    print("Survived",stats.spearmanr(data_avida_varlength_slip_001_survived_nonzerofitness["FinalGenomeSize"],data_avida_varlength_slip_001_survived_nonzerofitness["FinalMuLethal"]))
    
    print("N=7, U=0.15 ",stats.mannwhitneyu(list(data_avida_varlength_slip_001_extinct[data_avida_varlength_slip_001_extinct["PopulationSize"] == 7]["FinalMuLethal"]),
                                          list(data_avida_varlength_slip_001_survived[data_avida_varlength_slip_001_survived["PopulationSize"] == 7]["FinalMuLethal"])))
    print("N=8, U=0.15 ",stats.mannwhitneyu(list(data_avida_varlength_slip_001_extinct[data_avida_varlength_slip_001_extinct["PopulationSize"] == 8]["FinalMuLethal"]),
                                          list(data_avida_varlength_slip_001_survived[data_avida_varlength_slip_001_survived["PopulationSize"] == 8]["FinalMuLethal"])))
    print("N=10, U=0.15 ",stats.mannwhitneyu(list(data_avida_varlength_slip_001_extinct[data_avida_varlength_slip_001_extinct["PopulationSize"] == 10]["FinalMuLethal"]),
                                          list(data_avida_varlength_slip_001_survived[data_avida_varlength_slip_001_survived["PopulationSize"] == 10]["FinalMuLethal"])))

    fig,(a,b) = plt.subplots(1,2,figsize=(7.5,2.5))
    
    a.plot(data_avida_varlength_slip_001_extinct_nonzerofitness["FinalGenomeSize"],data_avida_varlength_slip_001_extinct_nonzerofitness["FinalMuLethal"],"ko",markerfacecolor="w")
    a.plot(data_avida_varlength_slip_001_survived_nonzerofitness["FinalGenomeSize"],data_avida_varlength_slip_001_survived_nonzerofitness["FinalMuLethal"],"ko",alpha=0.5)

    extinct_lowmut = b.boxplot([list(data_avida_varlength_slip_001_extinct_nonzerofitness[data_avida_varlength_slip_001_extinct_nonzerofitness["PopulationSize"] == N]["FinalMuLethal"]) for N in [5,6,7,8,10]],
           positions=[4.75,5.75,6.75,7.75,9.75],patch_artist=True)
    survive_lowmut = b.boxplot([list(data_avida_varlength_slip_001_survived_nonzerofitness[data_avida_varlength_slip_001_survived_nonzerofitness["PopulationSize"] == N]["FinalMuLethal"]) for N in [7,8,10,15,20]],
           positions=[7.25,8.25,10.25,15.25,20.25],patch_artist=True)
    
    b.plot([6.75,7.25],[1.3,1.3],"k-")
    b.plot([6.75,6.75],[1.25,1.3],"k-")
    b.plot([7.25,7.25],[1.25,1.3],"k-")
    
    b.plot([7.75,8.25],[1.3,1.3],"k-")
    b.plot([7.75,7.75],[1.25,1.3],"k-")
    b.plot([8.25,8.25],[1.25,1.3],"k-")
    
    b.plot([9.75,10.25],[1.3,1.3],"k-")
    b.plot([9.75,9.75],[1.25,1.3],"k-")
    b.plot([10.25,10.25],[1.25,1.3],"k-")
    
    for patch, color in zip(extinct_lowmut['boxes'], ["white" for j in range(5)]):
        patch.set_facecolor(color)
    for patch, color in zip(survive_lowmut['boxes'], ["gray" for j in range(5)]):
        patch.set_facecolor(color)
        
    a.set_xlabel("Genome Size",fontsize=10)
    a.set_ylabel("Lethal Mutation Rate",fontsize=10)
    a.set_yticks([0.5,1.0,1.5])
    a.set_yticklabels(["0.5","1.0","1.5"],fontsize=8)
    a.set_xticks(range(0,1010,200))
    a.set_xticklabels([str(i) for i in range(0,1010,200)],fontsize=8)
    b.set_yticks([0.5,1.0,1.5])
    b.set_yticklabels(["0.5","1.0","1.5"],fontsize=8)
    b.set_xlim(4,21)
    b.set_xticks([5,6,7,8,10,15,20])
    b.set_xticklabels([5,6,7,8,10,15,20],fontsize=8)
    b.set_xlabel("Population Size",fontsize=10)

    
    plt.text(-16.25,1.55,"A",fontsize=12)
    plt.text(4.2,1.55,"B",fontsize=12)
    plt.text(6.6,1.35,"**")
    plt.text(7.6,1.35,"**")
    plt.text(9.6,1.35,"**")
    
    plt.savefig("figure3.eps",bbox_inches = "tight", dpi=600)
    plt.show()
    
def figure4():
    
    print("Constant",np.mean(data_largegenomes_analysis[data_largegenomes_analysis["ChangeGenomeSize"] == 0]["ChangeLethalMutationRate"]),
          1.96*np.std(data_largegenomes_analysis[data_largegenomes_analysis["ChangeGenomeSize"] == 0]["ChangeLethalMutationRate"])/math.sqrt(len(data_largegenomes_analysis[data_largegenomes_analysis["ChangeGenomeSize"] == 0]["ChangeLethalMutationRate"])))
    print("Increase",np.mean(data_largegenomes_analysis[data_largegenomes_analysis["ChangeGenomeSize"] > 0]["ChangeLethalMutationRate"]),
          1.96*np.std(data_largegenomes_analysis[data_largegenomes_analysis["ChangeGenomeSize"] > 0]["ChangeLethalMutationRate"])/math.sqrt(len(data_largegenomes_analysis[data_largegenomes_analysis["ChangeGenomeSize"] > 0]["ChangeLethalMutationRate"])))
    print("Decrease",np.mean(data_largegenomes_analysis[data_largegenomes_analysis["ChangeGenomeSize"] < 0]["ChangeLethalMutationRate"]),
          1.96*np.std(data_largegenomes_analysis[data_largegenomes_analysis["ChangeGenomeSize"] < 0]["ChangeLethalMutationRate"])/math.sqrt(len(data_largegenomes_analysis[data_largegenomes_analysis["ChangeGenomeSize"] < 0]["ChangeLethalMutationRate"])))
    
    print("Correlation",stats.spearmanr(data_largegenomes_analysis[data_largegenomes_analysis["ChangeGenomeSize"] != 0]["ChangeGenomeSize"],
       data_largegenomes_analysis[data_largegenomes_analysis["ChangeGenomeSize"] != 0]["ChangeLethalMutationRate"]))


    fig,(a,b) = plt.subplots(1,2,figsize=(7.5,2.5))
    a.boxplot([data_largegenomes_analysis[data_largegenomes_analysis["ChangeGenomeSize"] == 0]["ChangeLethalMutationRate"],
              data_largegenomes_analysis[data_largegenomes_analysis["ChangeGenomeSize"] > 0]["ChangeLethalMutationRate"],
              data_largegenomes_analysis[data_largegenomes_analysis["ChangeGenomeSize"] < 0]["ChangeLethalMutationRate"]])
    a.set_xticklabels(["Genome\nSize\nConstant","Genome\nSize\nIncreased","Genome\nSize\nDecreased"],fontsize=10)
    a.set_ylabel("Change in Lethal Mutation Rate",fontsize=10)
    a.set_ylim(-0.8,0.8)
    a.set_yticks([-0.5,0.0,0.5])
    a.set_yticklabels(["-0.5","0.0","0.5"],fontsize=8)
    b.plot(data_largegenomes_analysis[data_largegenomes_analysis["ChangeGenomeSize"] != 0]["ChangeGenomeSize"],
           data_largegenomes_analysis[data_largegenomes_analysis["ChangeGenomeSize"] != 0]["ChangeLethalMutationRate"],"ko",alpha=0.5)
    b.plot([0,0],[-1,1],"k--")
    b.plot([-900,400],[0,0],"k--")
    b.set_ylim(-0.8,0.8)
    b.set_xlim(-900,400)
    b.set_xlabel("Change in Genome Size",fontsize=10)
    b.set_yticks([-0.5,0.0,0.5])
    b.set_yticklabels(["-0.5","0.0","0.5"],fontsize=8)
    b.set_xticks(range(-750,300,250))
    b.set_xticklabels([str(i) for i in range(-750,300,250)],fontsize=8)
    plt.text(-2450,0.65,"A",fontsize=12)
    plt.text(-890,0.65,"B",fontsize=12)
    plt.savefig("figure4.eps",bbox_inches = "tight", dpi=600)
    plt.show()

def figure5():
    varlength_slip_001_extinct_counts = [len(data_avida_varlength_slip_001_extinct[data_avida_varlength_slip_001_extinct["PopulationSize"] == N]) for N in lowmut_range]
    fixedlength_001_extinct_counts = [len(data_avida_fixedlength_001_extinct[data_avida_fixedlength_001_extinct["PopulationSize"] == N]) for N in lowmut_range]
 
    varlength_slip_revertlethal_001_extinct_counts = [len(data_avida_revertlethal_varlength_slip_001_extinct[data_avida_revertlethal_varlength_slip_001_extinct["PopulationSize"] == N]) for N in lowmut_range]   
    fixedlength_revertlethal_001_extinct_counts = [len(data_avida_revertlethal_fixedlength_001_extinct[data_avida_revertlethal_fixedlength_001_extinct["PopulationSize"] == N]) for N in lowmut_range]
    
    varlength_slip_revertdel_001_extinct_counts = [len(data_avida_revertdel_varlength_slip_001_extinct[data_avida_revertdel_varlength_slip_001_extinct["PopulationSize"] == N]) for N in lowmut_range]
    fixedlength_revertdel_001_extinct_counts = [len(data_avida_revertdel_fixedlength_001_extinct[data_avida_revertdel_fixedlength_001_extinct["PopulationSize"] == N]) for N in lowmut_range]
    
    print("Stochastic Viable Larger Genomes Test",np.median(data_stochastic_revertlethal_allviable["GenomeSize"]),
          np.median(data_stochastic_revertlethal_someviable["GenomeSize"]),
          stats.mannwhitneyu(data_stochastic_revertlethal_allviable["GenomeSize"],data_stochastic_revertlethal_someviable["GenomeSize"]))
    
    print("PercentNonViable_RevertNone_extinct = Any",len(data_stochastic_revertnone_extinct["PercentNonViableAnalyses"]))
    print("PercentNonViable_RevertNone_extinct = 0",len(data_stochastic_revertnone_extinct[data_stochastic_revertnone_extinct["PercentNonViableAnalyses"] == 0]))
    print("PercentNonViable_RevertNone_extinct = 100",len(data_stochastic_revertnone_extinct[data_stochastic_revertnone_extinct["PercentNonViableAnalyses"] == 1.0]))
    print("PercentNonViable_RevertLethal_extinct = Any",len(data_stochastic_revertlethal_extinct["PercentNonViableAnalyses"]))
    print("PercentNonViable_RevertLethal_extinct = 0",len(data_stochastic_revertlethal_extinct[data_stochastic_revertlethal_extinct["PercentNonViableAnalyses"] == 0]))
    print("PercentNonViable_RevertLethal_extinct = 100",len(data_stochastic_revertlethal_extinct[data_stochastic_revertlethal_extinct["PercentNonViableAnalyses"] == 1.0]))
    
    fig,((a,b),(c,d)) = plt.subplots(2,2,figsize=(7.5,6))
    a.errorbar(lowmut_range,varlength_slip_001_extinct_counts,yerr=gen_CIs(varlength_slip_001_extinct_counts,100,10000),fmt="-ko",capsize=2.5)
    a.errorbar(lowmut_range,varlength_slip_revertdel_001_extinct_counts,yerr=gen_CIs(varlength_slip_revertdel_001_extinct_counts,100,10000),fmt="-ks",capsize=2.5)
    a.errorbar(lowmut_range,fixedlength_001_extinct_counts,yerr=gen_CIs(fixedlength_001_extinct_counts,100,10000),fmt="--ko",capsize=2.5)
    a.errorbar(lowmut_range,fixedlength_revertdel_001_extinct_counts,yerr=gen_CIs(fixedlength_revertdel_001_extinct_counts,100,10000),fmt="--ks",capsize=2.5)
    a.set_xlabel("Population Size",fontsize=10)
    a.set_ylabel("Extinct Populations",fontsize=10)
    a.set_xticks(range(5,25,5))
    a.set_xticklabels([str(i) for i in range(5,25,5)],fontsize=8)
    a.set_yticks(range(0,120,20))
    a.set_yticklabels([str(i) for i in range(0,120,20)],fontsize=8)
    
    b.errorbar(lowmut_range,varlength_slip_001_extinct_counts,yerr=gen_CIs(varlength_slip_001_extinct_counts,100,10000),fmt="-ko",capsize=2.5)
    b.errorbar(lowmut_range,varlength_slip_revertlethal_001_extinct_counts,yerr=gen_CIs(varlength_slip_revertlethal_001_extinct_counts,100,10000),fmt="-k^",capsize=2.5)
    b.errorbar(lowmut_range,fixedlength_001_extinct_counts,yerr=gen_CIs(fixedlength_001_extinct_counts,100,10000),fmt="--ko",capsize=2.5)
    b.errorbar(lowmut_range,fixedlength_revertlethal_001_extinct_counts,yerr=gen_CIs(fixedlength_revertlethal_001_extinct_counts,100,10000),fmt="--k^",capsize=2.5)
    b.set_xlabel("Population Size",fontsize=10)
    b.set_ylabel("Extinct Populations",fontsize=10)
    b.set_xticks(range(5,25,5))
    b.set_xticklabels([str(i) for i in range(5,25,5)],fontsize=8)
    b.set_yticks(range(0,120,20))
    b.set_yticklabels([str(i) for i in range(0,120,20)],fontsize=8)
    
    c.boxplot([data_stochastic_revertlethal_survived["PercentNonViableAnalyses"],data_stochastic_revertlethal_extinct["PercentNonViableAnalyses"],
               data_stochastic_revertnone_survived["PercentNonViableAnalyses"],data_stochastic_revertnone_extinct["PercentNonViableAnalyses"]], positions = [0,1.25,2.75,4])
    c.set_xticklabels(["Revert-\nLethal\nSurvived","Revert-\nLethal\nExtinct","Original\nSurvived","Original\nExtinct"],fontsize=8)
    c.set_ylabel("Non-viable Trials (%)",fontsize=10)
    c.set_yticks([i/100.0 for i in range(0,120,20)])
    c.set_yticklabels([str(i/100.0) for i in range(0,120,20)],fontsize=8)
    c.set_xlim(-1,5)
    
    d.boxplot([data_stochastic_revertlethal_allviable["GenomeSize"],data_stochastic_revertlethal_someviable["GenomeSize"]])
    d.set_xticklabels(["All\nViable","Stochastic\nViable"],fontsize=10)
    d.set_ylabel("Genome Size (x$10^2$)",fontsize=10)
    d.set_yticks(range(0,1200,200))
    d.set_yticklabels([str(i) for i in range(0,12,2)],fontsize=8)
    
    plt.text(-0.02,2300,"A")
    plt.text(2.4,2300,"B")
    plt.text(-0.02,1000,"C")
    plt.text(2.4,1000,"D")
    
    plt.savefig("figure5.eps",bbox_inches = "tight", dpi=600)
    plt.show()
    

figure1()
#figure2()
#figure3()
#figure4()
figure5()