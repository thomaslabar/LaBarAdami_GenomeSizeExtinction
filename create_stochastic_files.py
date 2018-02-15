import pandas as pd
import copy

def GetFileLines(file_str):
    
    #Input: string of file name
    #Output: List of lines in file
    
    file_name=str(file_str)
    try:
        f=open(file_name)
    except:
        raise StandardError(file_name+" does not exist")
    lines=f.read().splitlines()
    f.close()
    assert len(lines)!=0,"length of "+str(file_str)+"=0"
    return copy.deepcopy(lines)

data_finalgenotype = pd.read_csv("data_finalgenotypes.csv")

data_finalgenotype_varlength_slip = data_finalgenotype[data_finalgenotype["Treatment"] == "varlength_slip"]

data_finalgenotype_varlength_slip_N5 = data_finalgenotype_varlength_slip[data_finalgenotype_varlength_slip["PopulationSize"] == 5]
data_finalgenotype_varlength_slip_N5_mu001 = data_finalgenotype_varlength_slip_N5[data_finalgenotype_varlength_slip_N5["MutationRate"] == 0.01]
data_finalgenotype_varlength_slip_N5_mu001_revertlethal = data_finalgenotype_varlength_slip_N5_mu001[data_finalgenotype_varlength_slip_N5_mu001["Revert"] == "Lethal"]

data_finalgenotype_varlength_slip_N8 = data_finalgenotype_varlength_slip[data_finalgenotype_varlength_slip["PopulationSize"] == 8]
data_finalgenotype_varlength_slip_N8_mu001 = data_finalgenotype_varlength_slip_N8[data_finalgenotype_varlength_slip_N8["MutationRate"] == 0.01]
data_finalgenotype_varlength_slip_N8_mu001_revertnone = data_finalgenotype_varlength_slip_N8_mu001[data_finalgenotype_varlength_slip_N8_mu001["Revert"] == "None"]

final_genotypes = {}
final_genotypes["revertnone"] = list(data_finalgenotype_varlength_slip_N8_mu001_revertnone["Sequence"])
final_genotypes["revertlethal"] = list(data_finalgenotype_varlength_slip_N5_mu001_revertlethal["Sequence"])
assert len(final_genotypes["revertnone"]) == 100
assert len(final_genotypes["revertlethal"]) == 100

final_replicates = {}
final_replicates["revertnone"] = list(data_finalgenotype_varlength_slip_N8_mu001_revertnone["Replicate"])
final_replicates["revertlethal"] = list(data_finalgenotype_varlength_slip_N5_mu001_revertlethal["Replicate"])

for revert in ["revertnone","revertlethal"]:

    f = open("config_files/analyze_stochastic_"+str(revert)+".cfg","w")
    f.write("VERBOSE\n")
    f.write("PURGE_BATCH\n")
    for seq in final_genotypes[revert]:
        f.write("LOAD_SEQUENCE "+seq+"\n")
    f.write("RENAME "+str(final_replicates[revert][0])+"\n")
    f.write("RECALC num_trials 1000\n")
    f.write("DETAIL finaldominant_stochastic_analysis_"+str(revert)+".dat id prob_viable\n")
    f.close()
    
    g = open("sub_files/avida_a_stochastic_"+str(revert)+".sub","w")
    g.write("#!/bin/bash -login\n")
    g.write("#PBS -l walltime=1:00:00,nodes=1:ppn=1,mem=2gb\n")
    g.write("#PBS -j oe\n")
    g.write("cd ${HOME}/Documents/genomesize_extinction/avida/config_files\n")
    g.write("mkdir -p  ../data_stochastic/"+str(revert)+"/stochastic_analysis\n")
    g.write("./avida -a -set WORLD_X 5")
    g.write(" -set DATA_DIR ../data_stochastic")
    g.write(" -set ANALYZE_FILE analyze_stochastic_"+str(revert)+".cfg")
    g.write(" -set WORLD_Y 1 -set BIRTH_METHOD 4 -set COPY_MUT_PROB 0.0 -set DIV_MUT_PROB 0.01")
    g.write(" -set DIVIDE_INS_PROB 0.0 -set DIVIDE_DEL_PROB 0.0 -set DIVIDE_SLIP_PROB 0.0 -set REQUIRE_EXACT_COPY 1")
    g.close()
    