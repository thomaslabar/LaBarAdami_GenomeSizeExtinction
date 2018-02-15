import copy
import pandas as pd

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

large_data = pd.read_csv("data_largegenomes.csv")

g=open("data_largegenomes_analysis.csv","wb")
g.write("ChangeGenomeSize,ChangeLethalMutationRate\n")

for t in ["varlength_slip"]:
    
    for mu in ["001"]:
        
        if mu == "001":
            mutation_rate = 0.01
    
        for N in [20]:
            
            replicates = range(62474,62574)
                
            for replicate in replicates:
                
                rep_data = large_data[large_data["Replicate"] == replicate]
                rep_data_size = list(rep_data["GenomeSize"])
                rep_data_mut = list(rep_data["LethalMutationRate"])
                rep_data_fitness = list(rep_data["Fitness"])
                
                for i,size in enumerate(rep_data_size):
                    
                    if i == 0 or rep_data_fitness[i] == 0 or rep_data_fitness[i-1] ==  0:
                        continue
                    
                    change_size = size-rep_data_size[i-1]
                    
                    g.write(str(change_size)+","+str(rep_data_mut[i]-rep_data_mut[i-1])+"\n")

g.close()