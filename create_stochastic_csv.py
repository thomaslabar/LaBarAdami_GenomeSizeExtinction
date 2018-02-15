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

data_avida = pd.read_csv("data_avida.csv")

data_avida_varlength_slip = data_avida[data_avida["Treatment"] == "varlength_slip"]

data_avida_varlength_slip_N5 = data_avida_varlength_slip[data_avida_varlength_slip["PopulationSize"] == 5]
data_avida_varlength_slip_N5_mu001 = data_avida_varlength_slip_N5[data_avida_varlength_slip_N5["MutationRate"] == 0.01]
data_avida_varlength_slip_N5_mu001_revertlethal = data_avida_varlength_slip_N5_mu001[data_avida_varlength_slip_N5_mu001["Revert"] == "Lethal"]

data_avida_varlength_slip_N8 = data_avida_varlength_slip[data_avida_varlength_slip["PopulationSize"] == 8]
data_avida_varlength_slip_N8_mu001 = data_avida_varlength_slip_N8[data_avida_varlength_slip_N8["MutationRate"] == 0.01]
data_avida_varlength_slip_N8_mu001_revertnone = data_avida_varlength_slip_N8_mu001[data_avida_varlength_slip_N8_mu001["Revert"] == "None"]

extinct_replicates = {}
extinct_replicates["revertnone"] = list(data_avida_varlength_slip_N8_mu001_revertnone[data_avida_varlength_slip_N8_mu001_revertnone["Extinction"] == 1]["Replicate"])
extinct_replicates["revertlethal"] = list(data_avida_varlength_slip_N5_mu001_revertlethal[data_avida_varlength_slip_N5_mu001_revertlethal["Extinction"] == 1]["Replicate"])

g=open("data_stochastic.csv","wb")
g.write("RevertStatus,Replicate,GenomeSize,Extinction,PercentNonViableAnalyses\n")

for revert_status in ["revertnone","revertlethal"]:
    
    stochastic_data =  GetFileLines("data_stochastic/finaldominant_stochastic_analysis_"+revert_status+".dat")
    
    for line in stochastic_data:
        
        if len(line) > 0 and not line.startswith("#"):
            
            temp = line.split(" ")
            
            replicate = int(temp[0])
            
            percent_nonviable = 1.0-float(temp[1])
            
            if revert_status == "revertnone":
                rep_data = data_avida_varlength_slip_N8_mu001_revertnone[data_avida_varlength_slip_N8_mu001_revertnone["Replicate"] == replicate]
            elif revert_status == "revertlethal":
                rep_data = data_avida_varlength_slip_N5_mu001_revertlethal[data_avida_varlength_slip_N5_mu001_revertlethal["Replicate"] == replicate]
                
            assert len(rep_data) == 1
            
            genome_size = list(rep_data["FinalGenomeSize"])[0]
            extinction = list(rep_data["Extinction"])[0]
                
                
            g.write(str(revert_status)+","+str(replicate)+","+str(genome_size)+","+str(extinction))
            g.write(","+str(percent_nonviable)+"\n")
        
g.close()