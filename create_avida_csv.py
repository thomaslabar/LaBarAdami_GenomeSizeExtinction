import copy
import os

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

g=open("data_avida.csv","wb")
g.write("Treatment,PopulationSize,MutationRate,Revert,Replicate,Extinction,TimeToExtinction,FinalGenomeSize,FinalFitness,FinalGestationTime,FinalProbLethal,FinalMuLethal\n")

for t in ["varlength_slip","fixedlength"]:
        
    for mu in ["001","01","001_revertlethal","001_revertdel"]:
        
        if mu == "001":
            mutation_rate = 0.01
            revert_status = "None"
            popsize_list = [5,6,7,8,10,15,20]
        elif mu == "01":
            mutation_rate = 0.1
            revert_status = "None"
            popsize_list = [10,12,15,16,17,20,25]
        elif mu == "001_revertlethal":
            mutation_rate = 0.01
            revert_status = "Lethal"
            popsize_list = [5,6,7,8,10,15,20]
        elif mu == "001_revertdel":
            mutation_rate = 0.01
            revert_status = "Del"
            popsize_list = [5,6,7,8,10,15,20]
        
        for N in popsize_list:
        
            final_dom_data = GetFileLines("data_"+t+"_N"+str(N)+"_mu"+mu+"/finaldominant_analysis.dat")
            for line in final_dom_data:
                
                if len(line) > 0 and not line.startswith("#"):
                    
                    temp = line.split(" ")
                    
                    replicate = int(temp[0])

                    avg_data = GetFileLines("data_"+t+"_N"+str(N)+"_mu"+mu+"/replicate_"+str(replicate)+"/average.dat")
    
                    final_avg_line = avg_data[-1].split(" ")
                    final_generation = str(final_avg_line[12])
                    
                    final_size = str(temp[5])
                    final_fitness = str(temp[2])
                    final_gestation = str(temp[4])
                    final_prob_lethal = str(temp[9])
                    final_mu_lethal = str(mutation_rate*float(final_size)*float(final_prob_lethal))
                    
                    extinction = 1
                    if float(final_generation) >= 100000:
                        extinction = 0      
                    
                    g.write(t+","+str(N)+","+str(mutation_rate)+","+str(revert_status)+","+str(replicate)+","+str(extinction)+","+str(final_generation))
                    g.write(","+str(final_size)+","+str(final_fitness)+","+str(final_gestation)+","+str(final_prob_lethal)+","+str(final_mu_lethal)+"\n")

g.close()