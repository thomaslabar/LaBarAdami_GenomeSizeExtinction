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

g=open("data_largegenomes.csv","wb")
g.write("Treatment,PopulationSize,MutationRate,Replicate,Entry,Fitness,GenomeSize,LethalMutationRate\n")

for t in ["varlength_slip"]:
    
    for mu in ["001"]:
        
        if mu == "001":
            mutation_rate = 0.01
    
        for N in [20]:
            
            replicates = range(62474,62574)
                
            for replicate in replicates:
                
                mut_lod_data = GetFileLines("data_"+t+"_N"+str(N)+"_mu"+mu+"/replicate_"+str(replicate)+"/mut_lod.dat")
                
                ct = 0
                for i,entry in enumerate(mut_lod_data):
                    
                    if len(entry) > 0 and not entry.startswith("#"):
                        
                       
                        temp = entry.split(" ")
                        size = int(temp[1])
                        lethal_mu = float(temp[6])*mutation_rate*size
                        fitness = str(temp[2])
                    
                        g.write(t+","+str(N)+","+str(mutation_rate)+","+str(replicate)+","+str(ct)+","+str(fitness)+","+str(size))
                        g.write(","+str(lethal_mu)+"\n")
                        ct += 1

g.close()