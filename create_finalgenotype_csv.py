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

g=open("data_finalgenotypes.csv","wb")
g.write("Treatment,PopulationSize,MutationRate,Revert,Replicate,Sequence\n")

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
            
            
            replicates = [int(i[10:]) for i in os.listdir("data_"+t+"_N"+str(N)+"_mu"+mu) if i != "finaldominant_analysis.dat"]
            replicates.sort()
            
            assert len(replicates) == 100, "Treatment = "+t+" "+str(N)+" "+str(mu)+"; Replicates = "+str(len(replicates))
            assert replicates[-1] == replicates[0] + 99
        
            for replicate in replicates:
                
                final_dom_genotype_data = GetFileLines("data_"+t+"_N"+str(N)+"_mu"+mu+"/replicate_"+str(replicate)+"/final_dominant_genotype.dat")
                
                instructions = []
                for line in final_dom_genotype_data:
                    if len(line) > 0 and not line.startswith("#"):
                        instructions.append(str(line))
                        
                instruction_code = {"nop-A":"a",
                                    "nop-B":"b",
                                    "nop-C":"c",
                                    "if-n-equ":"d",
                                    "if-less":"e",
                                    "if-label":"f",
                                    "mov-head":"g",
                                    "jmp-head":"h",
                                    "get-head":"i",
                                    "set-flow":"j",
                                    "shift-r":"k",
                                    "shift-l":"l",
                                    "inc":"m",
                                    "dec":"n",
                                    "push":"o",
                                    "pop":"p",
                                    "swap-stk":"q",
                                    "swap":"r",
                                    "add":"s",
                                    "sub":"t",
                                    "nand":"u",
                                    "h-copy":"v",
                                    "h-alloc":"w",
                                    "h-divide":"x",
                                    "IO":"y",
                                    "h-search":"z"}
                
                sequence = ""
                for i in instructions:
                    sequence += str(instruction_code[i])
                
                g.write(t+","+str(N)+","+str(mutation_rate)+","+str(revert_status)+","+str(replicate)+","+str(sequence)+"\n")

g.close()