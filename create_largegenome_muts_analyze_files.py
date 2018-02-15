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

replicates = [int(i[10:]) for i in os.listdir("data_varlength_slip_N20_mu001") if i != "finaldominant_analysis.dat"]
replicates.sort()

f = open("config_files/analyze_varlength_slip_N20_mu001_largegenomes_muts.cfg","w")
f.write("VERBOSE\n")

for rep in replicates:

    lod_data = GetFileLines("data_varlength_slip_N20_mu001/replicate_"+str(rep)+"/lod.dat")
    
    size_change_seqs = []
    
    for i,line in enumerate(lod_data):
        
        if len(line) > 0 and not line.startswith("#"):
            
            temp = line.split(" ")
            if temp[2] == "wzcagczvfcaxgab":
                size_change_seqs.append(str(temp[2]))
                continue
            
            temp_prev = lod_data[i-1].split(" ")
    
            size = int(temp[1])
            prev_size = int(temp_prev[1])
            if size != prev_size:
                if str(temp_prev[2]) != size_change_seqs[-1]:
                    size_change_seqs.append(str(temp_prev[2]))
                size_change_seqs.append(str(temp[2]))
    
    f.write("PURGE_BATCH\n")
    for seq in size_change_seqs:
        f.write("LOAD_SEQUENCE "+str(seq)+"\n")
    f.write("RECALCULATE\n")
    f.write("DETAIL data_varlength_slip_N20_mu001/replicate_"+str(rep)+"/mut_lod.dat id length fitness frac_pos frac_neut frac_neg frac_dead sequence\n")
f.close()

g = open("sub_files/avida_a_varlength_slip_N20_mu001_largegenomes_muts.sub","w")
g.write("#!/bin/bash -login\n")
g.write("#PBS -l walltime=24:00:00,nodes=1:ppn=1,mem=2gb\n")
g.write("#PBS -j oe\n")
g.write("cd ${HOME}/Documents/genomesize_extinction/avida/config_files\n")
g.write("./avida -a -set ANALYZE_FILE analyze_varlength_slip_N20_mu001_largegenomes_muts.cfg -set DATA_DIR .. -set WORLD_X 8 ")
g.write(" -set WORLD_Y 1 -set BIRTH_METHOD 4 -set COPY_MUT_PROB 0.0 -set DIV_MUT_PROB 0.01")
g.write(" -set DIVIDE_INS_PROB 0.0 -set DIVIDE_DEL_PROB 0.0 -set DIVIDE_SLIP_PROB 0.0 -set REQUIRE_EXACT_COPY 1")
g.close()
