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

f = open("config_files/analyze_varlength_slip_N20_mu001_largegenomes.cfg","w")
f.write("VERBOSE\n")

for rep in replicates: 
    rep_files = [i for i in os.listdir("data_varlength_slip_N20_mu001/replicate_"+str(rep)) if "detail" in i]
    assert len(rep_files) == 1
    detail_file = str(rep_files[0])
    f.write("PURGE_BATCH\n")
    f.write("LOAD  ../data_varlength_slip_N20_mu001/replicate_"+str(rep)+"/"+detail_file+"\n")
    f.write("FIND_LINEAGE num_cpus\n")
    f.write("DETAIL data_varlength_slip_N20_mu001/replicate_"+str(rep)+"/lod.dat id length sequence\n")
f.close()

g = open("sub_files/avida_a_varlength_slip_N20_mu001_largegenomes.sub","w")
g.write("#!/bin/bash -login\n")
g.write("#PBS -l walltime=3:59:00,nodes=1:ppn=1,mem=2gb\n")
g.write("#PBS -j oe\n")
g.write("cd ${HOME}/Documents/genomesize_extinction/avida/config_files\n")
g.write("./avida -a -set ANALYZE_FILE analyze_varlength_slip_N20_mu001_largegenomes.cfg -set DATA_DIR .. -set WORLD_X 8 ")
g.write(" -set WORLD_Y 1 -set BIRTH_METHOD 4 -set COPY_MUT_PROB 0.0 -set DIV_MUT_PROB 0.01")
g.write(" -set DIVIDE_INS_PROB 0.0 -set DIVIDE_DEL_PROB 0.0 -set DIVIDE_SLIP_PROB 0.0 -set REQUIRE_EXACT_COPY 1")
g.close()
