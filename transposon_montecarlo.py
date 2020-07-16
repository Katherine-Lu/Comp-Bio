#Before running transposon_montecarlo,in transposon_finder line 72, set transposon_sim = 1
#transposon_montecarlo should output stats of randomly generated transposons in the
# (cont.) command line.
#Use these results to set transposon_sim to 'average transposon length' or 
# (cont.) 'maximum transposon length' before running transposon_finder to filter out
# (cont.) any conserved sequences that may appear by chance
from random import choices
import numpy as np
from transposon_finder import window_compiler, window_streamliner, transposon_identifier
nucleotides = ["A", "C", "T", "G"]
nt_placer_1 = choices(nucleotides, weights = [1, 1, 1, 1], k = 1200)
nt_placer_2 = choices(nucleotides, weights = [1, 1, 1, 1], k = 900)
#k can be changed depending on how large your test sequences are. 



def pseudo_chromosome(nt_placer_1):
	seq_joiner_1 = ""
	test_seq_1 = seq_joiner_1.join(nt_placer_1)
	return(test_seq_1)



def pseudo_plasmid(nt_placer_2):
	seq_joiner_2 = ""
	test_seq_2 = seq_joiner_2.join(nt_placer_2)
	return(test_seq_2)
	
	
	
def pseudo_transposon_stats(test_seq_1, test_seq_2):
	test_seq_2 = pseudo_plasmid(nt_placer_2)
	test_seq_1 = pseudo_chromosome(nt_placer_1)
	transposon_library = transposon_identifier(test_seq_1, test_seq_2)
	transposon_sizes = []
	for transposon in transposon_library:
		transposon_sizes.append(len(transposon))
	transposon_array = np.array(transposon_sizes)
	transposon_avg = np.mean(transposon_array, axis = 0)
	transposon_max = np.amax(transposon_array, axis = 0)
	transposon_stdev = np.std(transposon_array, axis = 0, ddof = len(transposon_library)-1)
	print("average transposon length:", transposon_avg)
	print("transposon length standard deviation:", transposon_stdev)
	print("maximum transposon length:", transposon_max)
	print("number of transposons:", len(transposon_library))

	
	
	
pseudo_transposon_stats(nt_placer_1, nt_placer_2)
	


	