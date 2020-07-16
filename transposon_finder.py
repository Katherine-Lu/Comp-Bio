import matplotlib
import matplotlib.pyplot as plt
import numpy as np

#The purpose of transposon_finder is to compare two DNA sequences, chromosomal and 
# (cont.)virulence plasmid DNA in dysentery bacteria Shigella flexneri strain 2a 2457T, 
#(cont.) to identify conserved sequences in hopes of finding transposons.
#Transposons are DNA sequences that move and attach to other places.
#In Shigella, transposons move from the chromosome to plasmids as a method of gene 
# (cont.) regulation and sometimes confer antibiotic resistance. The plasmid is 220k bp. 
# The chromosome is 4.6 million bp.
#I, Katherine Lu, will only be using the first 70k bp of the chromosome and the first 7k 
# (cont.) bp of the primer because it took 20 min to run only window_compiler using the
# (cont.) whole plasmid and chromosome.
#transposon_finder should output a heatmap showing the order of transposons in the
# (cont.)
#transposon_finder should also output a text file titled "" containing useful information
#about the DNA sequences, transposon statistics, and the order which transposons appear,
#and the transposons themselves.
#transposon_finder still takes 20 min to run. I recommend looking at transposon_test
# (cont.) and transposon_montecarlo first.
#By changing the strings below, transposon_finder can be used to compare other sequences.


export_title = "shigella_flexneri_conserved_sequences.txt"
#title of output text file
seq_type = "DNA"
seq_unit = "bp"
#type of molecule compared
database_1 = "NCBI"
database_2 = "NCBI"
#database these gene sequences are from
test_seq_1_fas = "Shigella_flexneri_chromosome.fasta"
test_seq_2_fas = "shigella_flexneri_plasmid.fasta"
#names of fasta files of sequences compared
test_seq_1_num = "AE014073"
test_seq_2_num = "AF348706.1"
#Some databases have directory numbers.
test_seq_1_loc = "chromosome"
test_seq_2_loc = "pWR501 plasmid"
#Where are your target sequences located?
organism_1 = "Shigella flexneri strain 2a 2457T"
organism_2 = "Shigella flexneri strain 2a 2457T"
#organisms containing the genes being compared


build_seq_1 = open(test_seq_1_fas, "r")
test_seq_1 = ""
for line in build_seq_1:
	if line.count(">") == 0:
		test_seq_1 = test_seq_1 + line
build_seq_1.close()
	
build_seq_2 = open(test_seq_2_fas, "r")
test_seq_2 = ""
for line in build_seq_2:
	if line.count(">") == 0:
		test_seq_2 = test_seq_2 + line
build_seq_2.close()
#For storing fasta files as sequences of nucleotides in transposon_finder.

"""
test_seq_1 = "AAAGGGCCCTTT"
test_seq_2 = "CCCAAATTTGGG"
"""
#To see fast and easy results from transposon_finder, free these simple sequences 
# (cont.) from triple quotes and put everything from lines 46 through  61 in triple quotes


def window_compiler(test_seq_1, test_seq_2):
#window_compiler returns a list of lists of indices that are likely to contain transposons
	transposon_sim = 1 
	#Before running transposon_test or transposon_montecarlo, set transposon_sim = 1
	#Before running transposon_finder, set transposon sim = transposon_avg or transposon_max
	start_points = int(len(test_seq_2)/transposon_sim) 
	# Change transposon_sim if you need.
	# See transposon_montecarlo.py
	window_pile = []
	window_shelf = []
	for test_point in range(0, start_points):
		y = int(transposon_sim/2)
		x = transposon_sim*test_point
		#x is the "center" that y increases around for each iteration
		for y in range(0, len(test_seq_2)):
			if x - y > 0:
				start = x - y 
			else:
				start = 0
			if x + y > len(test_seq_2):
				end = len(test_seq_2)
			else:
				end = x + y
			window_range = []
			if start == 0:
				window_range.append(0)
			if start != 0:
				window_range.append(start+1)
			if end == len(test_seq_2):
				window_range.append(end)
			if end < len(test_seq_2):
				window_range.append(end-1)
			repeat_window = window_pile.count(window_range)
			unclear_window = test_seq_1.count(test_seq_2[start+1:end-1])
			null_window = len(test_seq_2[(start+1):(end-1)])
			if repeat_window < 1 and unclear_window == 1 and null_window > transposon_sim:
			#transposon_sim is used here to filter out any transposons that are shorter
			# (cont.) than randomly generated conserved sequences that still might appear.
				window_pile.append(window_range)
			window_shelf = sorted(window_pile)
	return(window_shelf)	
	


def window_streamliner(test_seq_1, test_seq_2):
#window_streamliner merges windows with overlapping ranges
	window_shelf = window_compiler(test_seq_1, test_seq_2)
	window_library = [window_shelf[0]]
	for window in window_shelf:
		previous = window_library[-1]
		if window[0] < previous[1]:
			previous[1] = max(previous[1], window[1])
		else:
			window_library.append(window)
	return(window_library)

 	

def transposon_identifier(test_seq_1, test_seq_2):
#transposon_identifier returns a list of sequences that are likely transposons
	window_library = window_streamliner(test_seq_1, test_seq_2)
	transposon_library = []
	for window in window_library:
		transposon_seq = test_seq_2[window[0]:window[1]]
		transposon_library.append(transposon_seq)
	return(transposon_library)








def seq_2_order(test_seq_1, test_seq_2):
#plasmid_order returns a dict assigning transposons an integer in the order they
# (cont.) appear in test_seq_2
	transposon_index_2 = {}
	transposon_library = transposon_identifier(test_seq_1, test_seq_2)
	order = 0
	for transposon in transposon_library:
		transposon_index_2[transposon] = order + 1
		order = order + 1
	return(transposon_index_2)
	

		
def seq_1_order(test_seq_1, test_seq_2):
#seq_1_order returns a list of the order of transposons in test_seq_1 relative to their
# (cont.) order in test_seq_2
	transposon_index_2 = seq_2_order(test_seq_1, test_seq_2)
	transposon_library = transposon_identifier(test_seq_1, test_seq_2)
	order_pile = []
	for transposon in transposon_library:
		local_index = []
		local_index.append(test_seq_1.index(transposon))
		local_index.append(transposon)
		order_pile.append(local_index)	
	transposon_index_1 = []
	new_order = sorted(order_pile)
	for rank in new_order:
		transposon_index_1.append(transposon_index_2[rank[1]])
	return(transposon_index_1)
	
	
	
		




def export_explainer(test_seq_1, test_seq_2):
#export_explainer outputs a text file that explains what's going on.
	info_fill = open(export_title, "a")
	info_fill.write("Conserved sequences and potential transposons" + "\n")
	if organism_1 == organism_2:
		info_fill.write("Comparing" + " " + test_seq_1_loc + " " + "and" + " " + test_seq_2_loc
		+ " " + seq_type + " " + "from" + " " + organism_1 + "\n")
	else:
		info_fill.write("Comparing" + " " + test_seq_1_loc + " " + seq_type + " " + "from"
		+ " " + organism_1 + " " + "and" + " " + test_seq_2_loc + " " + seq_type + "\n")
	if database_1 == database_2:
		info_fill.write(seq_type + " " + "sequences from" + " " + database_1 + "\n")
	else:
		info_fill.write(seq_type + " " + "sequence for" + " " + test_seq_1_loc + " "
	    + "from" + " " + database_1 + "\n")
		info_fill.write(seq_type + " " + "sequence for" + " " + test_seq_2_loc + " " +
		"from" + " " + database_2 + "\n")
	info_fill.write(test_seq_1_loc + " " + "directory number:" + " " + test_seq_1_num + "\n")
	info_fill.write(test_seq_2_loc + " " + "directory number:" + " " + test_seq_2_num + "\n")
	info_fill.close()



def transposon_stats_exporter(test_seq_1, test_seq_2):
#transposon_stats_exporter adds statistics about the transposons found to the text file
#(cont.) created by export_explainer
	transposon_library = transposon_identifier(test_seq_1, test_seq_2)
	transposon_sizes = []
	for transposon in transposon_library:
		transposon_sizes.append(len(transposon))
	transposon_array = np.array(transposon_sizes)
	transposon_avg = np.mean(transposon_array, axis = 0)
	transposon_max = np.amax(transposon_array, axis = 0)
	transposon_stdev = np.std(transposon_array, axis = 0, ddof = len(transposon_library)-1)
	transposon_pct = (len(transposon_library)*transposon_avg/len(test_seq_2))*100
	stat_fill = open(export_title, "a")
	stat_fill.write("Average conserved sequence length:" + " " + str(transposon_avg) + " "
	+ seq_unit + "\n")
	stat_fill.write("Longest conserved sequence length:" + " " + str(transposon_max) + " " 
	+ seq_unit + "\n")
	stat_fill.write("Conserved sequence standard deviation:" + " " + str(transposon_stdev) + "\n")
	stat_fill.write("Total conserved sequences in" + " " + test_seq_2_loc + ":" + "\n" +
	str(transposon_pct) + "%" + "\n")
	stat_fill.close()
	
	
	
def order_exporter(test_seq_1, test_seq_2):
#order_exporter adds the order of transposons in test_seq_1 to the text file created by
# (cont.) export_explainer
	transposon_order = seq_1_order(test_seq_1, test_seq_2)
	order_fill = open(export_title, "a")
	order_fill.write("Order of conserved elements in" + " " + test_seq_1_loc + " " +
	"relative to order in" + " " + test_seq_2_loc + ":" + "\n")
	order_fill.write(str(transposon_order) + "\n")
	order_fill.close()


		
def transposon_exporter(test_seq_1, test_seq_2):
#transposon_exporter adds the stars of the show, transposons, to the text file created by
# (cont.) export_explainer
	transposon_library = transposon_identifier(test_seq_1, test_seq_2)
	transposon_fill = open(export_title, "a")
	transposon_fill.write("Conserved sequences in order they appear in" 
	+ " " + test_seq_2_loc + ":" + "\n")
	for transposon in transposon_library:
		transposon_fill.write(transposon + "\n")
	transposon_fill.close()
	
	
	
	
	
	
	
def heatmapper(test_seq_1, test_seq_2):
#heatmapper generates a labeled and titled heatmap of relative order of transposons
	index_2_integer = []
	transposon_index_1 = seq_1_order(test_seq_1, test_seq_2)
	window_library = window_streamliner(test_seq_1, test_seq_2)
	i = 1
	for rank in range(len(window_library)):
		index_2_integer.append(i)
		i = i + 1
	index_pile = []
	index_pile.append(index_2_integer)
	index_pile.append(transposon_index_1)
	index_array = np.array(index_pile)
	print(index_array)
	fig, ax = plt.subplots()
	im = ax.imshow(index_array)
	ax.set_xticks(np.arange(len(window_library)))
	ax.set_yticks(np.arange(2))
	ax.set_yticklabels([test_seq_2_loc, test_seq_1_loc])
	if organism_1 == organism_2:
		ax.set_title("Conserved sequence order between" + " " + organism_1 + " " + test_seq_2_loc 
		+ " " + "and" + " " + test_seq_1_loc + " " + seq_type)
	else:
		ax.set_title("Conserved sequence order between" + " " + organism_1 + " " + test_seq_2_loc 
		+ " " + "and" + " " + organism_2 + " " + test_seq_1_loc + " " + seq_type)
	fig.tight_layout()
	plt.show()
	
	
	
export_explainer(test_seq_1, test_seq_2)
transposon_stats_exporter(test_seq_1, test_seq_2)
order_exporter(test_seq_1, test_seq_2)
transposon_exporter(test_seq_1, test_seq_2)
heatmapper(test_seq_1, test_seq_2)
#Time to run functions. There are fewer functions here than written above because 
#(cont.) many of them are called within functions.


#Issues: 
#The "widen window around starting point" method used by window_compiler is not very 
# (cont.)accurate because the range of each window can only extend outward based on the
#(cont.) 'centers' generated by transposon_montecarlo.
#This could be worked around by running transposon_finder on the same sequences using
#(cont.) different transposon_sim each time.
#transposon_finder can work on any type of sequence, but transposon_montecarlo only
#(cont.) generates DNA sequences. Less accessable than I thought.
#I made the export function series because the heatmap becomes unreadable with longer 
# (cont.) test sequences. Too wide and there are so many transposons the colors are bad.
#This also takes a very long time to run.
#How to optimise?