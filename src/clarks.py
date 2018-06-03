import numpy as np
import csv
import pdb
import itertools as it
from progress.bar import Bar

# Load data as 2d int array
def load(file_name):
	data = []
	with open(file_name, newline="\n") as csvfile:
		reader = csv.reader(csvfile, delimiter=' ')
		for row in reader:
			data.append([int(val) for val in row])
	return data

# Phase haplotypes from unambiguous genotypes
# TODO: assign 0 or 1 based on observed frequency
def gen_2_haps(genotype):
	hap1, hap2 = [], []
	for gen in genotype:
		if gen == 0:
			hap1.append(0)
			hap2.append(0)
		if gen == 2:
			hap1.append(1)
			hap2.append(1)
		if gen == 1:
			hap1.append(0)
			hap2.append(1) 
	return hap1, hap2

# Return 2d array of all known haplotypes
def get_known_haps(haplotypes,window_len):
	num_snps = len(haplotypes)
	num_haps = len(haplotypes[0])
	known_haps = []
	for i in range(num_haps):
		haplotype = [row[i] for row in haplotypes]
		for j in range(0, num_snps, window_len):
			haplotype_segment = haplotype[j:j+window_len]
			if haplotype_segment.count(-1) == 0 and haplotype_segment not in known_haps:
				known_haps.append(haplotype_segment)
	return known_haps

# Phase all unambiguous haplotypes of a given length
def fill_known_haps(data,window_len):
	num_snps = len(data)
	num_indiv = len(data[0])
	#print("num snps, indv = {} {}".format(num_snps,num_indiv))
	haplotypes = np.zeros((num_snps,2*num_indiv),dtype=np.int)
	haplotypes.fill(-1)
	known_haps = []
	for i in range(num_indiv):
		genotype = [row[i] for row in data]
		for j in range(0,num_snps,window_len):
			genotype_segment = genotype[j:j+window_len]
			if genotype_segment.count(1) < 2:
				hap1, hap2 = gen_2_haps(genotype_segment)
				if hap1 not in known_haps:
					known_haps.append(hap1)
				if hap2 not in known_haps:
					known_haps.append(hap2)
				hap1_np = np.array(hap1)
				hap2_np = np.array(hap2)
				for k in range(window_len):
					haplotypes[j+k][2*i] = hap1[k]
					haplotypes[j+k][2*i+1] = hap2[k]
	return haplotypes, known_haps

# Helper function for adding two haplotypes
def add_haplotypes(hap1, hap2):
    add_result = [hap1[i] +hap2[i] for i in range(len(hap1))]
    return add_result

# Helper function for subtracting two haplotypes
def subtract_genhap(genotype,haplotype):
	sub_result = [genotype[i] - haplotype[i] for i in range(len(genotype))]
	return sub_result

# Makes a dictionary of all possible pair-wise combinations of known haplotypes
# TODO: Remove incompatible haplotypes from known set before calling func in Clark's
#		to reduce # of comparisons.
def make_hap_map(known_set):
    hap_map = {}
    all_combos = list(it.combinations(known_set, 2))
    for pair in all_combos:
        hap_map[str(add_haplotypes(pair[0], pair[1]))] = [pair[0], pair[1]]
    return hap_map

# Infer haplotype 2 if genotype and haplotype 1 make it unambiguous 
def get_hap2_from_known_hap1(genotype, known_set):
    def validHap(haplotype):
        valid = True
        for element in haplotype:
            if element is not 0 and element is not 1:
                valid = False
        return valid
    for h1 in known_set:
        h2 = subtract_genhap(genotype, h1)
        if validHap(h2):
            return h1, h2
    return [],[]

# Guess rest of haplotypes that Clark's algorithm couldn't phase
# TODO: Incorporate frequency so that we don't randomly guess
def guess_rest(data,haplotypes):
	num_snps = len(data)
	num_indiv = len(data[0])
	for i in range(num_indiv):
		for j in range(num_snps):
			if haplotypes[j][2*i] == -1 and haplotypes[j][2*i+1] == -1:
				if data[j][i] == 0:
					haplotypes[j][2*i], haplotypes[j][2*i+1] = 0,0
				elif data[j][i] == 1:
					haplotypes[j][2*i], haplotypes[j][2*i+1] = 0,1
				elif data[j][i] == 2:
					haplotypes[j][2*i], haplotypes[j][2*i+1] = 1,1
			elif haplotypes[j][2*i] == -1:
				if haplotypes[j][2*i+1] == 0:
					if data[j][i] == 0:
						haplotypes[j][2*i] = 0
					elif data[j][i] == 1:
						haplotypes[j][2*i] = 1
					elif data[j][i] == 2:
						#should never be called
						print("SOMETHINGS WRONG")
						haplotypes[j][2*i] = 1
				elif haplotypes[j][2*i+1] == 1:
					if data[j][i] == 0:
						#should never be called
						print("SOMETHINGS WRONG")
						haplotypes[j][2*i] = 0
					elif data[j][i] == 1:
						haplotypes[j][2*i] = 0
					elif data[j][i] == 2:
						haplotypes[j][2*i] = 1
			elif haplotypes[j][2*i+1] == -1:
				if haplotypes[j][2*i] == 0:
					if data[j][i] == 0:
						haplotypes[j][2*i+1] = 0
					elif data[j][i] == 1:
						haplotypes[j][2*i+1] = 1
					elif data[j][i] == 2:
						#should never be called
						print("SOMETHINGS WRONG")
						haplotypes[j][2*i+1] = 1
				elif haplotypes[j][2*i] == 1:
					if data[j][i] == 0:
						#should never be called
						print("SOMETHINGS WRONG")
						haplotypes[j][2*i+1] = 0
					elif data[j][i] == 1:
						haplotypes[j][2*i+1] = 0
					elif data[j][i] == 2:
						haplotypes[j][2*i+1] = 1
	return haplotypes

# Main algorithm to phase genotype data
def clarks(data, window_len):
	num_snps = len(data)
	num_indiv = len(data[0])
	haplotypes, known_haps = fill_known_haps(data,window_len)

	hap_map = make_hap_map(known_haps)
	for num_iter in range(10):
		curr_known_haps_size = len(known_haps)
		#bar_i = Bar('Phasing Individuals...', max=num_indiv)		
		for i in range(num_indiv):
			genotype = [row[i] for row in data]
			#bar_j = Bar('SNPs', max=(int)(num_snps/window_len))
			for j in range(0,num_snps,window_len):
				#print("Iteration {}, Indiv {}, Window range [{},{}]".format(num_iter,i,j,j+window_len))
				genotype_segment = genotype[j:j+window_len] 
				seg_key = str(genotype_segment)

				if seg_key in hap_map:
					hap1, hap2 = hap_map[seg_key]
					hap1_np = np.array(hap1)
					hap2_np = np.array(hap2)
					for k in range(window_len):
						haplotypes[j+k][2*i] = hap1[k]
						haplotypes[j+k][2*i+1] = hap2[k]
				else:
					hap1, hap2 = get_hap2_from_known_hap1(genotype_segment,known_haps)
					if len(hap1) > 0 and len(hap2) > 0:
						hap1_np = np.array(hap1)
						hap2_np = np.array(hap2)
						for k in range(window_len):
							haplotypes[j+k][2*i] = hap1[k]
							haplotypes[j+k][2*i+1] = hap2[k]
						known_haps.append(hap2)
						hap_map = make_hap_map(known_haps)
				#bar_j.next()
			#bar_j.finish()
			#bar_i.next()
		#decrease window size if we didn't add any to our known set
		if len(known_haps) - curr_known_haps_size == 0:
			if window_len == 30:
				window_len = 20
			elif window_len == 20:
				window_len = 15
			elif window_len == 15:
				window_len = 10
			elif window_len == 10:
				break
			known_haps = get_known_haps(haplotypes,window_len)
		#bar_i.finish()

	#print("# of -1 before guess =", np.count_nonzero(haplotypes==-1))
	final_haplotypes = guess_rest(data,haplotypes)
	#print("# of -1 before guess =", np.count_nonzero(final_haplotypes==-1))

	return final_haplotypes

# Divides up data into blocks and runs Clark's algorithm one block at at time.
if __name__ == "__main__":
	data = load("../data/example_data_1/example_data_1_180.txt")

	num_snps = len(data)
	num_indivs = len(data[0])

	# Hyperparameters for splitting up the dataset
	# TODO: Optimize and choose best values
	block_size = 180
	window_len = 30

	haplotypes = []

	if len(data) % block_size == 0:
		num_blocks = (int)(len(data)/block_size)
	else:
		num_blocks = (int)(len(data)/block_size) + 1
	#print((int)(len(data)/block_size) + math.ceil((num_snps % block_size) / block_size)

	print("********** Clark's Algorithm **********")
	print("# of SNP sites =", num_snps)
	print("# of Individuals =", num_indivs)
	print("Block Size =", block_size)
	print("Window Length =", window_len)
	print("# of Blocks =", num_blocks,"\n")
	print("Starting...")
	
	bar = Bar('Progress Bar', max=num_blocks,suffix = '%(percent).1f%% - %(eta)ds')
	for i in range(num_blocks):
		#print("Working on block #{} out of {}".format(i+1, num_blocks))
		data_block = data[i*block_size:i*block_size+block_size]
		if i == 0:
			haplotypes_block = clarks(data_block,window_len)
			haplotypes = haplotypes_block
		elif i == num_blocks-1 and len(data[i*block_size:]) > 0:
			final_data_block = data[i*block_size:]
			final_num_snps = len(final_data_block)
			final_num_indiv = len(final_data_block[0])
			final_haplotypes_block = np.zeros((final_num_snps,2*final_num_indiv),dtype=np.int)
			final_haplotypes_block.fill(-1)
			final_haplotypes_block = guess_rest(final_data_block, final_haplotypes_block)
			haplotypes = np.concatenate((haplotypes,final_haplotypes_block),axis=0)
		else:
			haplotypes_block = clarks(data_block,window_len)
			haplotypes = np.concatenate((haplotypes,haplotypes_block),axis=0)
		bar.next()
	bar.finish()
	print("********** Clark's Algorithm **********")

	#print("Final haps row, col = {} {}".format(len(haplotypes), len(haplotypes[0])))
	np.savetxt('../data/example_data_1/example_data_1_180_my_sol.txt', haplotypes, fmt='%i', delimiter = ' ')


	










	