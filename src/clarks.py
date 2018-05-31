import numpy as np
import csv
import pdb
import itertools as it
from progress.bar import Bar

# LOAD FILE AS INT ARRAY
def load(file_name):
	data = []
	with open(file_name, newline="\n") as csvfile:
		reader = csv.reader(csvfile, delimiter=' ')
		for row in reader:
			data.append([int(val) for val in row])
	return data

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

def fill_known_haps(data,window_len):
	num_snps = len(data)
	num_indiv = len(data[0])
	#print("num snps, indv = {} {}".format(num_snps,num_indiv))
	haplotypes = np.zeros((num_snps,2*num_indiv))
	haplotypes.fill(1)

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
					haplotypes[j+k][i] = hap1[k]
					haplotypes[j+k][2*i+1] = hap2[k]

	return haplotypes, known_haps


def add_haplotypes(hap1, hap2):
    add_result = [hap1[i] +hap2[i] for i in range(len(hap1))]
    return add_result

def subtract_genhap(genotype,haplotype):
	sub_result = [genotype[i] - haplotype[i] for i in range(len(genotype))]
	return sub_result

def make_hap_map(haplotype_list):
    final_map = {}
    all_combos = list(it.combinations(haplotype_list, 2))
    for pair in all_combos:
        final_map[str(add_haplotypes(pair[0], pair[1]))] = [pair[0], pair[1]]
    return final_map

def get_hap2_from_known_hap1(genotype, bank):
    def validHap(haplotype):
        value = True
        for element in haplotype:
            if element is not 0 and element is not 1:
                value = False
        return value

    for hap in bank:
        new_haplotype = subtract_genhap(genotype, hap)
        if validHap(new_haplotype):
            return hap, new_haplotype
    
    return [],[]

def clarks(data, window_len):
	num_snps = len(data)
	num_indiv = len(data[0])
	
	haplotypes, known_haps = fill_known_haps(data,window_len)
	print("********** Clark's Algorithm **********")
	hap_map = make_hap_map(known_haps)
	for num_iter in range(10):
		print("Iteration #", num_iter)
		curr_known_haps_size = len(known_haps)
		bar_i = Bar('Phasing Individuals...', max=num_indiv)		
		for i in range(num_indiv):
			genotype = [row[i] for row in data]
			#bar_j = Bar('SNPs', max=(int)(num_snps/window_len))
			for j in range(0,num_snps,window_len):
				genotype_segment = genotype[j:j+window_len] 
				seg_key = str(genotype_segment)

				if seg_key in hap_map:
					hap1, hap2 = hap_map[seg_key]
					hap1_np = np.array(hap1)
					hap2_np = np.array(hap2)
					for k in range(window_len):
						haplotypes[j+k][i] = hap1[k]
						haplotypes[j+k][2*i+1] = hap2[k]
				else:
					hap1, hap2 = get_hap2_from_known_hap1(genotype_segment,known_haps)
					if len(hap1) > 0 and len(hap2) > 0:
						hap1_np = np.array(hap1)
						hap2_np = np.array(hap2)
						for k in range(window_len):
							haplotypes[j+k][i] = hap1[k]
							haplotypes[j+k][2*i+1] = hap2[k]
						known_haps.append(hap2)
						hap_map = make_hap_map(known_haps)
				#bar_j.next()
			#bar_j.finish()
			bar_i.next()
		if len(known_haps) - curr_known_haps_size == 0:
			break
		bar_i.finish()
	print("\n********** Clark's Algorithm **********")
	return haplotypes

if __name__ == "__main__":
	data = load("../data/example_data_1/example_data_1_180.txt")
	#print("data row, col = {} {}".format(len(data), len(data[0])))

	window_len = 10

	haplotypes = clarks(data,window_len)
	#print("window size =", window_len)
	#print("haps row, col = {} {}".format(len(haplotypes), len(haplotypes[0])))


	np.savetxt('../data/example_data_1/example_data_1_180_my_sol.txt',haplotypes, fmt='%i', delimiter = ' ')


	










	