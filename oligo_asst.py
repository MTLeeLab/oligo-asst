"""
This file is part of Oligo-ASST, released under the 3-Clause BSD License.
See License.txt or go to https://opensource.org/licenses/BSD-3-Clause for
full details.

Copyright 2020 Miler T. Lee.


oligo_tile.py: Design antisense oligos to target RNA sequences
Author: Miler Lee
"""

import mtl_fasta

import math
import re
import argparse
import sys


canonical_bases = {
	('A',): 'A',
	('C',): 'C',
	('G',): 'G',
	('T',): 'T'
	}

iupac_wc = {'A': ('A',),
			'C': ('C',),
			'G': ('G',),
			'T': ('T',),
	'R': ('A', 'G'),
	'Y': ('C', 'T'),
	'S': ('C', 'G'),
	'W': ('A', 'T'),
	'K': ('G', 'T'),
	'M': ('A', 'C'),
	'B': ('C', 'G', 'T'),
	'D': ('A', 'G', 'T'),
	'H': ('A', 'C', 'T'),
	'V': ('A', 'C', 'G'),
	'N': ('A', 'C', 'G', 'T')}

iupac_wc_rev = {
	('A',): 'A',
	('C',): 'C',
	('G',): 'G',
	('T',): 'T',
	('A', 'G'): 'R',
	('C', 'T'): 'Y',
	('C', 'G'): 'S',
	('A', 'T'): 'W',
	('G', 'T'): 'K',
	('A', 'C'): 'M',
	('C', 'G', 'T'): 'B',
	('A', 'G', 'T'): 'D',
	('A', 'C', 'T'): 'H',
	('A', 'C', 'G'): 'V',
	('A', 'C', 'G', 'T'): 'N'}

#Source: Sugimoto 1995 (RNA/DNA)
# listed with respect to the RNA target
RNADNAdeltaH = {'AA': 7.8,
		  'AC': 5.9,
		  'AG': 9.1,
		  'AT': 8.3,
		  'CA': 9.0,
		  'CC': 9.3,
		  'CG': 16.3,
		  'CT': 7.0,
		  'GA': 5.5,
		  'GC': 8.0,
		  'GG': 12.8,
		  'GT': 7.8,
		  'TA': 7.8,
		  'TC': 8.6,
		  'TG': 10.4,
		  'TT': 11.5
		  }
		   
RNADNAdeltaS = {'AA': 0.0219,
		  'AC': 0.0123,
		  'AG': 0.0235,
		  'AT': 0.0239,
		  'CA': 0.0261,
		  'CC': 0.0232,
		  'CG': 0.0471,
		  'CT': 0.0197,
		  'GA': 0.0135,
		  'GC': 0.0171,
		  'GG': 0.0319,
		  'GT': 0.0216,
		  'TA': 0.0232,
		  'TC': 0.0229,
		  'TG': 0.0284,
		  'TT': 0.0364
		  }
		  

#Source: Sugimoto 1996 (DNA/DNA)
DNAdeltaH = {'AA': 8.0,
		  'TT': 8.0,
		  'AT': 5.6,
		  'TA': 6.6,
		  'CA': 8.2,
		  'TG': 8.2,
		  'GT': 9.4,
		  'AC': 9.4,
		  'CT': 6.6,
		  'AG': 6.6,
		  'GA': 8.8,
		  'TC': 8.8,
		  'CG': 11.8,
		  'GC': 10.5,
		  'GG': 10.9,
		  'CC': 10.9
		   }

#kcal/K*mol
DNAdeltaS = {'AA': 0.0219,
		  'TT': 0.0219,
		  'AT': 0.0152,
		  'TA': 0.0184,
		  'CA': 0.0210,
		  'TG': 0.0210,
		  'GT': 0.0255,
		  'AC': 0.0255,
		  'CT': 0.0164,
		  'AG': 0.0164,
		  'GA': 0.0235,
		  'TC': 0.0235,
		  'CG': 0.0290,
		  'GC': 0.0264,
		  'GG': 0.0284,
		  'CC': 0.0284
		  }

#Source: SantaLucia 1998
RNAdeltaH = {'AA': 6.82,
		  'AC': 11.4,
		  'AG': 10.48,
		  'AT': 9.38,
		  'CA': 10.44,
		  'CC': 13.39,
		  'CG': 10.64,
		  'CT': 10.48,
		  'GA': 12.44,
		  'GC': 14.88,
		  'GG': 13.39,
		  'GT': 11.4,
		  'TA': 7.69,
		  'TC': 12.44,
		  'TG': 10.44,
		  'TT': 6.82,
		  }
		   
RNAdeltaS = {'AA': 0.019,
		  'AC': 0.0295,
		  'AG': 0.0271,
		  'AT': 0.0267,
		  'CA': 0.0269,
		  'CC': 0.0327,
		  'CG': 0.0267,
		  'CT': 0.0271,
		  'GA': 0.0325,
		  'GC': 0.0369,
		  'GG': 0.0327,
		  'GT': 0.0295,
		  'TA': 0.0205,
		  'TC': 0.0325,
		  'TG': 0.0269,
		  'TT': 0.0190
		  }


RNA_DNA_HELIX_INIT = -3.1 #deltaG
DNA_HELIX_INIT = -3.4
RNA_HELIX_INIT = -4.09
R = 0.001987 #kcal/K*mol #gas constant


###string.maketrans('acgturyACGTURY', 'tgcaayrTGCAAYR')
DNA_RC = '\x00\x01\x02\x03\x04\x05\x06\x07\x08\t\n\x0b\x0c\r\x0e\x0f\x10\x11\x12\x13\x14\x15\x16\x17\x18\x19\x1a\x1b\x1c\x1d\x1e\x1f !"#$%&\'()*+,-./0123456789:;<=>?@TBGDEFCHIJKLMNOPQYSAAVWXRZ[\\]^_`tbgdefchijklmnopqysaavwxrz{|}~\x7f\x80\x81\x82\x83\x84\x85\x86\x87\x88\x89\x8a\x8b\x8c\x8d\x8e\x8f\x90\x91\x92\x93\x94\x95\x96\x97\x98\x99\x9a\x9b\x9c\x9d\x9e\x9f\xa0\xa1\xa2\xa3\xa4\xa5\xa6\xa7\xa8\xa9\xaa\xab\xac\xad\xae\xaf\xb0\xb1\xb2\xb3\xb4\xb5\xb6\xb7\xb8\xb9\xba\xbb\xbc\xbd\xbe\xbf\xc0\xc1\xc2\xc3\xc4\xc5\xc6\xc7\xc8\xc9\xca\xcb\xcc\xcd\xce\xcf\xd0\xd1\xd2\xd3\xd4\xd5\xd6\xd7\xd8\xd9\xda\xdb\xdc\xdd\xde\xdf\xe0\xe1\xe2\xe3\xe4\xe5\xe6\xe7\xe8\xe9\xea\xeb\xec\xed\xee\xef\xf0\xf1\xf2\xf3\xf4\xf5\xf6\xf7\xf8\xf9\xfa\xfb\xfc\xfd\xfe\xff'

def rc(sequence):
	"""
	Reverse complement assuming DNA nucleotides
	"""

	result = sequence.translate(DNA_RC)
	return result[::-1]


def expand_wildcards(seq, disallow_non_iupac = False):
	"""
	Makes all combination of sequences with wildcards replaced
	by unambiguous bases. Returns a list
	"""
	bases = list(seq.upper())
	seqs = [bases]
	for i, base in enumerate(bases):
		if base not in iupac_wc and disallow_non_iupac:
			return None
	
		expansion = iupac_wc.get(base, ('X',))
		new_seqs = []
		for temp_seq in seqs:
			for x in expansion:
				new_seq = temp_seq[:]
				new_seq[i] = x
				new_seqs.append(new_seq)
		seqs = new_seqs
	seqs = [ ''.join(x) for x in seqs ]
		
	return seqs


def wildcard_stats(seq):
	"""
	Returns number of wildcard positions and number of possible
	expansions
	"""
	seq = seq.upper()
	wc_count = 0
	expansions = 1
	for base in seq.upper():
		if base not in ['A', 'C', 'G', 'T', 'X']:
			wc_count += 1
			expansions *= len(iupac_wc.get(base, ()))
	return seq.count('X'), wc_count, expansions


def gc_count(seq):
	seq = seq.upper()
	count = seq.count('C') + seq.count('G') + seq.count('S')
	count += 0.5 * (seq.count('R') + seq.count('Y') + seq.count('M') + seq.count('K') + seq.count('N'))
	count += 0.67 * (seq.count('B') + seq.count('V'))
	count += 0.33 * (seq.count('D') + seq.count('H'))
	return count/len(seq)
					 

def get_consensus_seq(seqs, wildcard = True):
	"""
	Given aligned sequences, generates a consensus sequence. IUPAC wildcards
	are used. For a gapped position, 'X' is the consensus base.
	"""

	if wildcard:
		return ''.join([iupac_wc_rev.get(tuple(sorted(set(bases))), 'X') for bases in zip(*[list(x) for x in seqs])])
	else:
		return ''.join([canonical_bases.get(tuple(set(bases)), 'X') for bases in zip(*[list(x) for x in seqs])])	

	
def get_consensus_from_fasta(fasta_file, wildcard = True):
	"""
	Given a multiple sequence alignment in FASTA format, returns the consensus
	sequence
	"""
	seqs = mtl_fasta.read_fasta(fasta_file)
	return get_consensus_seq([x[1].upper() for x in seqs], wildcard)


def mask_homopolymers(seq, min_polymer_len = 4):
	"""
	Detects homopolymers of length > min_polymer_len and masks part of them with 'X'
	to disrupt the homopolymer, e.g. AAAAA --> AAXAA
	"""
	
	mask_coords = []
	pad = min_polymer_len - 1
	
	for match in re.finditer(r'((\w)\2{%d,})' % (min_polymer_len - 1), seq):
		start, end = match.span()
		mask_start = start + pad
		mask_end = end - pad
		if mask_start >= mask_end:
			left_adjust = math.floor((mask_start - mask_end)/2)	+ 1
			right_adjust = math.ceil((mask_start - mask_end)/2)
			mask_start -= left_adjust
			mask_end += right_adjust
		mask_coords.append((mask_start, mask_end))

	seq = list(seq)
	for start, end in mask_coords:
		seq[start:end] = ['X'] * (end - start)
	seq = ''.join(seq)
	return seq
	

##################
# Tm calculation
##################

def get_tm_params(molec = 'dnarna'):
	"""
	Return parameters for dH, dS, and helix_init based on molecule type
	"""
	if molec == 'dnarna':
		return RNADNAdeltaH, RNADNAdeltaS, RNA_DNA_HELIX_INIT
	elif molec == 'rna':
		return RNAdeltaH, RNAdeltaS, RNA_HELIX_INIT
	else:
		return DNAdeltaH, DNAdeltaS, DNA_HELIX_INIT


def calculate_tm(dh, ds, helix_init, oligo_concen, na_concen):
	return round( (dh + helix_init)/ (ds + R * math.log((1.0/oligo_concen), math.e)) - 273.15 + 16.6*math.log(na_concen, 10) , 1)


def do_tm(seq, C = 0.00000005, Na_concen = 0.2, molec = 'dnarna'):
	"""
	Formula given by oligocalc
	http://biotools.nubic.northwestern.edu/OligoCalc.html
	(they use a K to C conversion of -272.9,
	so the value will not exactly match)
	(Also, they use the DNA helix initiation energy for RNA)
	
	Na_concen in mol/L
	C in mol/L concen of oligo
	"""
	dH, dS, helix_init = get_tm_params(molec)
		
	sum_dh = 0
	sum_ds = 0

	seq = seq.upper()
	for i in range(len(seq)-1):
		dnt = seq[i:i+2]
		sum_dh += dH[dnt]
		sum_ds += dS[dnt]

	return calculate_tm(sum_dh, sum_ds, helix_init, C, Na_concen)


def do_tm_wc(seq, C = 0.00000005, Na_concen = 0.2, molec = 'dnarna'):
	tms = []
	seqs = expand_wildcards(seq)
	for s in seqs:
		tms.append(do_tm(s, C = C, Na_concen = Na_concen, molec = molec))
	tms.sort()
	return (tms[0], tms[-1], len(seqs))



def iterative_tms_wc(seq, start_at = 0, min_len = 40, max_len = 50, C = 0.00000005, Na_concen = 0.2, molec = 'dnarna', min_tm = 70, max_tm = 80, strict = True, max_wildcard_expansions = 1):

	"""
	Calculates Tm for all prefixes, retaining only the
	ones with acceptable length. Support for wildcard bases,
	which induces a Tm range for each of the possible bases
	Returns a list of Tms: (start_pos, end_pos, Tm_low, Tm_high)
	"""
		
	dH, dS, helix_init = get_tm_params(molec)
	sub_seq = seq[start_at:start_at+max_len]
		
	tms = []
	oligo_lens = []
	
	if len(sub_seq) < min_len:
		return []

	prev_params = [(0, 0, '')]
	for i in range(1, len(sub_seq)):
		dnt = sub_seq[i-1:i+1]
		dnts = expand_wildcards(dnt, disallow_non_iupac = True)
		if not dnts:
			break

		tm_range = []
		curr_params = []

		#Strategy: store separate sets of parameters for each wildcard expansion
		#up to the expansion limit
		for dnt in dnts:
			for sum_dh, sum_ds, base in prev_params:
				if not base or base == dnt[0]:
					sum_dh += dH[dnt]
					sum_ds += dS[dnt]
					curr_params.append((sum_dh, sum_ds, dnt[1]))
					tm_range.append(calculate_tm(sum_dh, sum_ds, helix_init, C, Na_concen))

		prev_params = curr_params
		if len(prev_params) > max_wildcard_expansions:
			break
	
		if i+1 >= min_len:
			tm_low = min(tm_range)
			tm_high = max(tm_range)
			if (tm_low >= min_tm and tm_high <= max_tm) or not strict:
				tms.append((start_at, start_at + i+1, tm_low, tm_high, len(tm_range)))###
				oligo_lens.append(i+1)
	return tms


def generate_tms(seq, min_tm = 70, max_tm = 80, min_len = 40, max_len = 50, C = 0.00000005, Na_concen = 0.2, molec = 'dnarna', strict = True, max_wildcard_expansions = 2):
	"""
	Calculates all possible tms iteratively starting at each position.

	If strict, only oligos with Tm within range are returned

	[ [(start, end, tm_low, tm_high), ...], ]
	
	"""

	tms = []

	for i in range(len(seq)):
		tm_list = iterative_tms_wc(seq, i, min_len, max_len, C, Na_concen, molec, min_tm, max_tm, strict = strict, max_wildcard_expansions = max_wildcard_expansions)
		tms.append(tm_list)

	return tms


##################
# Oligo finding
##################

def wc_trim(seq, trim_amt = 1):
	"""
	Detect non ACGT characters at the beginning or end and trims up to trim_amt off,
	starting with the 5' end
	"""
	
	orig_seq = seq
	max_query_positions = trim_amt * 2
	trim_count = 0
	five_prime = True
	left_offset = 0
	right_offset = 0

	seq = seq.upper()
	while max_query_positions and trim_count < trim_amt and len(seq) > 1:
		if five_prime:
			if (seq[0],) not in canonical_bases:
				seq = seq[1:]
				trim_count += 1
				left_offset += 1
		else:
			if (seq[-1],) not in canonical_bases:
				seq = seq[:-1]
				trim_count += 1
				right_offset -= 1
		five_prime = not five_prime
		max_query_positions -= 1
		
	return seq, left_offset, right_offset


def find_oligo(oligo, seq, reverse_complement = True, expand_wildcard_pos = False, ident = ''):
	"""
	Returns the (start,end) of the oligo in the seq.
	By default, the oligo is assumed to be antisense to a
	target in the seq (reverse_complement = True)

	If the oligo contains wildcards, set expand_wildcard_pos = True to allow all possible oligos
	to be tested
	"""
	
	seq = seq.upper()
	if expand_wildcard_pos:
		oligos = expand_wildcards(oligo)
	else:
		oligos = [oligo]

	locations = []		  
	for olig in oligos:
		olig = olig.upper()
		if reverse_complement:
			olig = rc(olig)

		i = 0
		while True:
			i = seq.find(olig, i)
			if i == -1:
				break
			locations.append((i,i+len(olig), ident, oligo))
			i += 1
	return locations


def best_oligo(tms, start, end, min_tm, max_tm, min_len, max_len, max_untiled_len, return_imperfect = False, right_end_check = False):
	"""
	Find the longest oligo in the window that begins between start and start+max_untiled_len+1.
	If return_imperfect, the oligo with nearest Tm to acceptable is returned if no acceptable one
	is found.
	"""
	
	best_imperfect = None
	wc_backups = {}
	
	for i in reversed(range(min(max_untiled_len+1, end-start-min_len+1))):
		oligo_start = start + i
		
		for oligo in reversed(sorted(tms[oligo_start], key = lambda x: x[1]-x[0])):
			o_start, o_end, tm_low, tm_high, n_wildcard_seqs = oligo
			if o_end > end or (right_end_check and end - o_end > max_untiled_len):
				continue
			
			if tm_low >= min_tm and tm_high <= max_tm:
				if n_wildcard_seqs == 1:  #Favor non wildcard oligos
					return oligo
				else:
					if not wc_backups.get(n_wildcard_seqs, None):
						wc_backups[n_wildcard_seqs] = oligo
			elif return_imperfect:
				score = max(min(abs(tm_low - min_tm), abs(tm_low - max_tm)), min(abs(tm_high - min_tm), abs(tm_high - max_tm)))
				if not best_imperfect or score < best_imperfect[0]:
					best_imperfect = (score, oligo)
	else:
		if return_imperfect and best_imperfect:	 #This depends on whether you care about Tm more or less than wildcards
			return best_imperfect[1]
		elif wc_backups:
			return sorted(wc_backups.items())[0][1]

		return ()


def repair_oligos(oligos, min_len, max_untiled_len, end_index):
	"""
	Hard parameters will sometimes result in not enough sequence
	at the end for a full oligo. In this case, shift oligos left
	until enough space is freed up. Tms will be fixed upon
	oligo refinement
	
	"""
	
	new_oligos_rev = [(end_index - min_len, end_index, -1, -1)]
	
	for i, oligo in enumerate(reversed(oligos)):
		start = oligo[0]
		end = oligo[1]
		if end > new_oligos_rev[-1][0]:
			end = new_oligos_rev[-1][0]
			start = end - min_len
			new_oligos_rev.append((start, end, -1, -1))
		else:
			break
	new_oligos_rev = new_oligos_rev[::-1]
	return oligos[:-1*i] + new_oligos_rev


def greedy_oligo(seq, tms, min_tm = 70, max_tm = 80, min_len = 40, max_len = 50, max_untiled_len = 25):
	"""
	Iterate:
	1. Calculate an optimal window size based on number of remaining nucleotides
	2. Call best_oligo
	"""

	#Optimal window: start with maximum oligo len + 2* max_untiled_len,
	#and shrink down to ensure
	#that window size is a multiple of the sequence length

	remain_len = len(seq)
	oligos = []
	start = 0
	
	while len(seq) - start > max_untiled_len:
		est_n_oligos = math.ceil((remain_len - max_untiled_len) / (max_len + max_untiled_len))
		window_size = math.ceil(remain_len / est_n_oligos)

		if window_size < min_len:
			window_size = min_len

		#Check if the oligo will go over the end - this may happen at the 3' end
		#Terminate oligo finding, and attempt to repair existing solution
		if start + window_size > len(seq):
			oligos = repair_oligos(oligos, min_len, max_untiled_len, len(seq))		
			break

		#Include right-end consistency check for last oligo to make sure there's enough sequence left over at the end
		oligo = best_oligo(tms, start, start + window_size, min_tm, max_tm, min_len, max_len, max_untiled_len, return_imperfect = True, right_end_check = (est_n_oligos == 1))

		if not oligo:
			#Make a placeholder that is a max_len oligo with 1/2 max_untiled_len on each side
			window_size = min(max_len + max_untiled_len, window_size)
			start = start + (window_size - max_len)//2
			end = start + max_len
			oligo = (start, end, -1, -1)
			
		oligos.append(oligo)
		start = oligo[1]
		remain_len = len(seq) - start

	return oligos


def refine_oligos(oligos, seq, tms, min_tm = 70, max_tm = 80, min_len = 40, max_len = 50, max_untiled_len = 25):
	"""
	Given an intermediate solution, attempts to center the oligo in each window
	"""
	
	oligos = oligos + [(len(seq), len(seq), -1, -1, -1)]

	new_oligos = []
	for i in range(len(oligos)-1):
		try:
			start = new_oligos[-1][1]
		except:
			start = 0
		end = oligos[i+1][0]
		new_oligo = best_oligo(tms, start, end, min_tm, max_tm, min_len, max_len, max_untiled_len, return_imperfect = True, right_end_check = True)
		if new_oligo:
			new_oligos.append(new_oligo)
		else:
			new_oligos.append(oligos[i])
		
	return new_oligos



def tile_oligos(seq, min_tm = 70, max_tm = 80, min_len = 40, max_len = 50, C = 0.00000005, Na_concen = 0.2, molec = 'dnarna', max_untiled_len = 25, max_wildcard_expansions = 1):
	"""
	Returns an oligo list tuple: (start, end, Tm)
	"""
	
	if max_untiled_len > len(seq):
		max_untiled_len = len(seq)-1
	
	#Check if a solution exists based on input parameters
	max_n_oligos = len(seq) // min_len
	remaining_nts = len(seq) % min_len
	possible_slop = max_untiled_len + max_n_oligos * (max_len - min_len + max_untiled_len)
	if possible_slop < remaining_nts: #No possible solutions
		return ()

	#SANITIZE SEQUENCE

	#Generate Tms for candidate oligos
	#COULD BE OPTIMIZED USING DP?
	tms = generate_tms(seq, min_tm, max_tm, min_len, max_len, C, Na_concen, molec, strict = False, max_wildcard_expansions = max_wildcard_expansions)


	#Attempt to find an optimal oligo set
	oligos = greedy_oligo(seq, tms, min_tm, max_tm, min_len, max_len, max_untiled_len)

	
	#Refine oligos by centering in the region
	oligos = refine_oligos(oligos, seq, tms, min_tm, max_tm, min_len, max_len, max_untiled_len)
	

	return oligos



def iterative_tile_oligos(seq, min_tm = 70, max_tm = 80, min_len = 40, max_len = 50, C = 0.00000005, Na_concen = 0.2, molec = 'dnarna', max_untiled_len = 25, max_wildcard_expansions = 2, oligo_list = [], reverse_complement = False):
	"""
	- Find previous oligos in full sequence
	- For each intervening region, call greedy_oligo
	- Stitch together the solution (adjust coordinates of the oligos)
	
	"""

	#From the list of provided oligos, find the parts of the sequence
	#that are already targeted

	#Front-pad covered_locs
	covered_locs = [(0,0, '', '')]
	
	for (ident, oligo) in oligo_list:
		if wildcard_stats(oligo)[2] <= max_wildcard_expansions:
			covered_locs += find_oligo(oligo, seq, reverse_complement = reverse_complement, expand_wildcard_pos = True, ident = ident)
	covered_locs.sort()
	
	#End-pad covered_locs
	covered_locs.append((len(seq), len(seq), '', ''))

	oligos = []
	#Iteratively call greedy_oligo for each intervening section
	for i, (start, end, covered_loc_ident, covered_loc_oligo_seq) in enumerate(covered_locs[:-1]):
		if i > 0:
			tm_low, tm_high, n_wildcard_seqs = do_tm_wc(seq[start:end])
			oligos.append((start, end, tm_low, tm_high, covered_loc_ident, covered_loc_oligo_seq))
		r_start = end
		r_end = covered_locs[i+1][0]

		if r_end-r_start > max_untiled_len:
			if r_end-r_start < min_len:
				pass
			#What should be the policy here?
#				 oligos.append((r_start, r_end, -1, -1))
			else:
				subseq = seq[r_start:r_end]
				tms = generate_tms(subseq, min_tm, max_tm, min_len, max_len, C, Na_concen, molec, strict = False, max_wildcard_expansions = max_wildcard_expansions)
				new_oligos = greedy_oligo(subseq, tms, min_tm, max_tm, min_len, max_len, max_untiled_len)
				if not new_oligos:
					oligos.append((r_start, r_end, -1, -1))
				else:
					#Refine oligos by centering in the region
					new_oligos = refine_oligos(new_oligos, subseq, tms, min_tm, max_tm, min_len, max_len, max_untiled_len)
					for new_start, new_end, new_tm_low, new_tm_high, new_wildcard_seq_count in new_oligos:
						oligos.append((new_start+r_start, new_end+r_start, new_tm_low, new_tm_high))
						
	return oligos

	
	
	
def consensus_tiling(idents, seqs, min_len = 39, max_len = 40, max_untiled_len = 15, max_wildcard_expansions = 1, exclude = [], seqs_aligned = True, name_prefix = 'shared', trim_terminal_wildcards = True):
	"""
	Find oligos for the consensus. Designed to work best with a sequence alignment. Otherwise,
	oligos are generated for the first sequence and propagated forward
	"""
	
	if seqs_aligned:
		con_seq = get_consensus_seq([x.upper() for x in seqs], wildcard = (max_wildcard_expansions > 0))
	else:
		con_seq = get_consensus_seq([seqs[0].upper()], wildcard = (max_wildcard_expansions > 0))

	consensus_oligos = []
	oligos = tile_oligos(con_seq, min_len = min_len, max_len = max_len, min_tm=70, max_tm=80, max_untiled_len = max_untiled_len, max_wildcard_expansions = max_wildcard_expansions)

	if trim_terminal_wildcards:
		trim_amt = 4
	else:
		trim_amt = 0
	for oligo in oligos:
		target = con_seq[oligo[0]:oligo[1]]		   
		if 'X' not in target:
			target, left_offset, right_offset = wc_trim(target, trim_amt)
			oligo_name = '%s_%d-%d' % (name_prefix, oligo[0]+1+left_offset, oligo[1]+right_offset)
			consensus_oligos.append((oligo_name, target))

	#Output individual oligo sets for each sequence
	all_oligos = []
		
	for ident, seq in zip(idents, seqs):
		seq = seq.replace('-', '')
		oligos = iterative_tile_oligos(seq, min_len = min_len, max_len = max_len, min_tm=70, max_tm=80, max_untiled_len = max_untiled_len, max_wildcard_expansions = max_wildcard_expansions, oligo_list = consensus_oligos)
		all_oligos.append((ident, seq, oligos))


	return all_oligos
		

################
# OUTPUT
################		

def output_final_oligos(oligo_list, oligo_name_prefix = '', sep = ',', do_header = True):
	"""
	Formats oligo output
	"""
	
	oligo_use_count = {}
	oligo_targets = {}
	oligo_properties = {}
	output_oligos = []
	
	for ident, seq, oligos in oligo_list:	
		for i, oligo in enumerate(oligos):
			target_seq = seq[oligo[0]:oligo[1]]
			if len(oligo) > 5:
				oligo_seq = rc(oligo[5]).upper()  #Shared oligo, possibly with wildcards
			else:
				oligo_seq = rc(target_seq).upper()

			#Cache the oligo with properties
			if oligo_seq not in oligo_use_count:
				tm = oligo[2]
				gc = gc_count(oligo_seq)			
				x_count, wc_count, expand_count = wildcard_stats(oligo_seq)
				oligo_properties[oligo_seq] = (tm, gc, wc_count, expand_count)
				oligo_use_count[oligo_seq] = 1
				oligo_targets[oligo_seq] = {}		
			else:
				oligo_use_count[oligo_seq] += 1
		
			oligo_targets[oligo_seq][target_seq] = 1
		
			#Calculate oligo position in the target sequence
			start, end, oligo_len = oligo[0]+1, oligo[1], oligo[1]-oligo[0]
		
			if i > 0:
				left_gap = oligo[0] - oligos[i-1][1]
			else:
				left_gap = oligo[0]
			try:
				right_gap = oligos[i+1][0] - oligo[1]
			except:
				right_gap = len(seq) - oligo[1]

			#Collate the oligo parameters for output
			output_oligos.append((target_seq, oligo_seq, ident, start, end, oligo_len, left_gap, right_gap))

	#Iterate through the oligos for output, checking use count
	oligo_ids = {}
	output = []
	if do_header:
		output.append(sep.join(['Name', 'Target', 'Start', 'End', 'Length', 'L_gap', 'R_gap', 'Tm', 'GC', 'Wcards', 'Expands', 'N_targets', 'Antisense_oligo', 'Target_seq']))
		
	for target_seq, oligo_seq, ident, start, end, oligo_len, left_gap, right_gap in output_oligos:
		tm, gc, wc_count, expand_count = oligo_properties[oligo_seq]
		num_targets = oligo_use_count[oligo_seq]

		#Consistency check to remove wildcards if they're not needed
		if wc_count > 0 and len(oligo_targets[oligo_seq]) == 1:
			oligo_seq = rc(target_seq)
			x_count, wc_count, expand_count = wildcard_stats(oligo_seq)
			gc = gc_count(oligo_seq)
			#recalculate tm
			tm = do_tm(oligo_seq)

		#Generate a new name or use a previous name if shared oligo
		if oligo_seq in oligo_ids:
			name = oligo_ids[oligo_seq]
		else:
			name = str(len(oligo_ids)+1).zfill(3)
			if oligo_name_prefix:
				name = oligo_name_prefix + '_' + name
			oligo_ids[oligo_seq] = name

		#Assemble output fields
		output.append(sep.join((name, ident, str(start), str(end), str(oligo_len), str(left_gap), str(right_gap), str(tm), '%.2f' % gc, str(wc_count), str(expand_count), str(num_targets), oligo_seq, target_seq)))
	
	return '\n'.join(output)

		
################
# RUN
################		

def do_all(idents, seqs, min_len = 39, max_len = 40, max_untiled_len = 30, do_consensus = False, max_wildcard_expansions = 1, oligo_name_prefix = '', trim_terminal_wildcards = True, seqs_aligned = True, sep =','):
	"""
	Top-level function.
	"""
	
	if do_consensus:
		oligo_list = consensus_tiling(idents, seqs, min_len, max_len, max_untiled_len, max_wildcard_expansions, seqs_aligned = seqs_aligned, name_prefix = 'shared', trim_terminal_wildcards = trim_terminal_wildcards)
		
	else:
		oligo_list = []
		for ident, seq in zip(idents, seqs):
			seq = seq.replace('-', '') #Remove gaps
			oligos = tile_oligos(seq, min_len = min_len, max_len = max_len, min_tm=70, max_tm=80, max_untiled_len = max_untiled_len)
			oligo_list.append((ident, seq, oligos))

	return output_final_oligos(oligo_list, oligo_name_prefix, sep=sep, do_header=True)
			

if __name__ == '__main__':
	my_parser = argparse.ArgumentParser(description='Oligo-ASST: Antisense oligo tiling to design reagents for RNA depletion.')

	my_parser.add_argument('input_fasta',
                       metavar='input_fasta',
                       type=str,
                       help='FASTA file of nucleotide sequence(s) to target.')

	my_parser.add_argument('-o', '--outfile', action = 'store', type= str,
                       help='Filename to write results to (default: <STDOUT>)')
	my_parser.add_argument('-s', '--shared', action = 'store_true',
                       help='Find shared oligos across input sequences')
	my_parser.add_argument('-a', '--aligned', action = 'store_true',
                       help='Input sequences are aligned (aligned FASTA format) (automatically activates -s mode)')
	my_parser.add_argument('-m', '--minlen', action = 'store', type = int, default='39',
                       help='Minimum oligo length (default: 39)')
	my_parser.add_argument('-M', '--maxlen', action = 'store', type = int, default='40',
                       help='Maximum oligo length (default: 40)')
	my_parser.add_argument('-U', '--maxuntiledlen', action = 'store', type = int, default='30',
                       help='Maximum untiled length (default: 30)')
	my_parser.add_argument('-n', '--name', action = 'store', type = str, default='oligo',
                       help='Oligo name prefix (default: oligo)')
	my_parser.add_argument('-W', '--max_wc_expand', action = 'store', type = int, default='1',
                       help='Maximum wildcard expansions (default: 1, i.e. no wildcards allowed. Ignored unless in aligned mode)')
	my_parser.add_argument('-d', '--disable_wc_trim', action = 'store_true',
                       help='Disable trimming terminal wildcards (default: enabled if maxwildcard > 1)')
	
	args = my_parser.parse_args()

	#parse the fasta file into idents, seqs
#	idents, seqs = zip(*mtl_fasta.read_fasta(args.input_fasta))
	idents_long, seqs = mtl_fasta.read_fasta(args.input_fasta, zip_entries = False)

	#Process the identifiers
	ident_hash = {}	
	idents = []
	for ident in idents_long:
		ident = ident.split()[0]
		if ident not in ident_hash:
			ident_hash[ident] = 1
		else:
			ident_hash[ident] += 1
			ident = ident + '(' + str(ident_hash[ident]) + ')'
		idents.append(ident)

	result = do_all(idents, seqs, min_len = args.minlen, max_len = args.maxlen, max_untiled_len = args.maxuntiledlen, do_consensus = args.shared | args.aligned, max_wildcard_expansions = args.max_wc_expand, oligo_name_prefix = args.name, trim_terminal_wildcards = not args.disable_wc_trim, seqs_aligned = args.aligned, sep='\t')
	
	if args.outfile:
		try:
			outf = open(args.outfile, 'w')
			outf.write(result)
			outf.close()
		except:
			print("Error: problem writing to file %s" % args.outfile)
	else:
		print(result)
