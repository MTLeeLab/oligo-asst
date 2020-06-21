"""
This file is part of Oligo-ASST, released under the 3-Clause BSD License.
See License.txt or go to https://opensource.org/licenses/BSD-3-Clause for
full details.

Copyright 2020 Miler T. Lee.


mtl_fasta.py: Basic FASTA file reading
Author: Miler Lee
"""


def read_fasta(fastafile, zip_entries = True):
	"""
	Wrapper around parse_fasta to work on a file.
	"""
	
	f = open(fastafile)
	return parse_fasta(f.read(), zip_entries)


def parse_fasta(fastalines, zip_entries = True):
	"""
	Interprets the input as a string encoding an
	entire FASTA file.

	If zip_entries then the return is a list of (id, seq) tuples.
	Otherwise return is two lists
	"""

	ids = []
	seqs = []

	curr_id = ''
	curr_seq = []

	lines = fastalines.splitlines()
	for line in lines:
		if line.startswith('>'):
			#Write current entry
			if curr_seq:
				if not curr_id:
					curr_id = 'UNKNOWN'
				ids.append(curr_id)
				seqs.append(''.join(curr_seq))
				curr_id = ''
				curr_seq = []
				
			curr_id = line[1:].strip()
		else:
			curr_line = line.strip()
			if curr_line:
				curr_seq.append(curr_line)
	
	
	#Finally
	if curr_seq:
		if not curr_id:
			curr_id = 'UNKNOWN'
		ids.append(curr_id)
		seqs.append(''.join(curr_seq))

	if zip_entries:
		return list(zip(ids, seqs))
	else:
		return ids, seqs
		

