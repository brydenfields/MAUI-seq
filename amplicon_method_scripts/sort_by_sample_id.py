import sys

with open(sys.argv[1]) as in_file:
	triples = []
	for header in in_file:
		seq = in_file.__next__()

		sample_id = int(header.strip().split('_')[1][6:])

		triples.append((sample_id, header, seq))

	# Sort the triples according to sample ID.
	sorted_triples = sorted(triples, key=lambda t: t[0])

	print(''.join('{}{}'.format(header, seq) for _, header, seq in sorted_triples))
