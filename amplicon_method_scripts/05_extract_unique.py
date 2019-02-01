import sys, re

id_to_seq_map = {}

with open(sys.argv[1]) as in_file:
	for header in in_file:
		seq = in_file.__next__()
		header_parts = header[1:].split("_")
		seq_id = header_parts[0]
		sample_id = re.sub('Sample', '', header_parts[1])
		
		if seq_id in id_to_seq_map:
			continue
		
		id_to_seq_map[seq_id] = [seq, sample_id]
		
for seq_id, seqData in id_to_seq_map.items():
	print(">{}\n{}".format(seq_id, seqData[0].strip()))
	#print(seq_id + '\t' + seqData[1])