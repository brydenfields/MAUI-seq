#!/usr/bin/env python

#This script reads sequences from a set of fastq files listed in file_list.txt.
#The first seqid_len bases are a random tag (seqid); the script keeps track of the number of unique seqids associated with each distinct sequence in each sample.
#The sequence associated with each seqid is the most frequent one, provided it has at least 2 reads with that seqid.
#Output file 1 is a fasta file with sequences in decreasing order of overall abundance, down to a definable theshold.
#Output file 2 is a tab-delimited table of the number of counts of each of these sequences in each of the sample files.
#Written by Peter Young. 

#The following parameters are specific for each data set and gene
working_folder = "/Users/peter/Documents/Sequence_data/MiSeq_2018-05-31_NCHAIN/assembled/nodD/"
seqid_len = 13
f_primer_len = 23
r_primer_len = 18
total_len = 308


all_samples_table = {}
total_counts = {}

with open(working_folder + "file_list.txt") as file_list:
    for fastq_filename in file_list:

        if "fastq" in fastq_filename:
            fastq_filepath = working_folder + fastq_filename.rstrip("\n")

            sample_ID = fastq_filename[0:fastq_filename.find("_")]
            seqid_dict = {}
            total_seqs = 0

            #Function find_match to split a sequence line into seqid and sequence (removing the primers);
            #   then make a dictionary of all the different seqids with counts of each of their sequences.
            def find_match(line,dic):
                seqid = line[0:seqid_len]
                sequence = line[(seqid_len + f_primer_len):(len(line) - r_primer_len)]
                if seqid in dic:
                    find_sequence(seqid,sequence,dic)
                else:
                    dic[seqid] = {sequence:1}
    
            def find_sequence(id,seq,dic):
                if seq in dic[id]:
                    dic[id][seq] += 1
                else:
                    dic[id][seq] = 1
    
            #Read in a fastq file one sequence record at a time (4 lines) and process the DNA sequence (2nd line) with find.match
            with open(fastq_filepath) as fastq_file:
                ctr = 0
                record = []
                for next_line in fastq_file:
                    record.append(next_line.rstrip("\n"))
                    ctr += 1
                    if ctr == 4:
                        if total_len -2 <= len(record[1]) <= total_len + 2: #Only process sequences that are expected length +/- 2 bases 
                            find_match(record[1],seqid_dict)
                            total_seqs += 1
                        record = []
                        ctr = 0
        
            sample_counts = {}  #keys will be sequences, items will be number of seqids for each sequence in the sample
            for seqid, matches in sorted(seqid_dict.items(), key=lambda item: sum(item[1].values()), reverse=True):
           
                for sequence in sorted(seqid_dict[seqid].items(), key=lambda item: item[1], reverse=True)[0:1]: #Choose the most abundant sequence for each seqid.
                    if sequence[1] >=2: #Only include sequences that have at least two reads with that seqid.                   
                        if sequence[0] in sample_counts:
                            sample_counts[sequence[0]] += 1
                        else:
                            sample_counts[sequence[0]] = 1
                        if sequence[0] in total_counts:
                            total_counts[sequence[0]] += 1
                        else:
                            total_counts[sequence[0]] = 1
        
        
            all_samples_table[sample_ID] = sample_counts

fasfilename = working_folder + "ranked_total_sequences_03c.fas"
fasfile = open(fasfilename, "w")
rank = 0
ranked_sequence_list = []

for sequence, seq_count in sorted(total_counts.items(), key=lambda item: item[1], reverse=True):
    if seq_count > 32:   #omit sequences with an average of less than one seqid per sample
        rank += 1
        fasfile.write(">seq_%d_%d" % (rank, seq_count))
        fasfile.write("\n")
        fasfile.write(sequence)
        fasfile.write("\n")
  
        ranked_sequence_list.append(sequence)
    
fasfile.close()   
        
outfilename = working_folder + "table_by_seqids_03c.tab"
outfile = open(outfilename, "w")
for sampleID in sorted(all_samples_table):
    outfile.write(sampleID+"\t")
    for sequence in ranked_sequence_list:
        if sequence in all_samples_table[sampleID]:
            seq_count = all_samples_table[sampleID][sequence]
        else:
            seq_count = 0
        outfile.write(str(seq_count)+"\t")
    outfile.write("\n")
    
outfile.close()
        
 