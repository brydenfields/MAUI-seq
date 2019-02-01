import glob
import sys

global_id_map = {}
global_id_counter = 1

for in_file_name in glob.iglob(sys.argv[1]):
    within_file_count_map = {}

    with open(in_file_name) as in_file:
        for l in in_file:
            seq, _, qual = in_file.__next__().strip(), in_file.__next__(), in_file.__next__().strip()

            if seq not in global_id_map:
                global_id_map[seq] = global_id_counter
                global_seq_id = global_id_counter
                global_id_counter += 1
            else:
                global_seq_id = global_id_map[seq]

            if (global_seq_id, seq) not in within_file_count_map:
                within_file_count_map[(global_seq_id, seq)] = 1
            else:
                within_file_count_map[(global_seq_id, seq)] += 1

    in_file_base_name = in_file_name[in_file_name.rindex("/")+1:]
    sample_id = in_file_base_name[:in_file_base_name.index("_")]
    with open("{}.unique".format(in_file_name), "w") as out_file:
        for (global_seq_id, seq), count in sorted(within_file_count_map.items(), key=lambda kv: kv[1], reverse=True):
            if count >= 100:
                out_file.write(">{}_Sample{}_{}\n{}\n".format(global_seq_id, sample_id, count, seq))
