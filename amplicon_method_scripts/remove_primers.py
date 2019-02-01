import re
import sys

forward_primers = [".*[ACT]G[CT]TCGCAGTGGTGGATGTT", ".*[ACT]CGAGAATGTTGTCGAGAT[CT]GAGACGA", ".*[ACT]ATGCGTTTTAAGGG[AC][CT]TGGATCT", ".*[ACT]CCGGATCT[GC]GAGGGGCT"]
reverse_primers = ["ATGCCGTTCCT[TC]GAAGACGG", "CAGGCGCT[TC]GAAATCACCGATA", "GCGGAATCTGACCGG[TC]GC", "TCATTGATCGAACGAAACGG[ATCG]TGC"]

num_changes = 0
with open(sys.argv[1]) as in_file:
    for name in in_file:
        name = name.strip()
        seq, _, qual = in_file.__next__().strip(), in_file.__next__(), in_file.__next__().strip()

        len_before = len(seq)
        new_seq = seq
        for p in forward_primers:
            new_seq = re.sub(p, "", new_seq)

        len_after = len(new_seq)
        num_removed = len(seq) - len(new_seq)
        if num_removed > 0:
            num_changes += 1
        else:
            continue

        new_qual = qual[num_removed:]

        for p in reverse_primers:
            m = re.search(p, new_seq)
            if m is not None:
                new_seq = new_seq[:m.start()]
                # new_qual = new_qual[:m.start()]
                break

        print("{}\n{}\n+\n{}".format(name, new_seq, new_qual))
        # print(">{}\n{}".format(name[1:], new_seq))

# print(num_changes)
