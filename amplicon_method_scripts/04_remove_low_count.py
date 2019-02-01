# -*- coding: utf-8 -*-
# Kør med "python 04_remove_low_count.py > nodA_min75.fasta"
import sys
from itertools import groupby

# Den her funktion kan læse en fasta fil, så man får header og sekvens.
def fasta_iter(fasta_file):
    """
    Given a FASTA file, yield tuples of header, sequence.
    Author: Brent Pedersen (https://www.biostars.org/p/710/#1412)
    """
    fh = open(fasta_file)
    # Ditch the boolean (x[0]) and just keep the header or sequence since we know they alternate.
    faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
    for header in faiter:
        # Drop the ">".
        header = header.next()[1:].strip()
        # Join all sequence lines to one in case it spans multiple lines.
        seq = "".join(s.strip() for s in faiter.next())
        yield header, seq


# Åben filen (erstat med dit filnavn) og læs hver header-sekvens par.
for header, seq in fasta_iter("nodA.fasta."):
    # Split headeren på "_" og læs den 3. værdi (tæl fra 0) som et heltal.
    count = int(header.split("_")[2])
    # Hvis count er mindre end 500, så gå til næste sekvens.
    if count < 75:
        continue

    # Ellers skriver vi headeren og sekvensen ud på skærmen.
    sys.stdout.write(">{}\n{}\n".format(header, seq))
