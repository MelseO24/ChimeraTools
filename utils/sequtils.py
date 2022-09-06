import sys
from Bio import SeqIO

class sequence_info:
    def __init__(self, input_fasta):
        """Saves information from fasta file
         Provide input: fasta file, containing one or more sequences"""
        self.seqcounter = 0
        self.input_fasta = input_fasta
        self.fasta_seq = SeqIO.parse(open(input_fasta), 'fasta')

        self.seqdata = {}
        for fasta in self.fasta_seq:
            self.seqcounter += 1
            if fasta.id in self.seqdata.keys():
                sys.stderr.write(f"FATAL: multiple sequences with id: '{fasta.id}' found!.\n")
                sys.exit(1)
            self.seqdata[fasta.id] = fasta.seq

        print(f"Analyzed {self.seqcounter} sequences.")

    def get_seq(self, seqid):
        return self.seqdata[seqid]

    def get_length(self, seqid):
        return len(self.seqdata[seqid])

    def get_all_ids(self):
        return self.seqdata.keys()

    def __iter__(self):
        yield