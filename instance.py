import random
import math

class Oligonucleotide:
    """Sequence spectrum element"""
    def __init__(self, seq):
        self.seq = seq  # oligonucleotide sequence

    def __str__(self):
        return self.seq


class Spectrum:
    """Classic spectrum of SBH - set of oligonucleotide sequences"""
    def __init__(self, sequence, k):
        """Creates classic sequence spectrum without errors"""
        self.sequence = sequence   # original sequence
        self.n = len(sequence)     # len of original sequence
        self.k = k                 # len of oligonucleotide
        self.spectrum = set()      # classic spectrum of nucleotides

        for i in range(self.n - k + 1):
            kmer = sequence[i:i+k]
            oligo = Oligonucleotide(kmer)
            if not self.is_present(kmer):
                self.spectrum.add(oligo)

    def is_present(self, kmer):
        """Checks if oligonucleotide is present in spectrum, and if so, increases its amount"""
        return any(oligo.seq == kmer for oligo in self.spectrum)

    def add_positive_errors(self, percentage, substitutions):
        """Adds % of positive errors into the spectrum (modifies oligonucleotides)"""
        errors_num = math.ceil(percentage * len(self.spectrum) / 100)
        to_alter = random.sample(self.spectrum, errors_num)
        for oligo in to_alter:
            oligo.seq = mutate(oligo.seq, substitutions)

    def add_negative_errors(self, percentage):
        """Adds % of negative errors into the spectrum (removes oligonucleotides)"""
        errors_num = math.ceil(percentage * len(self.spectrum) / 100)
        to_delete = random.sample(self.spectrum, errors_num)
        self.spectrum = [oligo for oligo in self.spectrum if oligo not in to_delete]


class Chromosome:
    """Chromosome is a permutation of oligonucleotides from the spectrum"""
    def __init__(self, spectrum, fitness=0):
        self.seq = spectrum     # oligonucleotides permutation
        self.fitness = fitness

    def __str__(self):
        return self.seq

    def show(self):
        print(self.fitness, end=" ")
        for oligo in self.seq:
            print(oligo, end=" ")
        print()


def mutate(seq, x):
    """Inserts x random point substitutions into the sequence"""
    nucleotides = ['A', 'C', 'G', 'T']
    mutated = list(seq)

    for i in range(x):
        position = random.randint(0, len(seq)-1)
        to_mutate = [n for n in nucleotides if n != mutated[position]]
        n = random.choice(to_mutate)
        mutated[position] = n

    return ''.join(mutated)


def read_sequences(file_name):
    """Reads sequences from a file into a list"""
    sequences = []
    with open(file_name) as file:
        s = []
        file.readline()
        for line in file:
            if line.startswith(">"):
                sequences.append("".join(s))
                s.clear()
            else:
                s.append(line.upper().strip())
        sequences.append("".join(s))

    return sequences