import random
import math
import numpy as np
import instance as inst

# POPULATION
def create_initial_population(spectrum, size):
    """Creates an initial population of chromosomes of the given size"""
    initial_population = [inst.Chromosome(spectrum.copy()) for _ in range(size)]
    return initial_population

def greedy_create_initial_population(spectrum, overlap_len):
    """Creates an initial population of chromosomes of the given size. 
    Chromosomes are created in a semi-deterministic way"""
    initial_population = []

    for i in range(len(spectrum)):
        chromosome = []
        used = [False for _ in range(len(spectrum))]
        chromosome.append(spectrum[i])
        used[i] = True
        start = i
        while not all(used):
            next_oligo = greedy_best_overlap(spectrum[start], spectrum, used, overlap_len)
            chromosome.append(spectrum[next_oligo])
            used[next_oligo] = True
            start = next_oligo

        ch = inst.Chromosome(chromosome)
        initial_population.append(ch)

    return initial_population

def greedy_best_overlap(first, spectrum, used, overlap_len):
    """Helper function for greedy_create_initial_population"""
    max_overlap = 0
    best = 0
    for i in range(len(spectrum)):
        if not used[i]:
            overlap = check_overlap(first.seq, spectrum[i].seq, overlap_len)
            if overlap >= max_overlap:
                max_overlap = overlap
                best = i
    return best


# EVALUATION
def check_overlap(first, second, overlap_len):
    """Returns the length of the prefix overlap on the suffix 
    with the minimum length of overlap_len or 0 if no overlap"""
    max_overlap = min(len(first), overlap_len)
    for overlap_index in range(max_overlap, 0, -1):
        if second.startswith(first[-overlap_index:]):
            return overlap_index
    return 0

def fitness(chromosome, n, overlap_len):
    """Evaluates a chromosome, returning normalized fitness value"""
    sequence = [chromosome.seq[0].seq]
    seq_len = len(sequence[0])
    k = len(sequence[0])

    for i in range(len(chromosome.seq) - 1):
        overlap = check_overlap(chromosome.seq[i].seq, chromosome.seq[i + 1].seq, overlap_len)
        to_append = chromosome.seq[i + 1].seq[overlap:]
        sequence.append(to_append)
        seq_len += len(to_append)

    value = len(sequence) / (n - k + 1) * 1.5
    fitness_value = round(value, 4)
    chromosome.fitness = fitness_value
    return fitness_value

def evaluate(population, n, overlap_len):
    """Assigns a fitness value to all chromosomes in the population, 
    returning result as a list of tuples: (chromosome, fitness)"""
    evaluated = [(chromosome, fitness(chromosome, n, overlap_len)) for chromosome in population]
    return evaluated

def best_chromosome(evaluated):
    """Returns the chromosome with the best fitness value 
    as a tuple: (chromosome, fitness)"""
    return max(evaluated, key=lambda tup: tup[1])


# SELECTION
def selection(population, evaluation, pop_size):
    """Selects elements of the population using the roulette method"""
    fitness_values = np.array([x[1] for x in evaluation])
    probabilities = fitness_values / np.sum(fitness_values)

    selected_indices = np.random.choice(len(population), pop_size, p=probabilities)
    selected_population = [population[i] for i in selected_indices]
    np.random.shuffle(selected_population)
    return selected_population


# CROSSOVER (partially deterministic)
def is_present(oligo, chromosome):
    """Checking the presence of an oligonucleotide in a chromosome"""
    return any(o and o.seq == oligo.seq for o in chromosome)

def best_overlap(first, spectrum, child, overlap_len):
    """Returns the best overlap"""
    max_overlap = 0
    best = None
    for oligo in spectrum:
        if not is_present(oligo, child):
            overlap = check_overlap(first.seq, oligo.seq, overlap_len)
            if overlap >= max_overlap:
                max_overlap = overlap
                best = oligo
    return best

def greedy_crossover(parent1, parent2, spectrum, overlap_len):
    """Greedy crossover operator using overlap evaluation"""
    spectrum_size = len(spectrum)
    l = len(parent1.seq)
    locus = random.randint(0, l-1)
    oligo = random.choice(spectrum)

    child1 = [None] * l
    child1[locus] = oligo
    prev = locus

    while not all(child1):
        p1_oligo = parent1.seq[locus]
        p2_oligo = parent2.seq[locus]

        if (is_present(p1_oligo, child1)) or (is_present(p2_oligo, child1)):
            best_oligo = best_overlap(child1[prev], spectrum, child1, overlap_len)
            child1[locus] = best_oligo
        else:
            p1_overlap = check_overlap(child1[prev].seq, p1_oligo.seq, overlap_len)
            p2_overlap = check_overlap(child1[prev].seq, p2_oligo.seq, overlap_len)
            child1[locus] = p1_oligo if p1_overlap >= p2_overlap else p2_oligo

        prev = locus
        locus = (locus + 1) % l  # wrap around if necessary

    ch1 = inst.Chromosome(child1)
    return ch1

def enhanced_crossover(population, spectrum, percent, overlap_len):
    """Crossover of chromosomes from a population by greedy crossover method"""
    random.shuffle(population)
    indices = list(range(len(population)))
    random.shuffle(indices)

    uneven = len(indices) % 2
    s = len(indices) - 1 if uneven else len(indices)
    pairs = [(population[indices[i]], population[indices[i+1]]) for i in range(0, s, 2)]

    decision = [0, 1]
    parents = []
    children = []
    for pair in pairs:
        p = percent / 100
        probability = [1-p, p]
        if np.random.choice(decision, 1, replace=False, p=probability):
            child1 = greedy_crossover(pair[0], pair[1], spectrum, overlap_len)
            children.append(child1)
        else:
            parents.extend(list(pair))

    if uneven:
        parents.append(population[indices[-1]])
    return parents, children


# MUTATION
def swap_mutation(chromosome):
    """Mutation of a random chromosome involving random switching of two oligonucleotides"""
    i1, i2 = sorted(random.sample(range(len(chromosome.seq)), 2))
    chromosome.seq[i1], chromosome.seq[i2] = chromosome.seq[i2], chromosome.seq[i1]

def mutate(children, mut_percent):
    """Performs mutations on population elements with given probabilities"""
    num_of_mutating_chromosomes = math.floor(len(children) * mut_percent / 100)
    for _ in range(num_of_mutating_chromosomes):
        index = random.randrange(0, len(children))
        swap_mutation(children[index])


# RESULT SEQUENCE
def get_sequence(chromosome, n):
    """Returns result sequence of the chromosome"""
    sequence = []
    seq_len = 0
    overlap_len = 5

    k = len(chromosome.seq[0].seq)
    sequence.append(chromosome.seq[0].seq)
    seq_len += k

    for i in range(0, len(chromosome.seq)):
        overlap = check_overlap(chromosome.seq[i].seq, chromosome.seq[i+1].seq, overlap_len)
        to_append = chromosome.seq[i + 1].seq[overlap:]
        sequence.append(to_append)
        seq_len += len(to_append)

        if seq_len == n:
            return "".join(sequence)
        elif seq_len > n:
            sequence.pop()
            return "".join(sequence)

