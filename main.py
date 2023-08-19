import instance as inst
import operators as op
import comparison as comp

if __name__ == '__main__':

    BEST_GEN = 1         # number of generation with best fitness
    BEST_CH = (None, 0)  # best chromosome and its fitness

    # READ SEQUENCE
    file_name = "sequences.txt"
    SEQUENCES = inst.read_sequences(file_name)
    sequence = SEQUENCES[0]
    n = len(sequence)

    # SET PARAMETERS
    MAX_GEN = 60
    init_pop_size = 10000
    population_size = 60     # population size has to be even
    cross_percent = 90       # crossover probability
    mutation_percent = 2     # mutation probability
    k = 10                   # k-mer length
    overlap_len = 7          # length of overlap for hybridization
    substitutions = 4        # max = k
    negative = 5             # % of negative errors
    positive = 5             # % of positive errors

    # GENETIC ALGORITHM
    S = inst.Spectrum(sequence.upper(), k)
    S.add_positive_errors(positive, substitutions)
    S.add_negative_errors(negative)

    generation = 1
    init_population = op.create_initial_population(S.spectrum, init_pop_size)
    init_evaluation = op.evaluate(init_population, n, overlap_len)

    best_ch = op.best_chromosome(init_evaluation)
    if best_ch[1] > BEST_CH[1]:
        BEST_CH = best_ch

    selected = op.selection(init_population, init_evaluation, population_size)
    parents, children = op.enhanced_crossover(selected, S.spectrum, cross_percent, overlap_len)
    op.mutate(children, mutation_percent)
    population = parents + children
    evaluation = op.evaluate(population, n, overlap_len)

    best_ch = op.best_chromosome(evaluation)
    if best_ch[1] > BEST_CH[1]:
        BEST_CH = best_ch
        BEST_GEN = generation
    generation += 1

    while generation <= MAX_GEN:
        print(f"Generation {generation}")
        selected = op.selection(population, evaluation, population_size)
        parents, children = op.enhanced_crossover(selected, S.spectrum, cross_percent, overlap_len)
        op.mutate(children, mutation_percent)
        population = parents + children
        evaluation = op.evaluate(population, n, overlap_len)

        best_ch = op.best_chromosome(evaluation)
        if best_ch[1] > BEST_CH[1]:
            BEST_CH = best_ch
            BEST_GEN = generation
        generation += 1

    result = op.get_sequence(BEST_CH[0], n)

    print("\noriginal:", sequence)
    print("result:  ", result)

    print("\nPairwise Alingment:")
    al1, al2 = comp.pairwise_alignment(sequence, result)
    print(al1)
    print(al2)

    distance = comp.levenshtein_distance(sequence, result)
    similarity = comp.similarity(distance, sequence, result)
    print(f"best_gen={BEST_GEN} n={n} k={k} levenshtein distance={distance} similarity={similarity}% fitness={BEST_CH[1]}\n")
