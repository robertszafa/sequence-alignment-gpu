import sys
import random

DELETION_CHANCE = 0.05
INSERTION_CHANCE = 0.02
SUBSTITUTION_CHANCE = 0.05

if __name__ == "__main__":
    fname = str(sys.argv[-1])
    proteinOrDna = str(sys.argv[1])

    # Default dna.
    if len(sys.argv) == 1:
        proteinOrDna = 'dna'

    m_fname = fname.split('/')
    m_fname[-1] = "mutated_" + m_fname[-1]
    m_fname = ''.join([text + '/' for text in m_fname])[:-1]


    DNA_ALPHABET =  ['A', 'T', 'C', 'G']
    PROTEIN_ALPHABET =  ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P',
                         'S', 'T', 'W', 'Y', 'V', 'B', 'Z', 'X']

    alphabet = DNA_ALPHABET if 'd' in proteinOrDna or len(sys.argv) == 1 else PROTEIN_ALPHABET

    mutated_seq = '>> mutation of ' + fname + ' by mutate.py\n'
    mutated_seq += f'>> DELETION_CHANCE = {DELETION_CHANCE}\n'
    mutated_seq += f'>> INSERTION_CHANCE = {INSERTION_CHANCE}\n'
    mutated_seq += f'>> SUBSTITUTION_CHANCE = {SUBSTITUTION_CHANCE}\n\n'

    num_del = 0
    num_inserts = 0
    num_sub = 0
    with open(fname, 'r') as f:
        for line in f:
            # FASTA encoded files.
            if line.lstrip()[0] == '>':
                mutated_seq += line
                continue

            for c in line:
                c = c.upper()

                # 5% chance of deletion
                if random.random() < DELETION_CHANCE:
                    num_del += 1
                    pass
                # 3% chance of insertion
                elif random.random() < INSERTION_CHANCE:
                    num_inserts += 1
                    mutated_seq += random.choice(alphabet)
                # 5% chance of substitution
                elif random.random() < SUBSTITUTION_CHANCE:
                    num_sub += 1
                    alphabet_without_c = [l for l in alphabet if l != c]
                    mutated_seq += random.choice(alphabet_without_c)
                else:
                    mutated_seq += c

    print("Number of deletion:\t", num_del)
    print("Number of insertions:\t", num_inserts)
    print("Number of substitutions:", num_sub)

    print("\nSaving to", m_fname)

    with open(m_fname, 'w') as f:
        f.write(mutated_seq)
