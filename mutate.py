import sys
import random

if __name__ == "__main__":
    proteinOrDna = str(sys.argv[1])

    fname = str(sys.argv[2])
    m_fname = fname.split('/')
    m_fname[-1] = "mutated_" + m_fname[-1]
    m_fname = ''.join([text + '/' for text in m_fname])[:-1]


    DNA_ALPHABET =  ['A', 'T', 'C', 'G']
    PROTEIN_ALPHABET =  ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P',
                         'S', 'T', 'W', 'Y', 'V', 'B', 'Z', 'X']

    alphabet = DNA_ALPHABET if 'd' in proteinOrDna else PROTEIN_ALPHABET

    mutated_seq = ''

    num_del = 0
    num_inserts = 0
    num_sub = 0
    with open(fname, 'r') as f:
        for line in f:
            for c in line:
                c = c.upper()

                # 5% chance of deletion
                if random.random() < 0.05:
                    num_del += 1
                    pass
                # 3% chance of insertion
                elif random.random() < 0.03:
                    num_inserts += 1
                    mutated_seq += random.choice(alphabet)
                # 5% chance of substitution
                elif random.random() < 0.05:
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
