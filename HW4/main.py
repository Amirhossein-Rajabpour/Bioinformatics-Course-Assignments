import numpy as np

MAX_SEQUENCE_LENGTH = 10
PSEUDOCOUNT = 2


def find_amino_acids(sequences):
    amino_acids = []
    for seq in sequences:
        for char in seq:
            if char not in amino_acids:
                amino_acids.append(char)
    return amino_acids


def find_num_occurences(amino, index, sequences):
    counter = 0
    for seq in sequences:
        if amino == seq[index]:
            counter += 1
    return counter


def create_PSSM_matrix(sequences):
    num_of_seqs = len(sequences)
    seq_length = len(sequences[0])
    amino_acids = find_amino_acids(sequences)

    profile = np.zeros((len(amino_acids), seq_length))

    for amino_i in range(len(amino_acids)):
        for i in range(seq_length):
            profile[amino_i][i] = find_num_occurences(amino_acids[amino_i], i, sequences)

    profile = profile + PSEUDOCOUNT
    profile = profile / (num_of_seqs + len(amino_acids) * PSEUDOCOUNT)

    for row_i in range(profile.shape[0]):
        profile[row_i] = profile[row_i] / (profile[row_i].sum() / seq_length)

    profile = np.log2(profile)

    return profile, amino_acids


def calculate_score(sequence, profile, amino_acids):
    score = 0
    for char_i in range(len(sequence)):
        amino_index = amino_acids.index(sequence[char_i])
        score += profile[amino_index][char_i]

    return score


def add_one_gap(sequences):
    new_seqs = []
    for seq in sequences:
        for j in range(len(seq)):
            new_seq = seq[:j] + "-" + seq[j:]
            if new_seq not in new_seqs:
                new_seqs.append(new_seq)

        new_seq = seq + "-"
        if new_seq not in new_seqs:
            new_seqs.append(new_seq)

    return new_seqs


def generate_gaps(sequence, num_gaps):
    final_seqs = []
    sequences = [sequence]
    while True:
        generated_seqs = add_one_gap(sequences)
        if len(generated_seqs[0]) == len(sequence) + num_gaps:
            final_seqs += generated_seqs
            break
        sequences = generated_seqs

    return final_seqs


def find_best_subsequence(profile, sequence, amino_acids):
    max_score = float('-inf')
    best_subseq = ""
    seq_length = profile.shape[1]
    for sl in range(seq_length, 0, -1):
        already_checked = []
        for i in range(len(sequence) - sl + 1):
            tmp_seq = sequence[i:i+sl]

            if tmp_seq not in already_checked:
                # here I create different versions with gap
                if len(tmp_seq) < seq_length:
                    tmps_with_gap = generate_gaps(tmp_seq, seq_length - len(tmp_seq))
                    for seq in tmps_with_gap:
                        score_with_gap = calculate_score(seq, profile, amino_acids)
                        if score_with_gap > max_score:
                            best_subseq = seq
                            max_score = score_with_gap

                else:
                    score = calculate_score(tmp_seq, profile, amino_acids)
                    if score > max_score:
                        best_subseq = tmp_seq
                        max_score = score

            already_checked.append(tmp_seq)

    return best_subseq


if __name__ == '__main__':
    number_of_sequences = int(input())
    msa_sequences = []
    for nos in range(number_of_sequences):
        msa_sequences.append(input())

    target_sequence = input()

    PSSM_matrix, amino_acids_array = create_PSSM_matrix(msa_sequences)
    best_subsequence = find_best_subsequence(PSSM_matrix, target_sequence, amino_acids_array)

    print(best_subsequence)
