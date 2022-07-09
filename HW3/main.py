import copy

S_MATCH = 3
S_MISMATCH = -1
S_GAP = -2


def global_align(x, y, s_match, s_mismatch, s_gap):
    A = []
    for i in range(len(y) + 1):
        A.append([0] * (len(x) + 1))

    for i in range(len(y) + 1):
        A[i][0] = s_gap * i

    for i in range(len(x) + 1):
        A[0][i] = s_gap * i

    for i in range(1, len(y) + 1):
        for j in range(1, len(x) + 1):
            A[i][j] = max(
                A[i][j - 1] + s_gap,
                A[i - 1][j] + s_gap,
                A[i - 1][j - 1] + (s_match if (y[i - 1] == x[j - 1] and y[i - 1] != '-') else 0) + (
                    s_mismatch if (y[i - 1] != x[j - 1] and y[i - 1] != '-' and x[j - 1] != '-') else 0) + (
                    s_gap if (y[i - 1] == '-' or x[j - 1] == '-') else 0)
            )

    align_X = ""
    align_Y = ""

    i = len(x)
    j = len(y)

    while i > 0 or j > 0:
        current_score = A[j][i]
        if i > 0 and j > 0 and (
                ((x[i - 1] == y[j - 1] and y[j - 1] != '-') and current_score == A[j - 1][i - 1] + s_match) or
                ((y[j - 1] != x[i - 1] and y[j - 1] != '-' and x[i - 1] != '-') and current_score == A[j - 1][
                    i - 1] + s_mismatch) or
                ((y[j - 1] == '-' or x[i - 1] == '-') and current_score == A[j - 1][i - 1] + s_gap)
        ):

            align_X = x[i - 1] + align_X
            align_Y = y[j - 1] + align_Y

            i = i - 1
            j = j - 1

        elif i > 0 and (current_score == A[j][i - 1] + s_gap):
            align_X = x[i - 1] + align_X
            align_Y = "-" + align_Y
            i = i - 1

        else:
            align_X = "-" + align_X
            align_Y = y[j - 1] + align_Y
            j = j - 1

    return align_X, align_Y, A[len(y)][len(x)]


def calculate_pairwise_scores(sequences):
    num_of_sequences = len(sequences)
    # calculate pairwise scores
    sum_scores = []
    sequences_tmp = sequences + sequences  # just to make the loop easier
    for i in range(num_of_sequences):
        seq_scores = []
        for j in range(i + 1, i + num_of_sequences):
            seq_scores.append(global_align(sequences_tmp[i], sequences_tmp[j], S_MATCH, S_MISMATCH, S_GAP)[2])
        sum_scores.append(sum(seq_scores))
    return sum_scores


def find_center(sum_scores, sequences):
    # find center sequence
    center_sequence_index = sum_scores.index(max(sum_scores))
    center_sequence = sequences[center_sequence_index]
    return center_sequence, center_sequence_index


def extract_seqs(aligned_sequences):
    seqs = []
    for i in aligned_sequences.values():
        seqs.append(i["sequence"])
    return seqs


def calculate_score(sequence1, sequence2):
    score = 0
    length = len(sequence1)
    for i in range(length):
        if sequence1[i] == sequence2[i] == '-':  # gap, gap
            pass
        elif (sequence1[i] == '-') or (sequence2[i] == '-'):  # gap
            score -= 2
        elif sequence1[i] == sequence2[i]:  # match
            score += 3
        elif sequence1[i] != sequence2[i]:  # mismatch
            score -= 1
    return score


def calculate_scores(aligned_sequences):
    score = 0
    num_of_sequences = len(aligned_sequences)
    for i in range(num_of_sequences):
        for j in range(i + 1, num_of_sequences):
            score += calculate_score(aligned_sequences[i], aligned_sequences[j])
    return score


def insert_new_gaps(sequence, array_of_gaps):
    modified_sequence = ""
    for i in range(len(sequence)):
        if i in array_of_gaps:
            repeated_gap = array_of_gaps.count(i)
            if repeated_gap > 1:
                for rp in range(repeated_gap-1):
                    modified_sequence += "-"

            modified_sequence += "-"
            modified_sequence = modified_sequence + sequence[i]
        else:
            modified_sequence = modified_sequence + sequence[i]

    for gap in array_of_gaps:
        if gap >= len(sequence):
            modified_sequence += "-"

    return modified_sequence


def find_gaps(old_center, new_center):
    new_gap_indexes = []
    len_o = len(old_center)
    len_n = len(new_center)
    o = 0
    for n in range(len_n):
        if n >= len_o:
            break
        else:
            if old_center[o] != '-' and new_center[n] == '-':
                new_gap_indexes.append(o)
            else:
                o += 1

    gap_at_end = (len_n - len_o) - len(new_gap_indexes)
    for i in range(gap_at_end):
        index = len_o + i
        new_gap_indexes.append(index)

    return new_gap_indexes


def star_alignment(sequences, center_index):
    # doing pairwise alignment between center and other sequences and save aligned centers
    center_scores_with_other_seqs = {}
    for i in range(0, len(sequences)):
        if i != center_index:
            seq_i, aligned_center, score_with_center = global_align(sequences[i], sequences[center_index], S_MATCH, S_MISMATCH, S_GAP)
            center_scores_with_other_seqs[i] = score_with_center
    center_scores_with_other_seqs = {k: v for k, v in
                                     sorted(center_scores_with_other_seqs.items(), key=lambda item: item[1],
                                            reverse=True)}

    aligned_sequences = {}
    for i in range(len(sequences)):
        aligned_sequences[i] = {}
        aligned_sequences[i]["sequence"] = sequences[i]

    # sequences should be aligned in order of their distance (or score) from the center
    # find indexes that are not center
    other_seq_indexes = list(center_scores_with_other_seqs.keys())
    once_aligned_seqs = []
    for i in other_seq_indexes:

        # align with new sequence
        old_center = copy.deepcopy(aligned_sequences[center_index]["sequence"])
        aligned_sequences[i]["sequence"], aligned_sequences[center_index]["sequence"], _ = global_align(aligned_sequences[i]["sequence"], aligned_sequences[center_index]["sequence"], S_MATCH, S_MISMATCH, S_GAP)
        gaps = find_gaps(old_center, aligned_sequences[center_index]["sequence"])

        # update sequences that are already aligned
        for j in once_aligned_seqs:
            # just insert gaps
            aligned_sequences[j]["sequence"] = insert_new_gaps(aligned_sequences[j]["sequence"], gaps)

        once_aligned_seqs.append(i)
        once_aligned_seqs = set(once_aligned_seqs)
        once_aligned_seqs = list(once_aligned_seqs)

    array_final_seqs = extract_seqs(aligned_sequences)
    new_score = calculate_scores(array_final_seqs)

    return new_score, array_final_seqs


def modify_alignment(sequences):

    # score initialization
    new_final_score = calculate_scores(sequences)

    # find index of columns that can't be in a block
    col_not_in_block = []
    for i in range(len(sequences[0])):
        same_char_col = True
        check_char = sequences[0][i]
        for j in range(len(sequences)):
            if sequences[j][i] != check_char:
                same_char_col = False
        if same_char_col:
            col_not_in_block.append(i)
    # print("columns that can't be in a block: ", col_not_in_block)

    # find length of blocks. we need blocks with length of greater than 1
    block_lengths = []
    col_not_in_block.insert(0, -1)
    col_not_in_block.append(len(sequences[0]))
    for i in range(1, len(col_not_in_block)):
        block_lengths.append(col_not_in_block[i] - col_not_in_block[i - 1] - 1)
    # print("block lengths: ", block_lengths)

    # find range of each block. their start and end indexes
    blocks_to_be_modified = []  # (start, end) of each block which size is greater than 1. it includes the start and the end index themselves
    for i in range(len(block_lengths)):
        if block_lengths[i] > 1:
            block_range = (col_not_in_block[i] + 1, col_not_in_block[i + 1] - 1)
            blocks_to_be_modified.append(block_range)
    # print("blocks ranges: ", blocks_to_be_modified)

    # extract each block from their original sequences and then apply MSA and check the new score
    for r in blocks_to_be_modified:
        # print("----------new block----------")
        start = r[0]
        end = r[1]
        seqs_in_blocks = []
        for s in sequences:
            seqs_in_blocks.append(s[start:end + 1])
        # print("sequences in block: ", seqs_in_blocks)
        seqs_in_block_with_gaps = copy.deepcopy(seqs_in_blocks)  # just to check later for detecting changes in block

        # remove their gaps
        for i in range(len(seqs_in_blocks)):
            seqs_in_blocks[i] = seqs_in_blocks[i].replace("-", "")

        # find their center
        array_sum_scores = calculate_pairwise_scores(seqs_in_blocks)
        block_center, block_center_index = find_center(array_sum_scores, seqs_in_blocks)

        # apply MSA
        block_score, modified_block_seqs = star_alignment(seqs_in_blocks, block_center_index)
        # print("modified block: ", modified_block_seqs)

        # check if there is any change in the block. if block is unchanged, no need for further calculation
        seqs_in_block_changed = False
        for i in range(len(seqs_in_blocks)):
            if seqs_in_block_with_gaps[i] != modified_block_seqs[i]:
                seqs_in_block_changed = True

        # compare new MSA score of the modified sequences with previous MSA score
        if seqs_in_block_changed:
            prev_score = calculate_scores(seqs_in_block_with_gaps)
            new_score = calculate_scores(modified_block_seqs)
            # print("prev score: ", prev_score, " *****  modified score: ", new_score)

            if new_score > prev_score:
                # replace the modified block
                for s in range(len(sequences)):
                    sequences[s] = sequences[s][:start] + modified_block_seqs[s] + sequences[s][end+1:]

                new_final_score = calculate_scores(sequences)

    return new_final_score, sequences


if __name__ == '__main__':

    num_of_sequences = int(input())
    input_sequences = []
    for i in range(num_of_sequences):
        input_sequences.append(input())

    sum_scores = calculate_pairwise_scores(input_sequences)
    center_sequence, center_sequence_index = find_center(sum_scores, input_sequences)

    star_score, final_sequences = star_alignment(input_sequences, center_sequence_index)

    modified_score, modified_sequences = modify_alignment(final_sequences)
    print(modified_score)
    for fs in modified_sequences:
        print(fs)
