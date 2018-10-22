import numpy as np

def get_abundance_on_contig(G, con_profile, vertex, align_dict):
    # different algorith to get the abundance
    # get locations for each contig in the window

    height = G.in_degree(vertex) + G.out_degree(vertex) + 1
    width = len(fa_dict[vertex])
    A = np.zeros((height, width))
    B = np.zeros(width)
    C = np.zeros((height, width))
    B[:] += 1
    C[0][:] = 1
    align_locs = []

    map_idx = {}
    map_idx[vertex] = 0
    align_locs.append((1, width, 1, width))
    idx = 1
    pos_start = 1  # blastn starts from 1 to the length of the contig
    pos_end = width
    A[0][pos_start - 1:pos_end] = con_profile[vertex][pos_start - 1:pos_end]
    for succ_node in G.successors(vertex):
        map_idx[succ_node] = idx
        start1, end1, con_len1, start2, end2, con_len2 = G[vertex][succ_node]['align']
        align_seq1, align_seq2 = G[vertex][succ_node]['profile']
        con_align_start1 = pos_start + start1 - 1  # align position on X-Y
        con_align_end1 = end1 - start1 + con_align_start1
        B[start1 - 1:end1] += 1
        C[idx][start1 - 1:end1] = 1
        if (end1 - start1) != (end2 - start2):  # insertion or deletion
            align_contig_profile = get_align_profile(con_profile[succ_node], start2, end2, align_seq2, align_seq1)
            A[idx][con_align_start1 - 1:con_align_end1] = align_contig_profile
            # succ_abund = np.mean(con_profile[succ_node][start2-1:end2])
            # A[idx][con_align_start1-1:con_align_end1] = succ_abund
        else:
            A[idx][con_align_start1 - 1:con_align_end1] += con_profile[succ_node][start2 - 1:end2]
        idx += 1
        align_locs.append((start1, end1, start2, end2))

    for pre_node in G.predecessors(vertex):
        map_idx[pre_node] = idx
        start2, end2, con_len2, start1, end1, con_len1 = G[pre_node][vertex]['align']
        align_seq2, align_seq1 = G[pre_node][vertex]['profile']
        con_align_start1 = pos_start + start1 - 1  # align position on X-Y
        con_align_end1 = end1 - start1 + con_align_start1
        B[start1 - 1:end1] += 1
        C[idx][start1 - 1:end1] = 1
        if (end1 - start1) != (end2 - start2):
            # print pre_node
            align_contig_profile = get_align_profile(con_profile[pre_node], start2, end2, align_seq2, align_seq1)
            A[idx][con_align_start1 - 1:con_align_end1] = align_contig_profile
            # pre_abund = np.mean(con_profile[pre_node][start2-1:end2])
            # A[idx][con_align_start1-1:con_align_end1] = pre_abund
        else:
            A[idx][con_align_start1 - 1:con_align_end1] += con_profile[pre_node][start2 - 1:end2]
        idx += 1
        align_locs.append((start1, end1, start2, end2))

    A_abund = np.zeros((height, width))
    for i in range(height):
        A_abund[i] = np.convolve(A[i], np.ones((100,)) / 100.0, mode='same')
    A_sum = np.sum(A, axis=0)
    A_abund = A_abund / A_sum
    align_contigs = [0] * len(map_idx)
    for con, idx in map_idx.iteritems():
        align_contigs[idx] = con
    print align_contigs

    # get all potential windows
    win_start, win_end = 1, 1
    win_results = []
    for i in range(width - 1):
        if B[i] != B[i + 1]:
            if B[i] > 1:
                win_end = i + 1
                win_tmp = Window(vertex, win_start, win_end, B[i])
                # win_tmp.coverage = A[0][win_start-1:win_end]
                # win_tmp.abundance = A_abund[0][win_start-1:win_end]
                win_tmp.coverage = A
                win_tmp.abundance = A_abund
                win_tmp.aligns = align_contigs
                win_tmp.profile = C
                win_tmp.align_locs = align_locs
                # win_tmp.abund_mean = win_tmp.get_abund_mean()
                win_results.append(win_tmp)

            if B[i + 1] > 1:
                win_start = i + 2
    if win_start > win_end and B[width - 1] > 1:
        win_end = width
        win_tmp = Window(vertex, win_start, win_end, B[width - 1])
        # win_tmp.coverage = A[0][win_start-1:win_end]
        # win_tmp.abundance = A_abund[0][win_start-1:win_end]
        win_tmp.coverage = A
        win_tmp.abundance = A_abund
        win_tmp.aligns = align_contigs
        win_tmp.profile = C
        win_tmp.align_locs = align_locs
        # win_tmp.abund_mean = win_tmp.get_abund_mean()
        win_results.append(win_tmp)
    return win_results, A_abund, map_idx, align_locs

