import pdb
import numpy as np
from lib import *
from ReadFiles import *


def get_align_profile(profile, start, end, align_seq1, align_seq2):
    assert len(align_seq1)==len(align_seq2), "Align sequences error!"
    align_len = len(align_seq1)
    A = np.zeros(align_len)
    con_profile = profile[start-1:end]
    idx = 0
    rm_idx = []
    for i in range(align_len):
        if align_seq1[i]!='-':
            A[i] = con_profile[idx]
            idx+=1
        if align_seq2[i]=='-':
            rm_idx.append(i)
    A = np.delete(A, rm_idx)
    return A

def get_abundance_on_contig(G, con_profile, vertex, fa_dict, align_dict):
    # get locations for each contig in the window
    height = G.in_degree(vertex) + G.out_degree(vertex) + 1
    width = len(fa_dict[vertex])
    A = np.zeros((height, width))
    B = np.zeros(width) # Array to record the aligned contigs depth
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
        con_align_start1 = pos_start + start1 - 1  # align position on X-Y, same as on reference contig
        con_align_end1 = end1 - start1 + con_align_start1
        B[start1 - 1:end1] += 1
        C[idx][start1 - 1:end1] = 1
        if abs(end1 - start1) != abs(end2 - start2):  # insertion or deletion
            align_contig_profile = get_align_profile(con_profile[succ_node], start2, end2, align_seq2, align_seq1)
            A[idx][con_align_start1 - 1:con_align_end1] = align_contig_profile
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
        win_tmp.coverage = A
        win_tmp.abundance = A_abund
        win_tmp.aligns = align_contigs
        win_tmp.profile = C
        win_tmp.align_locs = align_locs
        win_results.append(win_tmp)
    return win_results, A_abund, map_idx, align_locs

"""
EM clustering
"""

def get_distribution(cov_array):
    #if len(cov_array) == 0: # empty group
    #    return np.array([])
    length = max(cov_array)
    hist = np.zeros(int(length)+1)
    for cov in cov_array:
        hist[int(cov)]+=1
    #hist = np.convolve(hist, np.ones((100,))/100.0, mode='same')
    hist = hist/np.sum(hist)
    return hist

def cal_prob(cov_array, distribute):
    #if len(distribute) == 0 or len(cov_array) == 0: # empty group
    #    return 0
    prob_array = []
    for cov in cov_array:
        if cov < len(distribute):
            prob = distribute[int(cov)]
            prob_array.append(prob)
        else:
            #prob_array.append(0)
            prob_array.append(distribute[-1])
    return sum(prob_array)/float(len(prob_array))

def assign_subcontig(con, distribution, bin_num):
    #cov_array = con.coverage
    cov_array = bin_num*con.abundance
    max_pro = 0
    max_idx = 0
    idx = 0
    for distribute in distribution:
        prob = cal_prob(cov_array, distribute)
        if max_pro<prob:
            max_pro = prob
            max_idx = idx
        idx+=1
    return max_idx

def assign_subcontig_by_posterior(con, cluster):
    con_probs = [x/sum(con.probs) for x in con.probs]

    priors = [0]*len(cluster)
    for group in cluster:
        for con_tmp in group:
            if con_tmp.contig == con.contig and (con_tmp.start!=con.start and con_tmp.end!=con.end):
                con_tmp_probs = [x/sum(con_tmp.probs) for x in con_tmp.probs]
                for i in range(len(con_tmp.probs)):
                    #pdb.set_trace()
                    priors[i] += con_tmp_probs[i]
    if sum(priors)==0:
        #if con.contig=='20485|126':
        #    pdb.set_trace()
        priors = [1]*len(cluster)
    else:
        priors = [x/sum(priors) for x in priors]
    return priors

def get_cluster(contig_win, k):
    cluster = []
    for i in range(k):
        cluster.append([])

    for contig in contig_win:
        cons = contig_win[contig]
        for con in cons:
            cluster[con.label].append(con)
    return cluster

def adjust_label(contig_win, k):
    cluster = []
    for i in range(k):
        cluster.append([])

    for contig in contig_win:
        cons = contig_win[contig]
        con_label = {}
        for con in cons:
            if not con.label in con_label:
                con_label[con.label] = con.length
            else:
                con_label[con.label] += con.length
        max_weight = 0
        max_label = 0
        #pdb.set_trace()
        for label in con_label:
            if max_weight<con_label[label]:
                max_label = label
                max_weight = con_label[label]
        for con in contig_win[contig]:
            con.label = max_label
            cluster[con.label].append(con)
    return cluster

def adjust_label_by_prob(contig_win, k, distribution, bin_num):
    cluster = []
    for i in range(k):
        cluster.append([])

    for contig in contig_win:
        cons = contig_win[contig]
        con_label = {}
        for con in cons:
            con_array = bin_num*con.abundance
            probs = [cal_prob(con_array, distribute) for distribute in distribution]
            if not con.label in con_label:
                con_label[con.label] = con.length*probs[con.label]
            else:
                con_label[con.label] += con.length*probs[con.label]
        max_weight = 0
        max_label = 0
        #pdb.set_trace()
        for label in con_label:
            if max_weight<con_label[label]:
                max_label = label
                max_weight = con_label[label]
        for con in contig_win[contig]:
            con.label = max_label
            cluster[con.label].append(con)
    return cluster

def sample_label(posters):
    return np.random.choice(np.arange(5), p=posters)

def format_cluster(cluster):
    cluster_new = []
    for group in cluster:
        group_new = set([])
        for con in group:
            con_format = con.contig+'_'+str(con.start)+'_'+str(con.end)
            group_new.add(con_format)
        cluster_new.append(group_new)
    return cluster_new

def is_same_cluster(cluster1, cluster2):
    if len(cluster1)!=len(cluster2):
        return 0
    flag = 1
    cluster1_new = format_cluster(cluster1)
    cluster2_new = format_cluster(cluster2)
    for i in range(len(cluster1)):
        if cluster1_new[i]!=cluster2_new[i]:
            flag = 0
    return flag

def EM_initialization(win_list):
    pass

def EM_cluster_gibbs(win_list, k, bin_num):
    # preprocessing
    tmp_win_list = []
    contig_win = {}
    for win in win_list:
        if len(win) < k - 1 or len(win) > k + 1 or win[0].length < 10:
            # if len(win)!=k or win[0].length<10:
            continue
        tmp_win_list.append(win)
    print len(tmp_win_list)

    # initialization
    for win in tmp_win_list:
        win = sorted(win, key=lambda x: np.mean(x.abundance))
        if len(win) <= k:
            for i in range(len(win)):
                win[i].label = i
        else:
            for i in range(k):
                win[i].label = i
            for i in range(k, len(win)):
                win[i].label = k - 1

        for con in win:
            if not con.contig in contig_win:
                contig_win[con.contig] = [con]
            else:
                contig_win[con.contig].append(con)

    # adjust label for subcontig on the same contig
    cluster = [[]] * k
    cluster = get_cluster(contig_win, k)

    new_cluster = cluster[:]
    cluster = []
    iteration = 0
    distribution = []
    adjust_cluster = []

    while (iteration < 100 and (not is_same_cluster(new_cluster, cluster))):
        cluster = new_cluster[:]

        # Expectation
        distribution = []
        for group in cluster:
            cov_array = np.array([])
            for con in group:
                cov_array = np.append(cov_array, bin_num * con.abundance)
            if len(cov_array) == 0:
                pdb.set_trace()
            distribute_tmp = get_distribution(cov_array)  # function
            distribution.append(distribute_tmp)

        for group in cluster:  # update likelihood
            for con in group:
                probs = np.array([])
                for distribute in distribution:
                    prob = cal_prob(bin_num * con.abundance, distribute)
                    probs = np.append(probs, prob)
                con.probs = probs

        # Maximization
        contig_win.clear()
        for group in cluster:
            for con in group:
                priors = assign_subcontig_by_posterior(con, cluster)  # function
                con.priors = priors
                # con.priors = np.array([1,1,1,1,1])
                posters = [0] * k
                for i in range(len(priors)):
                    posters[i] = priors[i] * con.probs[i]
                posters = [x / sum(posters) for x in posters]
                con.label = posters.index(max(posters))
                if not con.contig in contig_win:
                    contig_win[con.contig] = [con]
                else:
                    contig_win[con.contig].append(con)
        new_cluster = get_cluster(contig_win, k)
        adjust_cluster = adjust_label_by_prob(contig_win, k, distribution, bin_num)
        iteration += 1
    print 'Iteration:', iteration
    return new_cluster, distribution, adjust_cluster


def evaluate_bp(cluster, align_dict):
    TP_length = dict(zip(ground_truth, [0]*len(ground_truth)))
    ground_truth_length = {}
    pred_length = {}
    idx = 0
    for group in cluster:
        for con in group:
            ref = align_dict[con.contig]
            ref = ref.split('_')[0]
            pred_ref = ground_truth[idx]
            ground_truth_length[ref] = ground_truth_length.get(ref, 0) + con.length
            pred_length[pred_ref] = pred_length.get(pred_ref, 0) + con.length

            if ref==pred_ref:
                TP_length[ref] = TP_length.get(ref, 0) + con.length
        idx += 1

    precision, recall = {}, {}
    for ref in ground_truth:
        if not ref in pred_length or (not ref in ground_truth_length):
            precision[ref] = 0
            recall[ref] = 0
        else:
            precision[ref] = float(TP_length[ref])/pred_length[ref]
            recall[ref] = float(TP_length[ref])/ground_truth_length[ref]
    return precision, recall

def evaluate_bp_whole(cluster, align_dict):
    #TP_length = {'HXB2':0, 'YU2':0, '89.6':0, 'NL43':0, 'JRCSF':0}
    TP_length = dict(zip(ground_truth, [0]*len(ground_truth)))
    ground_truth_length = {}
    pred_length = {}
    idx = 0
    contig_len_whole = 0
    tp_len_whole = 0
    for group in cluster:
        for con in group:
            contig_len_whole += con.length
            ref = align_dict[con.contig]
            ref = ref.split('_')[0]
            pred_ref = ground_truth[idx]
            ground_truth_length[ref] = ground_truth_length.get(ref, 0) + con.length
            pred_length[pred_ref] = pred_length.get(pred_ref, 0) + con.length

            if ref==pred_ref:
                TP_length[ref] = TP_length.get(ref, 0) + con.length
                tp_len_whole += con.length
        idx += 1

    acc = float(tp_len_whole)/contig_len_whole

    return acc

