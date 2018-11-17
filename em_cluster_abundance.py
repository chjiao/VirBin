import re,sys,pdb
import matplotlib
import pygraphviz
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import networkx as nx
import numpy as np
from lib import *
from sklearn.cluster import KMeans
import math

# plot the alignment results for blastn
# update with BFS search 

def get_abundance(G, con_profile):
    # G: a connected graph
    height = len(G)
    start_nodes = [N for N in G if G.in_degree(N)==0]
    vertex = start_nodes[0]
    pos_start, pos_end = 0,0
    while G.out_degree(vertex)>0:
        succ_node = G.successors(vertex)[0]
        start1, end1, con_len1, start2, end2, con_len2 = G[vertex][succ_node]['align']
        con_align_start1 = pos_start + start1
        con_align_end1 = end1-start1+con_align_start1
        con_start1 = con_align_start1 + (1-start1)
        con_end1 = con_start1+con_len1 - 1
        if pos_end<con_end1:
            pos_end = con_end1
        con_start2 = con_align_start1 + (1-start2)
        con_end2 = con_start2 + con_len2 -1
        if pos_end<con_end2:
            pos_end = con_end2
        vertex = succ_node
    
    width = pos_end
    A = np.zeros((height, width))
    vertex = start_nodes[0]
    index = 0
    while G.out_degree(vertex)>0:
        succ_node = G.successors(vertex)[0]
        start1, end1, con_len1, start2, end2, con_len2 = G[vertex][succ_node]['align']
        con_align_start1 = pos_start + start1
        con_align_end1 = end1-start1+con_align_start1
        con_start1 = con_align_start1 + (1-start1)
        con_end1 = con_start1+con_len1 - 1
        profile1 = con_profile[vertex]
        #pdb.set_trace()
        A[height-index-1][con_start1-1:con_end1] = profile1
        con_start2 = con_align_start1 + (1-start2)
        con_end2 = con_start2 + con_len2 -1
        profile2 = con_profile[succ_node]
        vertex = succ_node
        index+=1
        if G.out_degree(succ_node)==0: # last node
            A[height-index-1][con_start2-1:con_end2] = profile2
    for i in range(height):
        A[i] = np.convolve(A[i], np.ones((100,))/100.0, mode='same')
    A_sum = np.sum(A, axis=0)
    A = A/A_sum
    return A

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

ref_dict={'HXB2':0.1087, 'JRCSF':0.2979, 'YU2':0.1262, '89.6':0.2168, 'NL43':0.2518}
def get_ref(ref_dict, mean):
    con_ref = ''
    min_distance=10
    for ref in ref_dict:
        if abs(mean-ref_dict[ref])<min_distance:
            min_distance = abs(mean-ref_dict[ref])
            con_ref = ref
    return con_ref

def get_abundance_contig(G, con_profile, vertex, align_dict):
    height = G.in_degree(vertex)+G.out_degree(vertex)+1
    width = len(fa_dict[vertex])
    A = np.zeros((height, width))
    B = np.zeros(width)
    B[:]+=1

    map_idx = {}
    map_idx[vertex] = 0
    idx = 1
    pos_start = 1  # blastn starts from 1 to the length of the contig
    pos_end = width
    A[0][pos_start-1:pos_end] = con_profile[vertex][pos_start-1:pos_end]
    for succ_node in G.successors(vertex):
        map_idx[succ_node] = idx
        start1, end1, con_len1, start2, end2, con_len2 = G[vertex][succ_node]['align']
        align_seq1, align_seq2 = G[vertex][succ_node]['profile']
        con_align_start1 = pos_start + start1 -1 # align position on X-Y
        con_align_end1 = end1-start1+con_align_start1
        B[start1-1:end1]+=1
        if (end1-start1)!=(end2-start2): # insertion or deletion
            print succ_node
            align_contig_profile = get_align_profile(con_profile[succ_node], start2, end2, align_seq2, align_seq1)
            A[idx][con_align_start1-1:con_align_end1] += align_contig_profile
        else:
            A[idx][con_align_start1-1:con_align_end1] += con_profile[succ_node][start2-1:end2]
        idx+=1
        
    for pre_node in G.predecessors(vertex):
        map_idx[pre_node] = idx
        start2, end2, con_len2, start1, end1, con_len1 = G[pre_node][vertex]['align']
        align_seq2, align_seq1 = G[pre_node][vertex]['profile']
        con_align_start1 = pos_start+ start1 -1 # align position on X-Y
        con_align_end1 = end1-start1+con_align_start1
        B[start1-1:end1] +=1
        if (end1-start1)!=(end2-start2):
            print pre_node
            align_contig_profile = get_align_profile(con_profile[pre_node], start2, end2, align_seq2, align_seq1)
            A[idx][con_align_start1-1:con_align_end1] += align_contig_profile
        else:
            A[idx][con_align_start1-1:con_align_end1] += con_profile[pre_node][start2-1:end2]
        idx +=1
    
    A_abund = np.zeros((height, width))
    for i in range(height):
        A_abund[i] = np.convolve(A[i], np.ones((100,))/100.0, mode='same')
    A_sum = np.sum(A, axis=0)
    A_abund = A_abund/A_sum

    # get all potential windows
    win_start, win_end = 1,1
    win_results = []
    for i in range(width-1):
        if B[i]!=B[i+1]:
            if B[i]>1:
                win_end = i+1
                win_tmp = Window(vertex, win_start, win_end, B[i])
                win_tmp.coverage = A[0][win_start-1:win_end]
                win_tmp.abundance = A_abund[0][win_start-1:win_end]
                win_tmp.abund_mean = win_tmp.get_abund_mean()
                win_results.append(win_tmp)

            if B[i+1]>1:
                win_start = i+2
    if win_start>win_end and B[width-1]>1:
        win_end = width
        win_tmp = Window(vertex, win_start, win_end, B[width-1])
        win_tmp.coverage = A[0][win_start-1:win_end]
        win_tmp.abundance = A_abund[0][win_start-1:win_end]
        win_tmp.abund_mean = win_tmp.get_abund_mean()
        win_results.append(win_tmp)
    return win_results,A_abund, map_idx

def get_abundance_contig2(G, con_profile, vertex, align_dict):
    # different algorith to get the abundance

    height = G.in_degree(vertex)+G.out_degree(vertex)+1
    width = len(fa_dict[vertex])
    A = np.zeros((height, width))
    B = np.zeros(width)
    C = np.zeros((height, width))
    B[:]+=1
    C[0][:] = 1
    

    map_idx = {}
    map_idx[vertex] = 0
    idx = 1
    pos_start = 1  # blastn starts from 1 to the length of the contig
    pos_end = width
    A[0][pos_start-1:pos_end] = con_profile[vertex][pos_start-1:pos_end]
    for succ_node in G.successors(vertex):
        map_idx[succ_node] = idx
        start1, end1, con_len1, start2, end2, con_len2 = G[vertex][succ_node]['align']
        align_seq1, align_seq2 = G[vertex][succ_node]['profile']
        con_align_start1 = pos_start + start1 -1 # align position on X-Y
        con_align_end1 = end1-start1+con_align_start1
        B[start1-1:end1]+=1
        C[idx][start1-1:end1] = 1
        if (end1-start1)!=(end2-start2): # insertion or deletion
            #align_contig_profile = get_align_profile(con_profile[succ_node], start2, end2, align_seq2, align_seq1)
            succ_abund = np.mean(con_profile[succ_node][start2-1:end2])
            A[idx][con_align_start1-1:con_align_end1] = succ_abund 
        else:
            A[idx][con_align_start1-1:con_align_end1] += con_profile[succ_node][start2-1:end2]
        idx+=1
        
    for pre_node in G.predecessors(vertex):
        map_idx[pre_node] = idx
        start2, end2, con_len2, start1, end1, con_len1 = G[pre_node][vertex]['align']
        align_seq2, align_seq1 = G[pre_node][vertex]['profile']
        con_align_start1 = pos_start+ start1 -1 # align position on X-Y
        con_align_end1 = end1-start1+con_align_start1
        B[start1-1:end1] +=1
        C[idx][start1-1:end1] = 1
        if (end1-start1)!=(end2-start2):
            #print pre_node
            #align_contig_profile = get_align_profile(con_profile[pre_node], start2, end2, align_seq2, align_seq1)
            pre_abund = np.mean(con_profile[pre_node][start2-1:end2])
            A[idx][con_align_start1-1:con_align_end1] = pre_abund
        else:
            A[idx][con_align_start1-1:con_align_end1] += con_profile[pre_node][start2-1:end2]
        idx +=1
    
    A_abund = np.zeros((height, width))
    for i in range(height):
        A_abund[i] = np.convolve(A[i], np.ones((100,))/100.0, mode='same')
    A_sum = np.sum(A, axis=0)
    A_abund = A_abund/A_sum
    align_contigs = [0]*len(map_idx)
    for con, idx in map_idx.iteritems():
        align_contigs[idx] = con
    print align_contigs

    # get all potential windows
    win_start, win_end = 1,1
    win_results = []
    for i in range(width-1):
        if B[i]!=B[i+1]:
            if B[i]>1:
                win_end = i+1
                win_tmp = Window(vertex, win_start, win_end, B[i])
                #win_tmp.coverage = A[0][win_start-1:win_end]
                #win_tmp.abundance = A_abund[0][win_start-1:win_end]
                win_tmp.coverage = A
                win_tmp.abundance = A_abund
                win_tmp.aligns = align_contigs
                win_tmp.profile = C
                #win_tmp.abund_mean = win_tmp.get_abund_mean()
                win_results.append(win_tmp)

            if B[i+1]>1:
                win_start = i+2
    if win_start>win_end and B[width-1]>1:
        win_end = width
        win_tmp = Window(vertex, win_start, win_end, B[width-1])
        #win_tmp.coverage = A[0][win_start-1:win_end]
        #win_tmp.abundance = A_abund[0][win_start-1:win_end]
        win_tmp.coverage = A
        win_tmp.abundance = A_abund
        win_tmp.aligns = align_contigs
        win_tmp.profile = C
        #win_tmp.abund_mean = win_tmp.get_abund_mean()
        win_results.append(win_tmp)
    return win_results,A_abund, map_idx


def cluster_windows(win_list):
    results = []
    wins1 = []
    wins2 = win_list[:]
    while len(wins1)!=len(wins2) and len(wins2)>0:
        wins1 = wins2[:]
        win_init = wins1[0]
        for win in wins1[1:]:
            if abs(win_init.abund_mean-win.abund_mean)<0.025:
                #pdb.set_trace()
                win_init.coverage = np.concatenate((win_init.coverage, win.coverage))
                win_init.abundance = np.concatenate((win_init.abundance, win.abundance))
                win_init.abund_mean = win_init.get_abund_mean()
                win_init.contig.extend(win.contig)
                wins2.remove(win)
        win_init.contig = list(set(win_init.contig))
        results.append(win_init)
        wins2.pop(0)
    return results

def cluster_windows_k(win_list, k, align_dict, ref_dict, f_out):
    X,y=[],[]
    for win in win_list:
        X.append(win.abund_mean)
        y.append(win.contig[0])
    X=np.asarray(X)
    X = X.reshape(-1,1)
    y_pred = KMeans(n_clusters=k).fit_predict(X)
    #pdb.set_trace()
    for j in range(5):
        f_out.write('>Cluster '+str(j)+'\n')
        for i in range(len(y_pred)):
            y_i = y_pred[i]
            if y_i==j:
                con = y[i]
                ref = align_dict[con]
                ref = ref.split('_')[0]
                f_out.write(str(y[i])+'\t'+str(X[i][0])+'\t'+ref+'\t'+str(ref_dict[ref])+'\n')
    f_out.close()

def merge_window(win1, win2):
    for con1 in win1:
        for con2 in win2:
            if con1.contig==con2.contig:
                con1.abundance = np.append(con1.abundance, con2.abundance)
                con1.coverage = np.append(con1.coverage, con2.coverage)
    return win1

def merge_same_windows(win_list):
    tmp_win_list = win_list[:]
    result = []
    while(tmp_win_list):
        win = tmp_win_list[0]
        win_cons = []
        for con in win:
            win_cons.append(con.contig)
        remove_idx = [0]
        for i in range(1, len(tmp_win_list)):
            win2 = tmp_win_list[i]
            win2_cons = []
            for con in win2:
                win2_cons.append(con.contig)
            if set(win_cons)==set(win2_cons):
                win = merge_window(win, win2)
                remove_idx.append(i)
        result.append(win)
        for idx in sorted(remove_idx, reverse=True):
            tmp_win_list.pop(idx)
    return result

def same_window(win1, win2):
    flag = 0
    for con in win2:
        if con.contig==win1[0].contig:
            if abs(con.start - win1[0].start)<50 and abs(con.end - win1[0].end)<50:
                flag = 1
    return flag

def merge_duplicate_windows(win_list):
    tmp_win_list = win_list[:]
    result = []
    while(tmp_win_list):
        win = tmp_win_list[0]
        if len(win)==0:
            pdb.set_trace()
        tmp_results = [win]
        win_cons = []
        for con in win:
            win_cons.append(con.contig)
        remove_idx = [0]
        for i in range(1, len(tmp_win_list)):
            win2 = tmp_win_list[i]
            win2_cons = []
            for con in win2:
                win2_cons.append(con.contig)
            if set(win_cons)==set(win2_cons) and same_window(win, win2):
                tmp_results.append(win2)
                #win = merge_window(win, win2)
                remove_idx.append(i)
        result.append(tmp_results)
        for idx in sorted(remove_idx, reverse=True):
            tmp_win_list.pop(idx)
    return result

def distance(win1, win2):
    dist = 0
    tmp_win1 = win1[:]
    tmp_win2 = win2[:]
    for con1 in win1:
        for con2 in win2:
            if con1.contig==con2.contig:
                tmp_win1.remove(con1)
                tmp_win2.remove(con2)
                break
    tmp_win1 = sorted(tmp_win1, key = lambda con:np.mean(con.abundance))
    tmp_win2 = sorted(tmp_win2, key = lambda con:np.mean(con.abundance))
    flag = 1
    if len(tmp_win1)!=len(tmp_win2):
        pdb.set_trace()
    for i in range(len(tmp_win1)):
        dist+= (tmp_win1[i].get_abund_mean() - tmp_win2[i].get_abund_mean())**2
        if tmp_win1[i].contig!= tmp_win2[i].contig:
            flag = 0
    return math.sqrt(dist), flag

def get_win_distance(win1, win2):
    dist = 0
    for i in range(len(win1)):
        con1 = win1[i]
        con2 = win2[i]
        con1_array = np.array([])
        con2_array = np.array([])
        for con in con1:
            con1_array = np.append(con1_array, con.abundance)
        for con in con2:
            con2_array = np.append(con2_array, con.abundance)

        dist+= (np.mean(con1_array) - np.mean(con2_array))**2
    return math.sqrt(dist)

    
def cluster_nearest_window(win_list, f):
    tmp_win_list = win_list[:]
    idx = 0
    res, total = 0,0
    while(len(tmp_win_list)>1):
        win = tmp_win_list[0]
        flag = 0
        for i in range(1, len(tmp_win_list)):
            win2 = tmp_win_list[i]
            if i==1:
                min_distance, flag = distance(win, win2)
                nearest_win = win2
            else:
                if min_distance>distance(win, win2):
                    min_distance, flag = distance(win, win2)
                    nearest_win = win2
                    idx = i
        tmp_win_list.pop(idx)
        tmp_win_list.pop(0)
        
        win = sorted(win, key=lambda con:np.mean(con.abundance), reverse=True)
        win2 = sorted(win2, key=lambda con:np.mean(con.abundance), reverse=True)
        f.write('>Pair+\t'+str(min_distance)+'\n')
        for j in range(len(win)):
            con1 = win[j]
            con2 = win2[j]
            f.write(con1.contig+'\t'+str(con1.get_abund_mean())+'\t'+con1.ref+'\t'+con2.contig+'\t'+str(con2.get_abund_mean())+'\t'+con2.ref+'\n')
        if flag: # correctly merged
            res+=1
        total +=1
    return res, total

def merge_two_window(win, win2):
    result = []
    for i in range(len(win)):
        con1 = win[i]
        con2 = win2[i]
        con1.extend(con2)
        result.append(con1)
    return result

def cluster_nearest_window2(win_list):
    #tmp_win_list = win_list[:]
    tmp_win_list = []
    for win in win_list:
        win = sorted(win, key=lambda con:np.mean(con.abundance), reverse=True)
        res = []
        for con in win:
            res.append([con])
        tmp_win_list.append(res)

    idx1, idx2 = 0, 0
    res, total = 0,0
    min_distance = 0
    while(len(tmp_win_list)>1):
        for i in range(0, len(tmp_win_list)-1):
            for j in range(i+1, len(tmp_win_list)):
                win = tmp_win_list[i]
                win2 = tmp_win_list[j]

                if i==0 and j==1:
                    min_distance = get_win_distance(win, win2)
                    idx1 = i
                    idx2 = j
                else:
                    if min_distance>get_win_distance(win, win2):
                        min_distance= get_win_distance(win, win2)
                        idx1 = i
                        idx2 = j
        
        win = merge_two_window(tmp_win_list[idx1], tmp_win_list[idx2])

        tmp_win_list.append(win)

        tmp_win_list.pop(idx1)
        tmp_win_list.pop(idx2)

        print len(tmp_win_list), min_distance

    return tmp_win_list

def get_distribution(cov_array):
    len = max(cov_array)
    hist = np.zeros(int(len)+1)
    for cov in cov_array:
        hist[int(cov)]+=1
    #hist = np.convolve(hist, np.ones((100,))/100.0, mode='same')
    hist = hist/np.sum(hist)
    return hist

def cal_prob(cov_array, distribute):
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

def EM_cluster(win_list, k, bin_num):
    #bin_num = 100
    # preprocessing
    tmp_win_list = []
    contig_win = {}
    for win in win_list:
        #if len(win)<k-1 or len(win)>k+1 or win[0].length<10:
        if len(win)!=k or win[0].length<10:
            continue
        tmp_win_list.append(win)
    print len(tmp_win_list)

    # initialization
    for win in tmp_win_list:
        # need change
        #win = sorted(win, key = lambda x:np.mean(x.coverage))
        win = sorted(win, key = lambda x:np.mean(x.abundance))
        if len(win)<=k:
            for i in range(len(win)):
                win[i].label = i
        else:
            for i in range(k):
                win[i].label = i
            for i in range(k, len(win)):
                win[i].label = k-1
        
        for con in win:
            if not con.contig in contig_win:
                contig_win[con.contig] = [con]
            else:
                contig_win[con.contig].append(con)

    # adjust label for subcontig on the same contig
    cluster = [[]]*k
    #cluster = adjust_label(contig_win, k)
    cluster = get_cluster(contig_win, k)
    #pdb.set_trace()
    
    new_cluster = cluster[:]
    cluster = []
    iteration = 0
    distribution = []
    adjust_cluster = []
    while(iteration<200 and new_cluster!=cluster):
        cluster = new_cluster[:]
        #pdb.set_trace()
        # Expectation
        distribution = []
        for group in cluster:
            cov_array = np.array([])
            for con in group:
                #cov_array = np.append(cov_array, con.coverage)
                cov_array = np.append(cov_array, bin_num*con.abundance)
            distribute_tmp = get_distribution(cov_array) # function
            distribution.append(distribute_tmp)
        
        #pdb.set_trace()
        # Maximization
        contig_win.clear()
        for group in cluster:
            for con in group:
                max_idx = assign_subcontig(con, distribution) # function
                con.label = max_idx
                if not con.contig in contig_win:
                    contig_win[con.contig] = [con]
                else:
                    contig_win[con.contig].append(con)
                #new_cluster[max_idx].append(con)
        #new_cluster = adjust_label(contig_win, k) # function
        new_cluster = get_cluster(contig_win, k)
        #adjust_cluster = adjust_label(contig_win, k)
        adjust_cluster = adjust_label_by_prob(contig_win, k, distribution, bin_num)
        iteration += 1
    print 'Iteration:', iteration
    return new_cluster, distribution, adjust_cluster

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

def EM_cluster_gibbs(win_list, k, bin_num):
    # preprocessing
    tmp_win_list = []
    contig_win = {}
    for win in win_list:
        #if len(win)<k-1 or len(win)>k+1 or win[0].length<10:
        if len(win)!=k or win[0].length<10:
            continue
        tmp_win_list.append(win)
    print len(tmp_win_list)

    # initialization
    for win in tmp_win_list:
        # need change
        #win = sorted(win, key = lambda x:np.mean(x.coverage))
        win = sorted(win, key = lambda x:np.mean(x.abundance))
        if len(win)<=k:
            for i in range(len(win)):
                win[i].label = i
        else:
            for i in range(k):
                win[i].label = i
            for i in range(k, len(win)):
                win[i].label = k-1
        
        for con in win:
            if not con.contig in contig_win:
                contig_win[con.contig] = [con]
            else:
                contig_win[con.contig].append(con)

    # adjust label for subcontig on the same contig
    cluster = [[]]*k
    cluster = get_cluster(contig_win, k)
    #pdb.set_trace()
    
    new_cluster = cluster[:]
    cluster = []
    iteration = 0
    distribution = []
    adjust_cluster = []
    #while(iteration<100 and new_cluster!=cluster):
    while(iteration<100 and (not is_same_cluster(new_cluster, cluster))):
        cluster = new_cluster[:]
        #pdb.set_trace()
        # Expectation
        distribution = []
        for group in cluster:
            cov_array = np.array([])
            for con in group:
                #cov_array = np.append(cov_array, con.coverage)
                cov_array = np.append(cov_array, bin_num*con.abundance)
            if len(cov_array)==0:
                pdb.set_trace()
            distribute_tmp = get_distribution(cov_array) # function
            distribution.append(distribute_tmp)

        for group in cluster: # update likelihood
            for con in group:
                probs = np.array([])
                for distribute in distribution:
                    prob = cal_prob(bin_num*con.abundance, distribute)
                    probs = np.append(probs, prob)
                con.probs = probs
        
        #pdb.set_trace()
        # Maximization
        contig_win.clear()
        for group in cluster:
            for con in group:
                priors = assign_subcontig_by_posterior(con, cluster) # function
                con.priors = priors
                #con.priors = np.array([1,1,1,1,1])
                posters = [0]*k
                for i in range(len(priors)):
                    posters[i] = priors[i]*con.probs[i]
                posters = [x/sum(posters) for x in posters]
                #con.label = sample_label(posters) # function
                #pdb.set_trace()
                con.label = posters.index(max(posters))
                if not con.contig in contig_win:
                    contig_win[con.contig] = [con]
                else:
                    contig_win[con.contig].append(con)
                #new_cluster[max_idx].append(con)
        #new_cluster = adjust_label(contig_win, k) # function
        new_cluster = get_cluster(contig_win, k)
        #adjust_cluster = adjust_label(contig_win, k)
        adjust_cluster = adjust_label_by_prob(contig_win, k, distribution, bin_num)
        iteration += 1
    print 'Iteration:', iteration
    return new_cluster, distribution, adjust_cluster

ground_truth = ['HXB2', 'YU2', '89.6', 'NL43', 'JRCSF']
def evaluate(cluster):
    #TP,FP, TN, FN = 0,0,0,0
    TP = [0]*5
    idx = 0
    ground_truth_dict = {}
    for hp in ground_truth:
        ground_truth_dict[hp] = idx
        idx += 1
    contig_label = {}
    ground_ref_dict = {}
    for group in cluster:
        for con in group:
            ref = align_dict[con.contig]
            ref = ref.split('_')[0]
            pred_ref = ground_truth[con.label]
            if not con.contig in contig_label:
                contig_label[con.contig] = ref
                if pred_ref == ref:
                    TP[con.label]+=1
                if not ref in ground_ref_dict:
                    ground_ref_dict[ref] = 1
                else:
                    ground_ref_dict[ref] += 1

    return float(TP)/(TP+FP), float(TP)/len(contig_label)

ground_truth = ['HXB2', 'YU2', '89.6', 'NL43', 'JRCSF']
def evaluate_bp(cluster):
    pred_length = {}
    TP_length = {'HXB2':0, 'YU2':0, '89.6':0, 'NL43':0, 'JRCSF':0}
    ground_truth_length = {}
    idx = 0
    for group in cluster:
        for con in group:
            ref = align_dict[con.contig]
            ref = ref.split('_')[0]
            pred_ref = ground_truth[idx]
            if not ref in ground_truth_length:
                ground_truth_length[ref] = con.length
            else:
                ground_truth_length[ref] += con.length

            if not pred_ref in pred_length:
                pred_length[pred_ref] = con.length
            else:
                pred_length[pred_ref] += con.length 

            if ref==pred_ref:
                if ref in TP_length:
                    TP_length[ref] += con.length
                else:
                    TP_length[ref] = con.length
        idx += 1

    precision, recall = {}, {}
    for ref in ground_truth:
        #pdb.set_trace()
        precision[ref] = float(TP_length[ref])/pred_length[ref]
        recall[ref] = float(TP_length[ref])/ground_truth_length[ref]
    return precision, recall

fa_file = sys.argv[1]
loc_file=sys.argv[2]
vcf_file = sys.argv[3]
ref_loc_file = '/mnt/home/chenjiao/research/Project-virus/Evaluation_comparison/DT_meta/Contig_abundance/nodes_sequence_5vm.blastn'

fa_dict = read_fa(fa_file)
G = get_graph_from_loc(loc_file)
con_profile = read_vcf_profile(vcf_file, fa_dict)
align_dict = get_aligned_contigs(ref_loc_file)

flag = 2
f_out = open('windows.txt','w')
f2 = open('contigs_windows.txt', 'w')
f3 = open('windows_all.txt','w')
f_out.write('#contig\treference\tref_abundance\twin_start\twin_end\twin_abundance\twhole_abundance\twhole_coverage\n')
total, correct = 0,0
all_windows = []
win_con_dict = {}
contig_windows = []
for con in fa_dict:
    win_results, A, map_idx, align_locs = get_abundance_contig3(G, con_profile, con, align_dict)
    output_windows(G, win_results, A, map_idx, align_dict, ref_dict, f2)
    win_results = sorted(win_results, key = lambda win:(win.length, np.mean(win.coverage)), reverse=True)
    all_windows.extend(win_results)
    contig_windows.append(win_results[0])
    #if win_results[0].contig[0]=='19802|313':
    #    pdb.set_trace()
    win_con_dict[con] = win_results[0].abund_mean
    #win_results = sorted(win_results, key = lambda win:np.mean(win.coverage), reverse=True)
    for win in win_results:
        ref = align_dict[win.contig[0]]
        ref = ref.split('_')[0]
        f_out.write(win.contig[0]+'\t'+align_dict[win.contig[0]]+'\t'+str(ref_dict[ref])+'\t'+str(win.start)+'\t'+str(win.end)+'\t'+str(round(win.get_abund_mean(),4))+'\t'+str(round(np.mean(A),4))+'\t'+str(round(np.mean(win.coverage),4))+'\n')
    total+=1
    con_ref = get_ref(ref_dict, win_results[0].get_abund_mean())
    if re.search(con_ref, align_dict[con]):
        correct +=1
f_out.write("Accuracy:\t"+str(round(float(correct)/total,4))+'\n')
f_out.close()
f2.close()

#all_windows = sorted(all_windows, key = lambda win:(win.length, np.mean(win.coverage)), reverse=True)
all_windows = sorted(all_windows, key = lambda win:(np.sum(win.coverage[0][win.start-1:win.end]), win.length), reverse=True)

#f3.write('#contig\treference\tref_abundance\twin_start\twin_end\twin_abundance\twin_coverage\n')
sort_ref = ['JRCSF', 'NL43', '89.6', 'YU2', 'HXB2']
for win in all_windows:
    con = win.contig[0]
    ref = align_dict[win.contig[0]]
    ref = ref.split('_')[0]
    tmp_results = []
    for i in range(len(win.aligns)):
        con = win.aligns[i]
        win_dep = np.sum(win.profile[i][win.start-1:win.end])
        if win_dep>1.0:
            ref = align_dict[con]
            ref = ref.split('_')[0]
            abund_mean = np.mean(win.abundance[i][win.start-1:win.end])
            cov_sum = np.sum(win.coverage[i][win.start-1:win.end])
            tmp_results.append([con, ref,abund_mean, cov_sum])
            #f3.write(win.aligns[i]+'\t'+str(round(np.mean(win.abundance[i][win.start-1:win.end]),4))+'\t'+ref+'\t'+str(ref_dict[ref])+'\n')
    tmp_results = sorted(tmp_results, key = lambda win: win[2], reverse=True)
    
    flag = 0
    if len(tmp_results)==len(sort_ref) and win.length>10:
        for i in range(len(tmp_results)):
            ref = tmp_results[i][1]
            if ref!=sort_ref[i]:
                flag=1

    if len(tmp_results)==len(sort_ref) and win.length>10 and not flag:
        f3.write('>Window\treference_contig: '+con+'\t'+str(win.depth)+'\t'+str(win.start)+'\t'+str(win.end)+'\t'+ref+'\t'+str(ref_dict[ref])+'\tCorrect'+'\n')
    else:
        f3.write('>Window\treference_contig: '+con+'\t'+str(win.depth)+'\t'+str(win.start)+'\t'+str(win.end)+'\t'+ref+'\t'+str(ref_dict[ref])+'\n')

    for res in tmp_results:
        f3.write(res[0]+'\t'+str(round(res[2],4))+'\t'+str(round(res[3], 2))+'\t'+res[1]+'\t'+str(ref_dict[res[1]])+'\n')
f3.close()

## clustering
windows_list = []
for win in all_windows:
    con = win.contig[0]
    ref = align_dict[win.contig[0]]
    ref = ref.split('_')[0]
    tmp_results = []
    start = win.start
    end = win.end
    for i in range(len(win.aligns)):
        con = win.aligns[i]
        win_dep = np.sum(win.profile[i][win.start-1:win.end])
        start1, end1, start2, end2 = win.align_locs[i]
        align_win_start = start - start1 + start2 # position on the contig for the window
        align_win_end = end2 - (end1 - end) 
        if win_dep>1.0:
            ref = align_dict[con]
            ref = ref.split('_')[0]
            abund_mean = np.mean(win.abundance[i][win.start-1:win.end])
            cov_sum = np.sum(win.coverage[i][win.start-1:win.end])
            tmp_subContig = subContig(con, align_win_start, align_win_end, ref)
            tmp_subContig.coverage = win.coverage[i][win.start-1:win.end]
            tmp_subContig.abundance = win.abundance[i][win.start-1:win.end]
            tmp_subContig.ref = ref
            tmp_results.append(tmp_subContig)
    if tmp_results:
        windows_list.append(tmp_results)

print len(windows_list)

windows_list = merge_duplicate_windows(windows_list)
print len(windows_list)

f_out = open('Duplicate_windows.txt', 'w')
# the windows that are considered same are saved in one group
non_dup_list = []
idx =0
for group in windows_list:
    longest_win = group[0]
    longest_con = longest_win[0].contig
    if len(longest_win)==0:
        pdb.set_trace()
    longest_len = longest_win[0].end-longest_win[0].start+1
    for win in group[1:]:
        for con in win:
            if con.contig == longest_con:
                win_len = con.end-con.start+1
                if longest_len<win_len:
                    longest_len = win_len
                    longest_win = win
    non_dup_list.append(longest_win)
    f_out.write('#Group\n')
    for win in group:
        f_out.write('>Window\n')
        for con in win:
            ref = align_dict[con.contig]
            ref = ref.split('_')[0]
            f_out.write(con.contig+'\t'+str(con.start)+'\t'+str(con.end)+'\t'+str(np.mean(con.abundance))+'\t'+ref+'\t'+str(ref_dict[ref])+'\n')
    idx+=1
f_out.close()

f_out = open('Non_duplicate_windows.txt', 'w')
for win in non_dup_list:
    f_out.write('>Window\n')
    for con in win:
        ref = align_dict[con.contig]
        ref = ref.split('_')[0]
        f_out.write(con.contig+'\t'+str(con.start)+'\t'+str(con.end)+'\t'+str(np.mean(con.abundance))+'\t'+ref+'\t'+str(ref_dict[ref])+'\n')
f_out.close()

f_out = open('Non_duplicate_subcontigs.txt', 'w')
all_subcontigs = []
dup_idx = []
for win in non_dup_list:
    for con in win:
        all_subcontigs.append(con)
for i in range(len(all_subcontigs)):
    con1 = all_subcontigs[i]
    for j in range(i+1, len(all_subcontigs)):
        con2 = all_subcontigs[j]
        if con1.contig == con2.contig:
            if con1.start<=con2.start and con1.end>=con2.end:
                dup_idx.append(j)
            elif con1.start>=con2.start and con1.end<=con2.end:
                dup_idx.append(i)

for idx in sorted(list(set(dup_idx)), reverse=True):
    all_subcontigs.pop(idx)


bin_num = 100
cluster, distribution, adjust_cluster = EM_cluster_gibbs(non_dup_list, 5, bin_num)
#cluster_1000, distribution, adjust_cluster_1000 = EM_cluster_gibbs(non_dup_list, 5, 1000)
#print is_same_cluster(cluster, cluster_1000)
#cluster_format = format_cluster(cluster)
#cluster_1000_format = format_cluster(cluster_1000)
#pdb.set_trace()
f_out = open('EM_clusters.txt', 'w')
idx = -1
#for group in adjust_cluster:
for group in cluster:
    idx+=1
    group1 = cluster[idx]
    abund_array = np.array([])
    for con in group1:
        abund_array = np.append(abund_array, con.abundance)

    f_out.write('>Group '+str(idx+1)+'\t'+str(np.mean(abund_array))+'\n')
    for con in group:
        ref = align_dict[con.contig]
        ref = ref.split('_')[0]
        probs = []
        cov_array = bin_num*con.abundance
        #cov_array = con.coverage
        for distribute in distribution:
            prob = cal_prob(cov_array, distribute)
            probs.append(prob)
        probs = [str(x) for x in probs]
        probs = '\t'.join(probs)
        priors = [str(x) for x in con.priors]
        priors = '\t'.join(priors)
        #f_out.write(con.contig+'\t'+str(con.start)+'\t'+str(con.end)+'\t'+str(np.mean(con.abundance))+'\t'+ref+'\t'+str(ref_dict[ref])+'\n')
        f_out.write(con.contig+'\t'+str(con.start)+'\t'+str(con.end)+'\t'+str(np.mean(con.abundance))+'\t'+ref+'\t'+str(ref_dict[ref])+'\t+\t'+probs+'\t+\t'+priors+'\n')
#precision, recall = evaluate(adjust_cluster)
#f_out.write("Evaluation\nPrecision: "+str(precision)+"\tRecall: "+str(recall)+'\n')
#"""
f_out.write("Evaluation\n")
precision, recall = evaluate_bp(cluster)
for ref in ground_truth:
    f_out.write(ref+'\t'+str(precision[ref])+'\t'+str(recall[ref])+'\n')
f_out.write("Evaluation\n")
#"""


f_out.close()

#"""
windows_5_list = [x for x in non_dup_list if len(x)==5]

f = open('Nearest_pair.txt', 'w')
cluster = cluster_nearest_window2(windows_5_list)
f.close()
#"""

"""
f_out = open('windows_clusters.txt', 'w')
win_clusters = cluster_windows(all_windows)
idx = 0
for win in win_clusters:
    idx+=1
    out_line=''
    for i in range(len(win.contig)):
        out_line+=win.contig[i]+'\t'+str(round(win_con_dict[win.contig[i]],4))+'\t'
    f_out.write('Cluster_'+str(idx)+'\n'+out_line+'\n')
    
    truth = []
    for con in win.contig:
        ref = align_dict[con].split('_')[0]
        truth.append(align_dict[con])
        truth.append(str(ref_dict[ref]))
    f_out.write('\t'.join(truth)+'\n\n')
f_out.close()
"""


#f_out = open('windows_clusters_k_means.txt', 'w')
#cluster_windows_k(contig_windows, 5, align_dict, ref_dict, f_out)
