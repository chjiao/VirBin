import pdb
import numpy as np
from lib import *

def get_subcontigs_from_windows(all_windows, align_dict):
    windows_list = []
    for win in all_windows:
        tmp_results = []
        start = win.start
        end = win.end
        for i in range(len(win.aligns)):
            con = win.aligns[i]
            win_dep = np.sum(win.profile[i][win.start - 1:win.end])
            start1, end1, start2, end2 = win.align_locs[i]
            align_win_start = start - start1 + start2  # position on the contig for the window
            align_win_end = end2 - (end1 - end)
            if win_dep > 1.0:
                ref = align_dict[con]
                ref = ref.split('_')[0]
                abund_mean = np.mean(win.abundance[i][win.start - 1:win.end])
                cov_sum = np.sum(win.coverage[i][win.start - 1:win.end])
                tmp_subContig = subContig(con, align_win_start, align_win_end, ref)
                tmp_subContig.coverage = win.coverage[i][win.start - 1:win.end]
                tmp_subContig.abundance = win.abundance[i][win.start - 1:win.end]
                tmp_subContig.ref = ref
                tmp_results.append(tmp_subContig)
        if tmp_results:
            windows_list.append(tmp_results)
    return windows_list

def same_window(win1, win2):
    """
    Same windows are defined as windows that have the same set of aligned contigs and
    the start position and end position are similar on the reference contig.
    :param win1:
    :param win2:
    :return: 1 if win1 is same to win2; 0 if different
    """
    flag = 0
    for con in win2:
        if con.contig==win1[0].contig:
            if abs(con.start - win1[0].start)<50 and abs(con.end - win1[0].end)<50:
                flag = 1
    return flag

def find_duplicate_windows(win_list):
    """
    Put the same windows in a list.
    :param win_list: a list of windows
    :return: a list of lists, each list contains windows that are considered same
    """
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
                remove_idx.append(i)
        result.append(tmp_results)
        for idx in sorted(remove_idx, reverse=True):
            tmp_win_list.pop(idx)
    return result

def get_longest_non_duplicate_windows(windows_list):
    """
    Get the longest window for each group of duplicate windows
    :param windows_list:
    :return:
    """
    non_dup_list = []
    for group in windows_list:
        longest_win = group[0]
        longest_con = longest_win[0].contig
        if len(longest_win) == 0:
            pdb.set_trace()
        longest_len = longest_win[0].end - longest_win[0].start + 1
        for win in group[1:]:
            for con in win:
                if con.contig == longest_con:
                    win_len = con.end - con.start + 1
                    if longest_len < win_len:
                        longest_len = win_len
                        longest_win = win
        non_dup_list.append(longest_win)
    return non_dup_list

def remove_duplicate_subcontigs(non_dup_list):
    all_subcontigs = []
    dup_idx = []
    for win in non_dup_list:
        for con in win:
            all_subcontigs.append(con)

    for i in range(len(all_subcontigs)):
        con1 = all_subcontigs[i]
        for j in range(i + 1, len(all_subcontigs)):
            con2 = all_subcontigs[j]
            if con1.contig == con2.contig:
                if con1.start <= con2.start and con1.end >= con2.end:
                    dup_idx.append(j)
                elif con1.start >= con2.start and con1.end <= con2.end:
                    dup_idx.append(i)

    for idx in sorted(list(set(dup_idx)), reverse=True):
        all_subcontigs.pop(idx)
    return all_subcontigs

