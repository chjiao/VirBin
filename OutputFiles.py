import numpy as np
import operator
from EM import cal_prob, evaluate_bp, evaluate_bp_whole


def output_windows_by_contig(G, win_results, A, map_idx, align_dict, ref_dict, f_out):
    idx_map = {}
    for con, idx in map_idx.items():
        idx_map[idx] = con
    H = A.shape[0]
    win_results = sorted(win_results, key = lambda win:(win.length, np.mean(win.coverage)), reverse=True)
    f_out.write('>Windows on contig: '+win_results[0].contig[0]+'\n')
    idx = 0
    for win in win_results:
        idx+=1
        f_out.write('Window'+str(idx)+'\t'+str(win.start)+'\t'+str(win.end)+'\t'+str(win.depth)+'\n')
        for i in range(H):
            con = idx_map[i]
            win_mean = np.mean(A[i][win.start-1:win.end])
            win_sum = np.sum(A[i][win.start-1:win.end])
            if win_sum>1.0:
                ref = align_dict[con]
                ref = ref.split('_')[0]
                f_out.write(con+'\t'+str(round(win_mean, 4))+'\t'+align_dict[con]+'\t'+str(ref_dict[ref])+'\n')


def output_windows(all_windows, align_dict, ref_dict, file_name):
    ref_sort = sorted(ref_dict.items(), key=operator.itemgetter(1))
    sort_ref = [x[0] for x in ref_sort]

    f_out = open(file_name, 'w')
    for win in all_windows:
        this_con = win.contig[0]
        this_ref = align_dict[win.contig[0]]
        this_ref = this_ref.split('_')[0]
        tmp_results = []
        for i in range(len(win.aligns)):
            con = win.aligns[i]
            win_dep = np.sum(win.profile[i][win.start - 1:win.end])
            if win_dep > 1.0:
                ref = align_dict[con]
                ref = ref.split('_')[0]
                abund_mean = np.mean(win.abundance[i][win.start - 1:win.end])
                cov_sum = np.sum(win.coverage[i][win.start - 1:win.end])
                tmp_results.append([con, ref, abund_mean, cov_sum])
        tmp_results = sorted(tmp_results, key=lambda win: win[2], reverse=True)  # desending

        flag = 0
        if len(tmp_results) == len(sort_ref) and win.length > 10:
            for i in range(len(tmp_results)):
                ref = tmp_results[i][1]
                if ref != sort_ref[i]:
                    flag = 1

        if len(tmp_results) == len(sort_ref) and win.length > 10 and not flag:
            f_out.write(
                '>Window\treference_contig: ' + this_con + '\t' + str(win.depth) + '\t' + str(win.start) + '\t' + str(
                    win.end) + '\t' + this_ref + '\t' + str(ref_dict[ref]) + '\tCorrect' + '\n')
        else:
            f_out.write(
                '>Window\treference_contig: ' + this_con + '\t' + str(win.depth) + '\t' + str(win.start) + '\t' + str(
                    win.end) + '\t' + this_ref + '\t' + str(ref_dict[ref]) + '\n')

        for res in tmp_results:
            f_out.write(res[0] + '\t' + str(round(res[2], 4)) + '\t' + str(round(res[3], 2)) + '\t' + res[1] + '\t' + str(
                ref_dict[res[1]]) + '\n')
    f_out.close()


def output_groups_windows(group_windows_list, align_dict, ref_dict, file_name):
    f_out = open(file_name, 'w')
    for group in group_windows_list:
        f_out.write('#Group\n')
        for win in group:
            f_out.write('>Window\n')
            for con in win:
                ref = align_dict[con.contig]
                ref = ref.split('_')[0]
                f_out.write(con.contig + '\t' + str(con.start) + '\t' + str(con.end) + '\t' + str(
                    np.mean(con.abundance)) + '\t' + ref + '\t' + str(ref_dict[ref]) + '\n')
    f_out.close()


def output_non_dup_windows(non_dup_list, align_dict, ref_dict, file_name):
    f_out = open(file_name, 'w')
    for win in non_dup_list:
        f_out.write('>Window\n')
        for con in win:
            ref = align_dict[con.contig]
            ref = ref.split('_')[0]
            f_out.write(con.contig + '\t' + str(con.start) + '\t' + str(con.end) + '\t' + str(
                np.mean(con.abundance)) + '\t' + ref + '\t' + str(ref_dict[ref]) + '\n')
    f_out.close()


def output_cluster(cluster, distribution, align_dict, ref_dict, bin_num, file_name):

    ref_sort = sorted(ref_dict.items(), key=operator.itemgetter(1))
    ground_truth = [x[0] for x in ref_sort]

    f_out = open(file_name, 'w')
    idx = -1
    for group in cluster:
        idx += 1
        group1 = cluster[idx]
        abund_array = np.array([])
        for con in group1:
            abund_array = np.append(abund_array, con.abundance)

        f_out.write('>Group ' + str(idx + 1) + '\t' + str(np.mean(abund_array)) + '\n')
        for con in group:
            ref = align_dict[con.contig]
            ref = ref.split('_')[0]
            probs = []
            cov_array = bin_num * con.abundance
            # cov_array = con.coverage
            for distribute in distribution:
                prob = cal_prob(cov_array, distribute)
                probs.append(prob)
            probs = [str(x) for x in probs]
            probs = '\t'.join(probs)
            priors = [str(x) for x in con.priors]
            priors = '\t'.join(priors)
            f_out.write(con.contig + '\t' + str(con.start) + '\t' + str(con.end) + '\t' + str(
                np.mean(con.abundance)) + '\t' + ref + '\t' + str(
                ref_dict[ref]) + '\t+\t' + probs + '\t+\t' + priors + '\n')
    f_out.write("Evaluation\n")
    precision, recall = evaluate_bp(cluster, align_dict, ref_dict)
    accuracy = evaluate_bp_whole(cluster, align_dict, ref_dict)
    for ref in ground_truth:
        f_out.write(ref + '\t' + str(precision[ref]) + '\t' + str(recall[ref]) + '\n')
    f_out.write("Evaluation\n")

    f_out.write("Accuracy\t" + str(accuracy) + '\n')
    f_out.close()

