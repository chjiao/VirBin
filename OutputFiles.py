
def output_windows(G, win_results, A, map_idx, align_dict, ref_dict, f_out):
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