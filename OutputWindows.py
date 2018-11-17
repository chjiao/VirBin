import sys,pdb
from ReadFiles import *
from Graph import *
from lib import *
from OutputFiles import *
from EM import *
from Align import *

ref_dict = {'FJ061':0.392531290961,
        'FJ061_h1':0.233236526254,
        'FJ061_h2':0.159603176643,
        'FJ061_h3':0.128119458223,
        'FJ066':0.0865095479195} # The ground truth abundance on each haplotype
sort_ref = ['FJ061', 'FJ061_h1', 'FJ061_h2','FJ061_h3', 'FJ066']


fa_file = sys.argv[1]  # contig file input
loc_file = sys.argv[2]
vcf_file = sys.argv[3]
ref_loc_file = '/mnt/home/chenjiao/research/Project-virus/Evaluation_comparison/DT_meta/Contig_abundance/HIV_simulate/contigs_simulate/hiv_simulate_5_contigs_ref.blastn'

fa_dict = read_fa(fa_file)
process_loc(loc_file, 'preprocessed.blastn')
G = get_graph_from_loc('preprocessed.blastn')
con_profile = read_vcf_profile(vcf_file, fa_dict)
align_dict = get_aligned_contigs(ref_loc_file)

flag = 2
f_out = open('windows.txt', 'w')
f2 = open('contigs_windows.txt', 'w')
f3 = open('windows_all.txt', 'w')
f_out.write('#contig\treference\tref_abundance\twin_start\twin_end\twin_abundance\twhole_abundance\twhole_coverage\n')
total, correct = 0, 0
all_windows = []
win_con_dict = {}
contig_windows = []
for con in fa_dict:
    win_results, A, map_idx, align_locs = get_abundance_on_contig(G, con_profile, con, fa_dict, align_dict)
    output_windows(G, win_results, A, map_idx, align_dict, ref_dict, f2)
    win_results = sorted(win_results, key=lambda win: (win.length, np.mean(win.coverage)), reverse=True)
    all_windows.extend(win_results)
    contig_windows.append(win_results[0])
    win_con_dict[con] = win_results[0].abund_mean

    for win in win_results:
        ref = align_dict[win.contig[0]]
        ref = ref.split('_')[0]
        f_out.write(win.contig[0] + '\t' + align_dict[win.contig[0]] + '\t' + str(ref_dict[ref]) + '\t' + str(win.start) + '\t' + str(win.end) + '\t' + str(round(win.get_abund_mean(), 4)) + '\t' + str(round(np.mean(A), 4)) + '\t' + str(round(np.mean(win.coverage), 4)) + '\n')
    #total += 1
    #con_ref = get_ref(ref_dict, win_results[0].get_abund_mean())
    #if re.search(con_ref, align_dict[con]):
    #    correct += 1
#f_out.write("Accuracy:\t" + str(round(float(correct) / total, 4)) + '\n')
f_out.close()
f2.close()

# all_windows = sorted(all_windows, key = lambda win:(win.length, np.mean(win.coverage)), reverse=True)
all_windows = sorted(all_windows, key=lambda win: (np.sum(win.coverage[0][win.start - 1:win.end]), win.length), reverse=True)

for win in all_windows:
    con = win.contig[0]
    ref = align_dict[win.contig[0]]
    ref = ref.split('_')[0]
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
            # f3.write(win.aligns[i]+'\t'+str(round(np.mean(win.abundance[i][win.start-1:win.end]),4))+'\t'+ref+'\t'+str(ref_dict[ref])+'\n')
    tmp_results = sorted(tmp_results, key=lambda win: win[2], reverse=True)

    flag = 0
    if len(tmp_results) == len(sort_ref) and win.length > 10:
        for i in range(len(tmp_results)):
            ref = tmp_results[i][1]
            if ref != sort_ref[i]:
                flag = 1

    if len(tmp_results) == len(sort_ref) and win.length > 10 and not flag:
        f3.write('>Window\treference_contig: ' + con + '\t' + str(win.depth) + '\t' + str(win.start) + '\t' + str(
            win.end) + '\t' + ref + '\t' + str(ref_dict[ref]) + '\tCorrect' + '\n')
    else:
        f3.write('>Window\treference_contig: ' + con + '\t' + str(win.depth) + '\t' + str(win.start) + '\t' + str(
            win.end) + '\t' + ref + '\t' + str(ref_dict[ref]) + '\n')

    for res in tmp_results:
        f3.write(res[0] + '\t' + str(round(res[2], 4)) + '\t' + str(round(res[3], 2)) + '\t' + res[1] + '\t' + str(
            ref_dict[res[1]]) + '\n')
f3.close()



