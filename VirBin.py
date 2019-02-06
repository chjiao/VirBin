import sys,pdb
from ReadFiles import *
from Graph import *
from Align import *
from lib import *
from OutputFiles import *
from EM import *
from Windows import *


ref_dict = {'FJ061':0.392531290961,
        'FJ061_h1':0.233236526254,
        'FJ061_h2':0.159603176643,
        'FJ061_h3':0.128119458223,
        'FJ066':0.0865095479195} # The ground truth abundance on each haplotype
sort_ref = ['FJ061', 'FJ061_h1', 'FJ061_h2', 'FJ061_h3', 'FJ066']  # desending


fa_file = sys.argv[1]  # contig file input
loc_file = sys.argv[2]
vcf_file = sys.argv[3]
ref_loc_file = sys.argv[4]

fa_dict = read_fa(fa_file)
process_loc(loc_file, 'preprocessed.blastn')
G = get_graph_from_loc('preprocessed.blastn')
con_profile = read_vcf_profile(vcf_file, fa_dict)
align_dict = get_aligned_contigs(ref_loc_file)

f_out = open('contigs_windows.txt', 'w')
total, correct = 0, 0
all_windows = []
#win_con_dict = {}
#contig_windows = []
for con in fa_dict:
    win_out =  get_abundance_on_contig(G, con_profile, con, fa_dict, align_dict)
    if win_out == None:
        continue
    win_results, A, map_idx, align_locs =  get_abundance_on_contig(G, con_profile, con, fa_dict, align_dict)

    if len(win_results) == 0:
        pdb.set_trace()
    output_windows_by_contig(G, win_results, A, map_idx, align_dict, ref_dict, f_out)
    win_results = sorted(win_results, key=lambda win: (win.length, np.mean(win.coverage)), reverse=True)
    all_windows.extend(win_results)
    #contig_windows.append(win_results[0])
    #win_con_dict[con] = win_results[0].abund_mean
f_out.close()

## output all windows
all_windows = sorted(all_windows, key=lambda win: (win.length, np.mean(win.coverage)), reverse=True)
output_windows(all_windows, align_dict, 'windows_all.txt')

## clustering
windows_list = get_subcontigs_from_windows(all_windows)
print len(windows_list)

groups_duplicate_windows_list = find_duplicate_windows(windows_list)
print len(groups_duplicate_windows_list)
#output_groups_windows(groups_duplicate_windows_list, align_dict, "Duplicate_windows.txt")

non_dup_list = get_longest_non_duplicate_windows(groups_duplicate_windows_list)
#output_non_dup_windows(non_dup_list, align_dict, 'Non_duplicate_windows.txt')

bin_num = 1000
cluster, distribution, adjust_cluster = EM_cluster_gibbs(non_dup_list, 5, bin_num)
output_cluster(cluster, 'EM_clusters.txt')



