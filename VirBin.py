import sys,pdb
from Graph import *
from OutputFiles import *
from EM import *
from Windows import *
import argparse


def main():

    ref_dict = {'FJ061': 0.392531290961,
                'FJ061-h1': 0.233236526254,
                'FJ061-h2': 0.159603176643,
                'FJ061-h3': 0.128119458223,
                'FJ066': 0.0865095479195}  # The ground truth abundance on each haplotype

    parser = argparse.ArgumentParser(description="VirBin: clustering of viral contigs for quasispecies")
    parser.add_argument('-contig', dest='input_contig', type=str, required=True,
                                help='input contig file in fasta format')
    parser.add_argument('-align', dest="contig_alignment", type=str, required=True,
                                help="alignment between contigs")
    parser.add_argument('-vcf', dest="vcf_file", type=str, required=True,
                                help="reads mapping profile in vcf format")
    parser.add_argument('-ref', dest="reference_alignment", type=str,
                                help="alignment on reference genomes")
    parser.add_argument('-bin', dest='bin_number', type=int, help="bin number for distribution of EM, default: 1000")

    if len(sys.argv[1:]) == 0:
        parser.print_help()
        parser.exit()

    args = parser.parse_args()

    if args.input_contig and args.contig_alignment and args.vcf_file and args.reference_alignment:
        fa_dict = read_fa(args.input_contig)
        process_loc(args.contig_alignment, 'preprocessed.blastn')
        G = get_graph_from_loc('preprocessed.blastn')
        con_profile = read_vcf_profile(args.vcf_file, fa_dict)
        align_dict = get_aligned_contigs(args.reference_alignment)

        if not args.bin_number:
            args.bin_number = 1000


    f_out = open('contigs_windows.txt', 'w')
    all_windows = []
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
    f_out.close()

    ## output all windows
    all_windows = sorted(all_windows, key=lambda win: (win.length, np.mean(win.coverage)), reverse=True)
    output_windows(all_windows, align_dict, ref_dict, 'windows_all.txt')

    ## clustering
    windows_list = get_subcontigs_from_windows(all_windows, align_dict)
    print len(windows_list)

    groups_duplicate_windows_list = find_duplicate_windows(windows_list)
    print len(groups_duplicate_windows_list)
    #output_groups_windows(groups_duplicate_windows_list, align_dict, "Duplicate_windows.txt")

    non_dup_list = get_longest_non_duplicate_windows(groups_duplicate_windows_list)
    #output_non_dup_windows(non_dup_list, align_dict, 'Non_duplicate_windows.txt')

    cluster, distribution, adjust_cluster = EM_cluster_gibbs(non_dup_list, 5, args.bin_number)
    output_cluster(cluster, distribution, align_dict, ref_dict, args.bin_number, 'EM_clusters.txt')

if __name__ == '__main__':
    main()


