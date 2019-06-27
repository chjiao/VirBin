import sys,pdb
from Graph import *
from OutputFiles import *
from EM import *
from Windows import *
import argparse
import numpy as np

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
        fa_dict = read_fa(args.input_contig)                                                             # return dict[title]=seq
        process_loc(args.contig_alignment, 'preprocessed.blastn')                                        
        G = get_graph_from_loc('preprocessed.blastn')
        con_profile = read_vcf_profile(args.vcf_file, fa_dict)
        align_dict = get_aligned_contigs(args.reference_alignment)

        if not args.bin_number:
            args.bin_number = 1000


    ## windows information on each contig
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



    file_in = open("windows_all.txt")
    cnt = 0
    tmp = []
    for line in file_in.readlines():
        if line[0] == ">":
            if cnt != 0:
                tmp.append(cnt)
            cnt=0
        else:
            cnt+=1

    counts = np.bincount(tmp)
    k = np.argmax(counts)
    print "number of cluter: " +str(k)
    
if __name__ == '__main__':
    main()
