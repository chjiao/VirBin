import re, sys, pdb
import ReadFiles
import operator
import numpy as np


def get_predict_ref(marker_gene, align_dict):
    # get the reference for marker genes
    result = []
    with open(marker_gene, 'r') as f:
        for line in f:
            gene = line.strip()
            assert gene in align_dict, "Error: reference not found "+ gene
            result.append(align_dict[gene])
    return result

def get_predict_ref_from_summary(cluster, align_dict, contig_dict):
    result = []
    for group in cluster:
        pred_dict = {}
        for con in group:
            ref = align_dict[con]
            con_seq = contig_dict[con]
            pred_dict[ref] = pred_dict.get(ref, 0) + len(con_seq)
        ground_truth_ref = max(pred_dict.iteritems(), key=operator.itemgetter(1))[0]
        result.append(ground_truth_ref)
    return result

def get_predict_ref_by_ranking(cluster, haplotype_list, abund_array):
    result = []
    indexes = [sorted(abund_array, reverse=True).index(x) for x in abund_array]
    #pdb.set_trace()
    for idx in indexes:
        result.append(haplotype_list[idx])
    return result


def read_maxbin_results(maxbin_summary):
    # read maxbin summary file to a list of groups
    result = []
    con_num = 0
    group = []
    cluster_abund = 0
    cluster_abund_array = []
    with open(maxbin_summary, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith("Bin ["):
                if group:
                    result.append(group)
                    con_num += len(group)
                    group = []
                cluster_abund = float(line.strip().split()[2])
                cluster_abund_array.append(cluster_abund)
            elif line.startswith("Bins without"):
                break
            else:
                if line:
                    group.append(line.split()[0])
        if group:
            result.append(group)
            con_num += len(group)
    return result, con_num, cluster_abund_array

def evaluate_bp(cluster, align_dict, contig_dict, ground_truth):
    # cluster: list of groups
    # align_dict: key: contig, value: reference haplotype
    # contig_dict: key: contig, value: contig sequence
    haplotypes = list(set(align_dict.values()))
    
    pred_length = {}
    TP_length = dict(zip(haplotypes, [0]*len(ground_truth)))
    ground_truth_length = {}

    for con in contig_dict.keys():
        ref = align_dict[con]
        ground_truth_length[ref] = ground_truth_length.get(ref, 0) + len(contig_dict[con])
        
    idx = 0
    for group in cluster:
        for con in group:
            ref = align_dict[con]
            #ref = ref.split('_')[0]
            predict_ref = ground_truth[idx]
            #ground_truth_length[ref] = ground_truth_length.get(ref, 0) + len(contig_dict[con])
            pred_length[predict_ref] = pred_length.get(predict_ref, 0) + len(contig_dict[con])

            if ref==predict_ref:
                TP_length[ref] = TP_length.get(ref, 0) + len(contig_dict[con])
        idx += 1

    precision, recall = {}, {}
    #for ref in ground_truth:
    for ref in list(set(align_dict.values())):
        if ref in pred_length and ref in ground_truth_length:
            precision[ref] = float(TP_length[ref])/pred_length[ref]
            recall[ref] = float(TP_length[ref])/ground_truth_length[ref]
        else:
            precision[ref] = 0
            recall[ref] = 0
    return precision, recall

def preprocess_align_dict(align_dict):
    result = {}
    for (key, value) in align_dict.items():
        result[key] = value.split('_')[0]
    return result

def main():
    maxbin_result = sys.argv[1]
    contig_file = sys.argv[2]
    align_ref = sys.argv[3]
    marker_gene_file = sys.argv[4]
    f_out = open(sys.argv[5], 'w')
    #haplotype_list = ['FJ061', 'FJ061-h1', 'FJ061-h2', 'FJ061-h3', 'FJ066']
    haplotype_list = ['JRCSF', 'NL43', '89.6', 'YU2', 'HXB2']

    clusters, con_num, abund_list = read_maxbin_results(maxbin_result)
    #pdb.set_trace()
    contig_dict = ReadFiles.read_fa(contig_file)
    align_dict = ReadFiles.get_aligned_contigs(align_ref)
    align_dict = preprocess_align_dict(align_dict)
    marker_references = get_predict_ref(marker_gene_file, align_dict)
    """
    if len(abund_list) == len(haplotype_list):
        marker_references = get_predict_ref_by_ranking(clusters, haplotype_list, abund_list)
    else:
        marker_references = get_predict_ref_from_summary(clusters, align_dict, contig_dict)
    """

    abund_list = [x/sum(abund_list) for x in abund_list]
    #abund_sort = sorted(abund_list)
    abund_str = [str(x) for x in abund_list]

    #pdb.set_trace()
    #f_out = open("maxbin_evaluation.txt", 'w')
    precision, recall = evaluate_bp(clusters, align_dict, contig_dict, marker_references)
    f_out.write('Evaluation\n')
    #for ref in marker_references:
    for ref in list(set(align_dict.values())):
        f_out.write(ref + '\t' + str(precision[ref]) + '\t' + str(recall[ref])+'\n')
        print(ref + '\t' + str(precision[ref]) + '\t' + str(recall[ref]))
    f_out.write('Evaluation\n')

    f_out.write('\n\n')
    f_out.write("Total contigs number: " + str(len(contig_dict))+'\n')
    unclassified_num = len(contig_dict) - con_num
    f_out.write("Unclassified contigs number: " + str(unclassified_num)
                + ' (' + str(100*float(unclassified_num)/len(contig_dict)) + '%)\n')
    f_out.write("Abundances:\n" + '\t'.join(marker_references) + '\n' + '\t'.join(abund_str) + '\n') 
    f_out.close()

main()

