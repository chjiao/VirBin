import re, sys, pdb
import ReadFiles


def get_predict_ref(marker_gene, align_dict):
    # get the reference for marker genes
    result = []
    with open(marker_gene, 'r') as f:
        for line in f:
            gene = line.strip()
            assert gene in align_dict, "Error: reference not found "+ gene
            result.append(align_dict[gene])
    return result

def read_maxbin_results(maxbin_summary):
    # read maxbin summary file to a list of groups
    result = []
    con_num = 0
    group = []
    with open(maxbin_summary, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith("Bin ["):
                if group:
                    result.append(group)
                    con_num += len(group)
                    group = []
            elif line.startswith("Bins without"):
                break
            else:
                if line:
                    group.append(line.split()[0])
        if group:
            result.append(group)
            con_num += len(group)
    return result, con_num

def evaluate_bp(cluster, align_dict, contig_dict, ground_truth):
    # cluster: list of groups
    # align_dict: key: contig, value: reference haplotype
    # contig_dict: key: contig, value: contig sequence
    haplotypes = list(set(align_dict.values()))
    
    pred_length = {}
    TP_length = dict(zip(haplotypes, [0]*len(ground_truth)))
    ground_truth_length = {}
    
    idx = 0
    for group in cluster:
        for con in group:
            ref = align_dict[con]
            #ref = ref.split('_')[0]
            predict_ref = ground_truth[idx]
            ground_truth_length[ref] = ground_truth_length.get(ref, 0) + len(contig_dict[con])
            pred_length[predict_ref] = pred_length.get(predict_ref, 0) + len(contig_dict[con])

            if ref==predict_ref:
                TP_length[ref] = TP_length.get(ref, 0) + len(contig_dict[con])
        idx += 1

    precision, recall = {}, {}
    for ref in ground_truth:
        if ref in pred_length:
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

    clusters, con_num = read_maxbin_results(maxbin_result)
    contig_dict = ReadFiles.read_fa(contig_file)
    align_dict = ReadFiles.get_aligned_contigs(align_ref)
    align_dict = preprocess_align_dict(align_dict)
    marker_references = get_predict_ref(marker_gene_file, align_dict)

    pdb.set_trace()
    f_out = open("maxbin_evaluation.txt", 'w')
    precision, recall = evaluate_bp(clusters, align_dict, contig_dict, marker_references)
    for ref in marker_references:
        f_out.write(ref + '\t' + str(precision[ref]) + '\t' + str(recall[ref])+'\n')
        print(ref + '\t' + str(precision[ref]) + '\t' + str(recall[ref]))

    f_out.write('\n\n')
    f_out.write("Total contigs number: " + str(len(contig_dict))+'\n')
    unclassified_num = len(contig_dict) - con_num
    f_out.write("Unclassified contigs number: " + str(unclassified_num)
                + ' (' + str(100*float(unclassified_num)/len(contig_dict)) + '%)')
    f_out.close()

main()

