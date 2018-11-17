import res, sys, pdb

def get_predict_ref(marker_gene):

def evaluate_bp(cluster, align_dict, contig_dict, ground_truth):
    # cluster: list of groups
    # align_dict: key: contig, value: reference haplotype
    # contig_dict: key: contig, value: contig sequence
    haplotypes = list(set(align_dict.values))
    
    pred_length = {}
    TP_length = dict(zip(haplotypes, [0]*len(ground_truth)))
    ground_truth_length = {}
    
    idx = 0
    for group in cluster:
        for con in group:
            ref = align_dict[con.contig]
            ref = ref.split('_')[0]
            pred_ref = ground_truth[idx]
            ground_truth_length[ref] = ground_truth_length.get(ref, 0) + len(contig_dict[con])
            pred_length[pred_ref] = pred_length.get(pred_ref, 0) + len(contig_dict[con])

            if ref==pred_ref:
                TP_length[ref] = TP_length.get(ref, 0) + con.length
        idx += 1

    precision, recall = {}, {}
    for ref in ground_truth:
        #pdb.set_trace()
        precision[ref] = float(TP_length[ref])/pred_length[ref]
        recall[ref] = float(TP_length[ref])/ground_truth_length[ref]
    return precision, recall

def main():
    maxbin_result = sys.argv[1]
    contig_file = sys.argv[2]
    align_ref = sys.argv[3]
    marker_gene_file = sys.argv[4]
    
    contig_dict = read_fa(contig_file)
    align_dict = 

