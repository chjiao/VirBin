from Align import *
import numpy as np


def read_fa(fa_file):
    """
    # key: title; value: sequence
    :param fa_file:
    :return:
    """

    seq_dict = {}
    title = ''
    seq = ''
    with open(fa_file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                if seq:
                    seq_dict[title] = seq
                title = line[1:-1]
                title = title.split()[0]
                seq = ''
            else:
                seq+=line.strip()
        if seq:
            seq_dict[title] = seq
    return seq_dict

def read_vcf_profile(vcf_file, contig_dict):
    con_profile_dict = {}
    for con in contig_dict:
        con_profile = np.zeros(len(contig_dict[con]))
        con_profile_dict[con] = con_profile

    with open(vcf_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            lmap=line.strip().split('\t')
            con,pos,info = lmap[0],lmap[1],lmap[7]
            m=re.search('DP=(\d+);',info)
            depth=int(m.group(1))
            con_profile_dict[con][int(pos)-1] = depth
    return con_profile_dict


def get_aligned_contigs(loc_file):
    """
    Get contigs aligned to references
    :param loc_file:
    :return:
    """

    lineno = 0
    simi_dict = {}
    align_dict = {}
    with open(loc_file, 'r') as f:
        for line in f:
            lineno += 1
            if lineno % 5 == 1:
                lmap1 = line.strip().split()
                con, ref, similarity = lmap1[0], lmap1[1], float(lmap1[2])
            elif lineno % 5 == 2:
                lmap2 = line.strip().split()
            elif lineno % 5 == 3:
                lmap3 = line.strip().split()
                con_len, align_start, align_end = int(lmap2[1]), int(lmap2[2]), int(lmap2[3])
                if align_start > align_end:
                    pdb.set_trace()
                align_len = abs(align_end - align_start + 1)
                if int(lmap2[1]) < 500 or float(align_len) / con_len < 0.5:
                    continue

                ## only look at contigs that are fully aligned
                if not con in simi_dict and float(align_len) / con_len > 0.7:
                    simi_dict[con] = similarity
                    align_dict[con] = ref
                elif align_len == con_len and simi_dict[con] < similarity:
                    simi_dict[con] = similarity
                    align_dict[con] = ref
    return align_dict

