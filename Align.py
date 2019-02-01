import re, sys, pdb 

def process_loc(loc_file, out_file):
    # preprocess the alignment file to ensure one-one alignment between two contigs
    align_dict = {}
    lineno = 0
    with open(loc_file,'r') as f:
        for line in f:
           lineno+=1
           if lineno%5==1:
               lmap1=line.strip().split()
               line1 = line.strip()
               if len(lmap1)<3:
                   pdb.set_trace()
               con1, con2, simi = lmap1[0], lmap1[1], float(lmap1[2])
           elif lineno%5==2:
               lmap2=line.strip().split()
               line2 = line.strip()
           elif lineno%5==3:
               lmap3=line.strip().split()
               line3 = line.strip()
               # contig_1_13812  89.6_9669:      99.9
               # contig_1_1381   13812   1       9650
               # 89.6    9669    2       9651
               con1, con_len1, align_start1, align_end1 = lmap2[0], int(lmap2[1]), int(lmap2[2]), int(lmap2[3])
               con2, con_len2, align_start2, align_end2 = lmap3[0], int(lmap3[1]), int(lmap3[2]), int(lmap3[3])
           elif lineno%5==4:
               align1 = line.strip()
           elif lineno%5==0:
               align2 = line.strip()
               #if not (align_end1-align_start1>con_len1/2 or align_end2-align_start2>con_len2/2):
               align_len1 = align_end1-align_start1 + 1
               align_len2 = align_end2-align_start2 + 1
               if not (abs(align_end1-align_start1)>100 and abs(align_end2-align_start2)>100):
                   continue

               align_key = tuple(sorted([con1, con2]))
               if not align_key in align_dict:
                   align_dict[align_key] = [line1, line2, line3, align1, align2, simi, align_len1, align_len2]
               else:
                   prev_simi, prev_len1, prev_len2 = align_dict[align_key][5], align_dict[align_key][6], align_dict[align_key][7]
                   if simi >= prev_simi and (align_len1 >= prev_len1 or align_len2 >= prev_len2):
                       align_dict[align_key] = [line1, line2, line3, align1, align2, simi, align_len1, align_len2]
    
    f_out = open(out_file, 'w')
    for align_key in align_dict:
        align = align_dict[align_key]
        f_out.write('\n'.join(align[:5]))
        f_out.write('\n')
    f_out.close()

def get_reference(loc_file, out_file):
    # find the reference for contigs
    align_dict = {}
    lineno = 0
    with open(loc_file,'r') as f:
        for line in f:
           lineno+=1
           if lineno%5==1:
               lmap1=line.strip().split()
               line1 = line.strip()
               if len(lmap1)<3:
                   pdb.set_trace()
               con1, con2, simi = lmap1[0], lmap1[1], float(lmap1[2])
           elif lineno%5==2:
               lmap2=line.strip().split()
               line2 = line.strip()
               con1, con_len1, align_start1, align_end1 = lmap2[0], int(lmap2[1]), int(lmap2[2]), int(lmap2[3])
               align_len = align_end1 - align_start1 + 1
               
               if not con1 in align_dict:
                   align_dict[con1] = [con2, simi, align_len, con_len1]
               elif simi > align_dict[con1][1]:
                   align_dict[con1] = [con2, simi, align_len, con_len1]
    
    f_out = open(out_file, 'w')
    for align in align_dict:
        f_out.write(align+'\t'+align_dict[align][0]+'\t'+str(align_dict[align][1])+'\t'+str(align_dict[align][3])+'\t'+str(align_dict[align][2])+'\n')
    f_out.close()

def main():
    loc_file = sys.argv[1]
    out_file = sys.argv[2]
    process_loc(loc_file, out_file)
    #get_reference(loc_file, out_file)

main()
