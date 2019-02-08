import re,sys,subprocess,pdb
from itertools import groupby

# Calculate the sequence similarity using Water from EMBOSS
# Group fasta files before pairwise alignment
# Compare the sequences from multiple sequences in single fasta files
# 2018.03.14

fa_file=sys.argv[1]
out_file=sys.argv[2]
simi = 85

fa_dict1={}
title_list1=[]
contig_index=0
with open(fa_file,'r') as f:
    faiter=(x[1] for x in groupby(f,lambda line:line[0]=='>'))
    for header in faiter:
        contig_index+=1
        header=header.next().strip()
        seq="".join(s.strip() for s in faiter.next())

        if header.startswith('>'):
            title=header[1:].split()[0]
            #title='contig_'+str(contig_index)+'_'+str(len(seq))
            fa_dict1[title]=seq
            title_list1.append(title)
        else:
            print "Error!"
            pdb.set_trace()

f_out=open(out_file,'w')
##
for i in range(len(title_list1)-1):
    for j in range(i+1, len(title_list1)):
        f1=open('01.fa','w')
        f2=open('02.fa','w')
        seq1=fa_dict1[title_list1[i]]
        seq2=fa_dict1[title_list1[j]]
        f1.write('>'+title_list1[i]+'\n'+seq1)
        f2.write('>'+title_list1[j]+'\n'+seq2)
        f1.close()
        f2.close()
        blastn_command="blastn -query 01.fa -subject 02.fa -out temp.blastn"

        subprocess.call(blastn_command,shell=True)  #call will wait until the command completes
        #pdb.set_trace()

        loc_list=[]
        align_dict={}
        similarity=''
        flag=0
        with open("temp.blastn",'r') as f:
            for line in f:
                line=line.strip()
                if line.startswith('Identities'):
                    m=re.search('\((.*)%\)',line.split(',',1)[0])
                    #pdb.set_trace()
                    if m:
                        if flag:
                            # not the first time meeting with "Identities", process the last one
                            if float(similarity)>=simi:
                                f_out.write(title_list1[i]+'\t'+title_list1[j]+'_'+str(len(seq2))+':'+'\t'+similarity+'\n')
                                print i,j
                                #if i==56 and j==1:
                                #    pdb.set_trace()
                                f_out.write(title_list1[i]+'\t'+str(len(seq1))+'\t'+loc_list[0][1]+'\t'+loc_list[-2][2]+'\n')
                                f_out.write(title_list1[j]+'\t'+str(len(seq2))+'\t'+loc_list[1][1]+'\t'+loc_list[-1][2]+'\n')
                                f_out.write(align_dict[loc_list[0][0]]+'\n')
                                f_out.write(align_dict[loc_list[1][0]]+'\n')

                            loc_list=[]
                            align_dict={}
                        else:
                            flag=1
                        similarity=m.group(1)
                    else:
                        print "Similarity not found!"
                
                ## 
                lmap=line.strip().split()
                if len(lmap)==4:
                    m1=re.match('[atcgATCG-]',lmap[2])
                    if m1:
                        loc_list.append([lmap[0],lmap[1],lmap[3]])
                        if not lmap[0] in align_dict:
                            align_dict[lmap[0]]=lmap[2]
                        else:
                            align_dict[lmap[0]]+=lmap[2]

        ## process the last alignment
        if similarity and float(similarity)>=simi:
            f_out.write(title_list1[i]+'\t'+title_list1[j]+'_'+str(len(seq2))+':'+'\t'+similarity+'\n')
            f_out.write(title_list1[i]+'\t'+str(len(seq1))+'\t'+loc_list[0][1]+'\t'+loc_list[-2][2]+'\n')
            f_out.write(title_list1[j]+'\t'+str(len(seq2))+'\t'+loc_list[1][1]+'\t'+loc_list[-1][2]+'\n')
            f_out.write(align_dict[loc_list[0][0]]+'\n')
            f_out.write(align_dict[loc_list[1][0]]+'\n')
        #pdb.set_trace()

f_out.close()
subprocess.call('rm 01.fa',shell=True)
subprocess.call('rm 02.fa',shell=True)
subprocess.call('rm temp.blastn',shell=True)

