
def DFS_transitive_reduction(G):
    for u in G.nodes():
        for v in G.successors(u):
            visited=set()
            stack=[v]
            while stack:
                vertex=stack.pop()
                if vertex not in visited:
                    if vertex in G[u] and vertex!=v:
                        G.remove_edge(u,vertex)
                    visited.add(vertex)
                    stack.extend(set(G.successors(vertex))-visited)

def get_graph_from_loc(loc_file):
    # get the directed graph from loc file and order the sequences
    G = nx.DiGraph()
    lineno = 0
    with open(loc_file,'r') as f:
        for line in f:
           lineno+=1
           if lineno%5==1:
               lmap1=line.strip().split()
           elif lineno%5==2:
               lmap2=line.strip().split()
           elif lineno%5==3:
               lmap3=line.strip().split()
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
               if not (align_end1-align_start1>100 and (align_end2-align_start2)>100):
                   continue
               if align_start1>align_start2:
                   G.add_edge(con1, con2, align = [align_start1, align_end1, con_len1, align_start2, align_end2, con_len2], profile=[align1, align2])
               elif align_start2>align_start1:
                   G.add_edge(con2, con1, align = [align_start2, align_end2, con_len2, align_start1, align_end1, con_len1], profile=[align2, align1])
               else:
                   if con_len1>con_len2:
                       G.add_edge(con1, con2, align = [align_start1, align_end1, con_len1, align_start2, align_end2, con_len2],profile=[align1, align2])
                   else:
                       G.add_edge(con2, con1, align = [align_start2, align_end2, con_len2, align_start1, align_end1, con_len1],profile=[align2, align1])
    #DFS_transitive_reduction(G)
    return G