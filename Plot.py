def plot_graph(G, figname):
    G_plot=nx.drawing.nx_agraph.to_agraph(G)
    G_plot.draw(figname,prog='dot')


def plot_alignment_contig(G, fa_dict, con_profile, align_dict, vertex):
    A, map_idx, B = get_abundance_contig(G, con_profile, con, align_dict, 0)
    con_start = 1
    height = 0
    plt.figure(figsize=(30, 20))
    con_len = len(fa_dict[vertex])
    len_shift = con_len / 20
    con_end = con_start + con_len - 1
    plt.plot((con_start, con_end), (height, height), 'g-', linewidth=4.0)
    plt.plot(range(con_start, con_end + 1), A[map_idx[vertex]] + height, 'r.')
    plt.text(con_start - len_shift, height + 0.1, str(con_start), horizontalalignment='center',
             verticalalignment='center', fontsize=12, color='g')  # left coordinate
    plt.text(con_end + len_shift, height + 0.1, str(con_end), horizontalalignment='center', verticalalignment='center',
             fontsize=12)  # right coordinate
    # pdb.set_trace()
    plt.text((con_start + con_end) / 2 - 50, height + 0.2, vertex + '(' + align_dict[vertex] + ')',
             horizontalalignment='center', verticalalignment='center', fontsize=12, color='b', weight='bold')
    height += 1
    for succ_node in G.successors(vertex):
        start1, end1, con_len1, start2, end2, con_len2 = G[vertex][succ_node]['align']
        con_align_start1 = con_start + start1 - 1
        con_align_end1 = end1 - start1 + con_align_start1
        con_start2 = con_align_start1 + (1 - start2)
        con_end2 = con_start2 + con_len2 - 1
        plt.plot((con_align_start1, con_align_end1), (height, height), 'b-', linewidth=4.0)
        plt.plot(range(con_align_start1, con_align_end1 + 1),
                 A[map_idx[succ_node]][con_align_start1 - 1:con_align_end1] + height, 'r.')
        plt.text(con_align_start1 - len_shift, height + 0.1, str(start2) + '(' + str(con_align_start1) + ')',
                 horizontalalignment='center', verticalalignment='center', fontsize=12, color='g')  # left coordinate
        plt.text(con_align_end1 + len_shift, height + 0.1, str(end2) + '(' + str(con_align_end1) + ')',
                 horizontalalignment='center', verticalalignment='center', fontsize=12)  # right coordinate
        plt.text((con_align_start1 + con_align_end1) / 2 - 50, height + 0.2,
                 succ_node + '(' + align_dict[succ_node] + ')', horizontalalignment='center',
                 verticalalignment='center', fontsize=12, color='b', weight='bold')
        height += 1

    for pre_node in G.predecessors(vertex):
        start2, end2, con_len2, start1, end1, con_len1 = G[pre_node][vertex]['align']
        con_align_start1 = con_start + start1 - 1
        con_align_end1 = end1 - start1 + con_align_start1
        con_start2 = con_align_start1 + (1 - start2)
        con_end2 = con_start2 + con_len2 - 1
        # plt.plot((con_start2, con_end2), (height, height), 'y-', linewidth=3.0)
        plt.plot((con_align_start1, con_align_end1), (height, height), 'b-', linewidth=4.0)
        plt.plot(range(con_align_start1, con_align_end1 + 1),
                 A[map_idx[pre_node]][con_align_start1 - 1:con_align_end1] + height, 'r.')
        plt.text(con_align_start1 - len_shift, height + 0.1, str(start2) + '(' + str(con_align_start1) + ')',
                 horizontalalignment='center', verticalalignment='center', fontsize=12, color='g')  # left coordinate
        plt.text(con_align_end1 + len_shift, height + 0.1, str(end2) + '(' + str(con_align_end1) + ')',
                 horizontalalignment='center', verticalalignment='center', fontsize=12)  # right coordinate
        plt.text((con_align_start1 + con_align_end1) / 2 - 50, height + 0.2,
                 pre_node + '(' + align_dict[pre_node] + ')', horizontalalignment='center', verticalalignment='center',
                 fontsize=12, color='b', weight='bold')
        height += 1
    figname = vertex + '_align.png'
    plt.savefig(figname, format='png', dpi=300)
    plt.close()


def plot_alignment(G, fa_dict, con_profile):
    subgraph_num = nx.number_weakly_connected_components(G)
    index = 1
    for G_sub in nx.weakly_connected_component_subgraphs(G):
        idx = 0
        height = 0
        plotted_dict = {}
        plt.figure(figsize=(30, 20))
        # A = get_abundance(G_sub, con_profile)
        # H, W = np.shape(A)
        plt.subplot(subgraph_num, 1, index)
        start_nodes = [N for N in G_sub if G_sub.in_degree(N) == 0]
        visited = set([])
        pos_start = 0
        queue = []
        put_start_nodes(G, queue, start_nodes)
        while (queue):
            vertex, con_start = queue.pop(0)
            visited.add(vertex)
            con_len = len(fa_dict[vertex])
            con_end = con_start + con_len
            plt.plot((con_start, con_end), (height, height), 'y-', linewidth=3.0)
            plt.text(con_start - 200, height + 0.1, str(con_start), horizontalalignment='center',
                     verticalalignment='center', fontsize=12, color='g')  # left coordinate
            plt.text(con_end + 200, height + 0.1, str(con_end), horizontalalignment='center',
                     verticalalignment='center', fontsize=12)  # right coordinate
            plt.text((con_start + con_end) / 2 - 50, height + 0.9, vertex, horizontalalignment='center',
                     verticalalignment='center', fontsize=12, color='b', weight='bold')
            for succ_node in G.successors(vertex):
                start1, end1, con_len1, start2, end2, con_len2 = G[vertex][succ_node]['align']
                con_align_start1 = con_start + start1
                con_align_end1 = end1 - start1 + con_align_start1
                plt.plot((con_align_start1, con_align_end1), (height, height), 'b-', linewidth=4.0)
                con_start2 = con_align_start1 + (1 - start2)
                con_end2 = con_start2 + con_len2 - 1
                if not succ_node in visited:
                    queue.append((succ_node, con_start2))
            height += 1
        index += 1
    figname = 'contig_alignments_all.png'
    plt.savefig(figname, format='png', dpi=300)
    plt.close()