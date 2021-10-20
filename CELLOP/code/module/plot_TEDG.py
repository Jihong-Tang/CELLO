import matplotlib as mpl
import matplotlib.pyplot as plt
import networkx as nx

def plot_TEDG(df_edge, df_node):
    """
    """

    seed = 423
    G = nx.from_pandas_edgelist(
        df_edge, 'geneA', 'geneB', create_using=nx.DiGraph())
    pos = nx.spring_layout(G, seed=seed)
    #nx.spiral_layout

    #node_count = [float(i) for i in ]
    node_sizes = [15 + 2 * float(i) for i in df_node.Occurence.to_list()]
    edge_weight = [int(i)/2.5 for i in df_edge.weight.to_list()]
    nodes = nx.draw_networkx_nodes(G, pos, node_size=node_sizes, node_color=range(
        G.number_of_nodes()), cmap=plt.cm.plasma)
    edges = nx.draw_networkx_edges(
        G,
        pos,
        node_size=node_sizes,
        arrowstyle="->",
        arrowsize=10,
        edge_color="grey",
        width=edge_weight,
        alpha=0.7
    )
    labels = nx.draw_networkx_labels(G, pos, font_size=12)
    #nx.draw(G, pos, node_size = node_sizes, with_labels = True)
    plt.show()
    return(0)
