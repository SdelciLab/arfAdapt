import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import networkx as nx

G = nx.read_edgelist('/home/upamanyu/PropaNet/result/intermediate_results/PropaNet.nwk.t1.trim', data=(('weight',float),),create_using=nx.DiGraph())
# nx.draw_networkx(DEGnet, arrows=True, with_labels=False, node_size=10, alpha=0.5, width=0.2)
pos = nx.spring_layout(G)
elarge = [(u, v) for (u, v, d) in G.edges(data=True) if d["weight"] > 0.5]
esmall = [(u, v) for (u, v, d) in G.edges(data=True) if d["weight"] <= 0.5]
labels = nx.get_edge_attributes(G, 'weight')

nx.draw_networkx(G, pos, arrows=True, with_labels=False, node_size=10, alpha=0.5, width=0.2, edge_color='gray')
nx.draw_networkx_edges(G, pos, edgelist=elarge, width=1)
# nx.draw_networkx_edges(G, pos, edgelist=esmall, width=1, alpha=0.5, edge_color="b", style="dashed")
# plt.savefig(<wherever>)
plt.savefig("t1.png", dpi=200)