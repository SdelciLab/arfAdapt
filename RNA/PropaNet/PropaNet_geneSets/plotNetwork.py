import os
import sys
import pandas as pd
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import plotly
import plotly.graph_objects as go
from plotly.subplots import make_subplots

def get_tfs():
    df = pd.read_csv('/home/upamanyu/PropaNet/data2/R6_temp_network.txt', sep=' ')
    tfs = list(df['TF'].unique())
    return tfs

def gen_bg_plot(G, pos):
    edge_x = []
    edge_y = []
    for edge in G.edges():
        x0, y0 = pos[edge[0]]
        x1, y1 = pos[edge[1]]
        edge_x.append(x0)
        edge_x.append(x1)
        edge_x.append(None)
        edge_y.append(y0)
        edge_y.append(y1)
        edge_y.append(None)

    edge_trace = go.Scatter(
        x=edge_x, y=edge_y,
        line=dict(width=0.3, color='gray'),
        opacity=0.2,
        hoverinfo='none',
        mode='lines')

    tfs = get_tfs()
    
    node_xy = {'tf':{'x':[], 'y':[], 'symbol':'diamond'}, 
               'no-tf':{'x':[], 'y':[], 'symbol':'circle'}}
    for node in G.nodes():
        x, y = pos[node]
        if node in tfs:
            k = 'tf'
        else:
            k = 'no-tf'
        node_xy[k]['x'].append(x)
        node_xy[k]['y'].append(y)

    node_traces = []
    for k, v in node_xy.items():
        node_trace = go.Scatter(
            x=v['x'], y=v['y'],
            mode='markers',
            opacity=0.2,
            marker=dict(
                color='gray',
                size=7,
                line_width=0,
                symbol=v['symbol'])
            )
        node_traces.append(node_trace)
    # node_text = []
    # for node in list(G.nodes()):
    #     node_text.append(node)
    # node_trace.text = node_text

    return node_traces, edge_trace

def gen_fg_plot(G, pos):
    edge_traces = []
    node_traces = []
    node_done = []
    tfs = get_tfs()
    for edge in G.edges():
        edge_x = []
        edge_y = []
        x0, y0 = pos[edge[0]]
        x1, y1 = pos[edge[1]]
        edge_x.append(x0)
        edge_x.append(x1)
        edge_x.append(None)
        edge_y.append(y0)
        edge_y.append(y1)
        edge_y.append(None)
        width = G.edges[edge[0], edge[1]]['weight']
        dash = 'solid'
        if width < 0:
            continue
        else:
            edge_traces.append(
                go.Scatter(
                    x=edge_x, y=edge_y,
                    line=dict(
                        width=1 + np.abs(width), 
                        color='black',
                        dash=dash),
                    opacity=0.4,
                    hoverinfo='none',
                    mode='lines')
            )

            node_xy = {'tf':{'x':[], 'y':[], 'symbol':'diamond'}, 
               'no-tf':{'x':[], 'y':[], 'symbol':'circle'}}
            node_x = []
            node_y = []
            for node in [edge[0], edge[1]]:
                if node in node_done:
                    continue
                
                if node in tfs:
                    k = 'tf'
                else:
                    k = 'no-tf'
                x, y = pos[node]
                node_xy[k]['x'].append(x)
                node_xy[k]['y'].append(y)
                
                node_done.append(node)

            for k, v in node_xy.items():
                if len(v['x']) != 0:
                    node_traces.append(go.Scatter(
                        x=v['x'], y=v['y'],
                        mode='markers',
                        hoverinfo='text',
                        opacity=0.5,
                        marker=dict(
                            color='blue',
                            size=7,
                            line_width=0.2,
                            symbol=v['symbol'])
                        ))
            
    return node_traces, edge_traces

def get_bg_net(path):
    G = nx.read_edgelist(path, delimiter='\t', 
        create_using=nx.DiGraph(), data=(("weight", float),))
    pos = nx.shell_layout(G)
    node_drop = []
    for adj in G.adjacency():
        if len(adj[1]) == 0:
            node_drop.append(adj[0])
    G.remove_nodes_from(node_drop)
    
    return G, pos

def get_fg_net(path, pos, glist):
    G = nx.read_edgelist(path, delimiter='\t', 
            create_using=nx.DiGraph(), data=(("weight", float),))
    G_sub = nx.DiGraph()
    n_sub = []
    e_sub = []
    for edge in G.edges():
        u = edge[0]
        v = edge[1]
        if (u in glist) or (v in glist):
            n_sub.append(u)
            n_sub.append(v)
            e_sub.append((u, v, G.edges[u, v]))
    G_sub.add_nodes_from(list(set(n_sub)))
    G_sub.add_edges_from(e_sub)
    G_sub_nodes = list(G_sub.nodes())
    pos_sub = {n:p for n,p in pos.items() if n in G_sub_nodes}

    return G_sub, pos_sub

def get_text_pos(x, y):
    if x < 0:
        if y < 0:
            return 'bottom left'
        elif y > 0:
            return 'top left'
        else:
            return 'left'
    elif x > 0:
        if y < 0:
            return 'bottom right'
        elif y > 0:
            return 'top right'
        else:
            return 'right'
    else:
        if y < 0:
            return 'bottom'
        else:
            return 'top'

def run_glist(path):
    with open(path, 'r') as f:
        glist = f.read().strip('\n').split('\n')
    glist = [g.strip(' ') for g in glist]

    traces = []

    G, pos = get_bg_net("nets/template_net.txt")
    bg_node_traces, bg_edge_trace = gen_bg_plot(G, pos)
    traces.append(bg_edge_trace)
    traces.extend(bg_node_traces)

    fig = go.Figure(layout=go.Layout(
                titlefont_size=16,
                showlegend=False,
                hovermode=False)
            )
    fig = make_subplots(1, 4, figure=fig, 
        vertical_spacing=0.02, horizontal_spacing=0.02)

    tfs = get_tfs()
    intersection = []
    for t in range(1, 5):
        G_sub, pos_sub = get_fg_net('nets/subnetwork.{}'.format(t), pos, glist)
        fg_node_traces, fg_edge_traces = gen_fg_plot(G_sub, pos_sub)
        traces.extend(fg_edge_traces)
        traces.extend(fg_node_traces)

        node_xy = {'tf':{'x':[], 'y':[], 'node_text':[], 'node_text_pos':[], 'symbol':'diamond'}, 
               'no-tf':{'x':[], 'y':[], 'node_text':[], 'node_text_pos':[], 'symbol':'circle'}}
        for node in G_sub.nodes():
            x, y = pos[node]
            if node in glist:
                if node in tfs:
                    k = 'tf'
                else:
                    k = 'no-tf'
                node_xy[k]['node_text'].append(node)
                node_xy[k]['node_text_pos'].append(get_text_pos(x, y))
                node_xy[k]['x'].append(x)
                node_xy[k]['y'].append(y)
                intersection.append(node)

        for k, v in node_xy.items():
            node_trace = go.Scatter(
                x=v['x'], y=v['y'],
                mode='markers+text',
                opacity=1,
                textfont_size=8,
                marker=dict(
                    color='Red',
                    size=12,
                    line_width=0.2,
                    symbol=v['symbol'])
                )
            node_trace.text = v['node_text']
            node_trace.textposition = v['node_text_pos']
            traces.append(node_trace)

        # fig.add_traces(data=traces, rows=1 + (t-1)//2, cols=1 + (t-1)%2,)
        fig.add_traces(data=traces, rows=1, cols=t,)
    
    fig.update_xaxes(showticklabels=False, showgrid=False, zeroline=False)
    fig.update_yaxes(showticklabels=False, showgrid=False, zeroline=False)
    fig.update_layout(width=2500, height=600, margin=dict(t=5, b=5, r=5, l=5))
    
    plotly.io.write_image(fig, path.split('.')[0]+'.pdf', format='pdf')
    
    with open(path.split('.')[0]+'_intersec.txt', 'w') as f:
        f.write('\n'.join(list(set(intersection))))
    # fig.update_layout(plot_bgcolor='white')
    # plotly.io.write_image(fig, 'base_network.png')
    # fig.show()
    
if __name__ == '__main__':
    # get_tfs()

    try:
        run_glist(sys.argv[1])
    except:
        for path in os.listdir('glists_withTFShapes'):
            if path.endswith('txt'):
                print('glists_withTFShapes/'+path)
                run_glist('glists_withTFShapes/'+path)
                # break