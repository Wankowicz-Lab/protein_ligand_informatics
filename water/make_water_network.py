import networkx as nx
import json as json
from collections import Counter

import numpy as np

from scipy import stats

def get_kinase_waters():
    with open('/.../allNonKinaseWaters.json', 'r') as f: # [ { PDB, waters: [] } ]
        kinase_waters = json.load(f)
    return kinase_waters

def get_distance(water1, water2):
    return ((water1['x'] - water2['x']) ** 2 + (water1['y'] - water2['y']) ** 2 + (water1['z'] - water2['z']) ** 2) ** 0.5


def calculate_power_loss_scaling(graph):
    degrees = [d for _, d in graph.degree()]
    
    if len(degrees) < 10:
        return 1.0
    
    degree_counts = Counter(degrees)
    x = np.array(list(degree_counts.keys()))
    y = np.array(list(degree_counts.values()))
    
    if len(x) < 3:
        return 1.0
    
    log_x = np.log10(x + 1) # +1 for degree 0
    log_y = np.log10(y)
    
    try:
        # Fit power law: log(y) = α*log(x) + β
        slope, intercept, r_value, p_value, std_err = stats.linregress(log_x, log_y)
        power_loss = 1.0 - abs(r_value)  # 0 = perfect fit, 1 = no fit
        
        if slope > 0:
            power_loss = 1.0
            
        return power_loss
        
    except:
        return 1.0



def get_metrics(graph):
    metrics = {}
    
    # Degree
    degrees = [d for _, d in graph.degree()]
    metrics['average_degree'] = np.mean(degrees) if degrees else 0
    metrics['total_degree'] = sum(degrees)
    
    # Centrality
    if graph.number_of_nodes() > 1:
        try:
            # sum(max_degree - degree_i) / ((n-1)*(n-2))
            max_deg = max(degrees)
            n = graph.number_of_nodes()
            if n > 2:
                centralization = sum(max_deg - d for d in degrees) / ((n-1)*(n-2))
                metrics['degree_centralization'] = centralization
            else:
                metrics['degree_centralization'] = 0
        except:
            metrics['degree_centralization'] = 0
    else:
        metrics['degree_centralization'] = 0
        
    # Modality with louvain
    try:
        communities = nx.community.louvain_communities(graph)
        modularity = nx.community.modularity(graph, communities)
        metrics['network_modality'] = modularity
    except:
        metrics['network_modality'] = 0

    # Extra Metrics just to see
    metrics['node_count'] = graph.number_of_nodes()
    metrics['edge_count'] = graph.number_of_edges()
    metrics['density'] = nx.density(graph)
    
    # Largest component ratio
    components = list(nx.connected_components(graph))
    largest_component = max(components, key=len) if components else set()
    metrics['largest_component_ratio'] = len(largest_component) / graph.number_of_nodes() if graph.number_of_nodes() > 0 else 0
    metrics['num_connected_components'] = len(components)

    # Power loss scaling
    metrics['power_loss_scaling'] = calculate_power_loss_scaling(graph)


    return metrics



if __name__ == "__main__":
    csv_output = ""
    kinase_waters = get_kinase_waters()
    out_metrics = []
    for kinase in kinase_waters:
        if len(kinase['waters']) == 0:
            continue
        graph = nx.Graph()
        for water in kinase['waters']:
            graph.add_node(water['serial'],
                x=water['x'], y=water['y'], z=water['z']
                )
            for neighbor in kinase['waters']:
                distance = get_distance(water, neighbor)
                if distance <= 3.5 and distance >= 2.4 and water['serial'] != neighbor['serial']:
                    graph.add_edge(water['serial'], neighbor['serial'], weight=1)
        out_metrics.append({
            'PDB': kinase['PDB'],
            'metrics': get_metrics(graph)
        })
    with open('/dors/wankowicz_lab/ellas/network_metrics.json', 'w') as f:
        json.dump(out_metrics, f, indent=4)
