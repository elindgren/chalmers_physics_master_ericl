import pandas as pd
import networkx as nx
import numpy as np
from scipy import sparse
import os
import torch
from torch_geometric.data import Data, DataLoader
from torch_geometric.data import InMemoryDataset, Dataset


def sparse_mx_to_torch_sparse_tensor(sparse_mx):
    """Convert a scipy sparse matrix to a torch sparse tensor."""
    sparse_mx = sparse_mx.tocoo().astype(np.float32)
    indices = torch.from_numpy(
        np.vstack((sparse_mx.row, sparse_mx.col)).astype(np.long))

    return indices


def processing_topology(ppi_path):
    ppi_edges = pd.read_csv(ppi_path, sep='\t', header=0)
    ppi_graph = nx.from_pandas_edgelist(ppi_edges, 'A', 'B')
    ppi_nodes = ppi_graph.nodes()
    ppi_adj = np.array(nx.adjacency_matrix(ppi_graph).todense())
    ppi_edge_index = sparse_mx_to_torch_sparse_tensor(sparse.csr_matrix(ppi_adj).tocoo())
    print(ppi_edge_index.shape)
    return ppi_nodes, ppi_edge_index


def construct_PPIDataset(patient_path, patient_label_path, ppi_nodes, ppi_edges):
    pt_list = []
    patient_gene_expression = pd.read_csv(patient_path, sep='\t', header=0)
    patient_label = pd.read_csv(patient_label_path, sep='\t', header=0)

    print(patient_gene_expression)
    patient_name = patient_gene_expression.columns.tolist()[1:]
    for pt in patient_name:
        if pt in patient_label['sample_id'].values:
            pt_x = []
            pt_exp = patient_gene_expression.loc[:, ['gene', pt]]
            # print(pt_exp)
            for node in ppi_nodes:
                node_exp = pt_exp[pt_exp['gene'] == node][pt].values
                pt_x.append(node_exp)
            pt_x = torch.Tensor(pt_x).to(torch.float)
            pt_edge_index = ppi_edges
            pt_y = patient_label[patient_label['sample_id'] == pt]['sample_type'].values
            print(pt_y, pt)
            if pt_y[0] == 'Primary Tumor':
                pt_y = torch.Tensor(np.array([1])).to(torch.long)
            elif pt_y[0] == 'Solid Tissue Normal':
                pt_y = torch.Tensor(np.array([0])).to(torch.long)

            pt_graph = Data(x=pt_x, y=pt_y, edge_index=pt_edge_index)
            pt_list.append(pt_graph)
        else:
            continue

    return pt_list


class BRCA_PPI_Train_Dataset(InMemoryDataset):
    def __init__(self, root, transform=None, pre_transform=None):
        super(BRCA_PPI_Train_Dataset, self).__init__(root, transform, pre_transform)
        self.data, self.slices = torch.load(self.processed_paths[0])

    @property
    def processed_file_names(self):
        return ['BRCA_PPI_Train_Dataset']

    def process(self):
        data, slices = self.collate(train_pt_list)
        torch.save((data, slices), self.processed_paths[0])


class BRCA_PPI_Test_Dataset(InMemoryDataset):
    def __init__(self, root, transform=None, pre_transform=None):
        super(BRCA_PPI_Test_Dataset, self).__init__(root, transform, pre_transform)
        self.data, self.slices = torch.load(self.processed_paths[0])

    @property
    def processed_file_names(self):
        return ['BRCA_PPI_Test_Dataset']

    def process(self):
        data, slices = self.collate(test_pt_list)
        torch.save((data, slices), self.processed_paths[0])


if __name__ == '__main__':
    PPI_path = 'selected_PPI_network'
    train_patient_path = 'train_pt_gene_expression'
    test_patient_path = 'test_pt_gene_expression'

    patient_label_path = 'pt_label'
    ppi_nodes, ppi_edges = processing_topology(PPI_path)
    train_pt_list = construct_PPIDataset(train_patient_path, patient_label_path, ppi_nodes, ppi_edges)
    test_pt_list = construct_PPIDataset(test_patient_path, patient_label_path, ppi_nodes, ppi_edges)

    BRCA_PPI_Train_Dataset(root='BRCA_PPI_Train_Dataset')
    BRCA_PPI_Test_Dataset(root='BRCA_PPI_Test_Dataset')

