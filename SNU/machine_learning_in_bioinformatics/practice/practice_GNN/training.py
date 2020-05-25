import torch
import torch.nn.functional as F
from torch_geometric.data import DataLoader
from torch_geometric.nn import GCNConv
from torch_geometric.nn import global_max_pool as gmp
from graph_construction import BRCA_PPI_Test_Dataset, BRCA_PPI_Train_Dataset

from sklearn.metrics import roc_auc_score
use_gpu = torch.cuda.is_available()
print(use_gpu)


def label_distribution(data_set):
    y_1 = 0
    y_0 = 0
    for d in data_set:
        if d.y == torch.Tensor([1]):
            y_1 += 1
        else:
            y_0 += 1
    print('#y_1: ', str(y_1), ';#y_0: ', str(y_0))
    return y_1 + y_0


def get_train_test(batch_size):
    train_set = BRCA_PPI_Train_Dataset(root='BRCA_PPI_Train_Dataset')
    test_set = BRCA_PPI_Test_Dataset(root='BRCA_PPI_Test_Dataset')

    train_loader = DataLoader(train_set, batch_size=batch_size, shuffle=True)
    test_loader = DataLoader(test_set, batch_size=batch_size, shuffle=False)

    return train_loader, test_loader


class PracticeNet(torch.nn.Module):
    def __init__(self):
        super(PracticeNet, self).__init__()
        self.conv1 = GCNConv(1, 32)
        self.conv2 = GCNConv(32, 32)

        self.lin2 = torch.nn.Linear(32, 2)

    def forward(self, data):
        x, edge_index, batch = data.x, data.edge_index, data.batch
        x = F.relu(self.conv1(x, edge_index))
        x = F.relu(self.conv2(x, edge_index))
        x = gmp(x, batch)
        x = self.lin2(x)
        x = F.softmax(x, dim=1)

        return x


def one_hot(label, batch_size, class_num):
    """
    scalar label to one hot encoding with batch.
    """
    label = label.resize_(batch_size, 1)
    m_zeros = torch.zeros(batch_size, class_num).to(device)
    onehot = m_zeros.scatter_(1, label, 1)  # (dim,index,value)

    return onehot.to(torch.long)


def train(train_loader):
    model.train()
    train_length = 0
    total_loss_train = 0
    for index, data in enumerate(train_loader):
        train_length += data.num_graphs
        data = data.to(device)
        optimizer.zero_grad()
        out = model(data)
        y_onehot = one_hot(data.y, batch_size=data.y.shape[0], class_num=2)
        loss = LOSS(out, y_onehot.to(torch.float))

        loss.backward()
        total_loss_train += loss.item() * data.num_graphs
        optimizer.step()

    train_loss = total_loss_train / train_length
    return train_loss


def test(loader):
    model.eval()
    pred_y = []
    data_y = []
    for data in loader:
        data = data.to(device)
        out = model(data)
        data_y.append(data.y)
        pred_y.append(out[:, 1])
    y_true = torch.cat(data_y, dim=0).cpu().detach().numpy()
    y_pred = torch.cat(pred_y, dim=0).cpu().detach().numpy()
    auc = roc_auc_score(y_true, y_pred)
    return auc


def train_main():
    train_loader, test_loader = get_train_test(batch_size=16)
    for epoch in range(20):
        train_loss = train(train_loader)
        train_results = test(train_loader)
        test_results = test(test_loader)
        print(('Epoch: {:03d}, train_loss: {:.4f}, Train_auc: {:.3f}, Test_auc: {:.3f}').format(
            epoch, train_loss, train_results, test_results))
        # print(train_results[0], test_results[0])
        # print(train_results[1], test_results[1])
        # print(train_results[2], test_results[2])
        #


if __name__ == '__main__':

    device = torch.device('cuda:2' if torch.cuda.is_available() else 'cpu')
    model = PracticeNet().to(device)
    optimizer = torch.optim.Adam(model.parameters(), lr=0.001)
    LOSS = torch.nn.BCELoss()
    train_main()

