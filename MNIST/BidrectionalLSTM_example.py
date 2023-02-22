#RuntimeError: For unbatched 2-D input, hx should also be 2-D but got 3-D tensor
#Imports
import torch
import torchvision
import torch.nn as nn
import torch.optim as optim
import torch.nn.functional as F
from torch.utils.data import DataLoader
import torchvision.datasets as datasets
import torchvision.transforms as transforms

device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

#MNIST is Nx1x28x28 N=64, or 28 time sequences, each with 28 features
#Hyperparameters
input_size = 28
hidden_size = 256
num_layers = 2
num_classes = 10
sequence_length = 28    #one row at a time for each time step
learning_rate = 0.001
batch_size = 64
num_epochs = 2

#Create a bidirectional LSTM
class BRNN(nn.Module):
    def __init__(self, input_size, hidden_size, num_layers, num_classes):
        super(BRNN, self).__init__()
        self.hidden_size = hidden_size
        self.num_layers = num_layers
        self.lstm = nn.LSTM(input_size, hidden_size, num_layers, batch_first=True, bidirectional = True)
        self.fc = nn.Linear(hidden_size*2, num_classes) #1 layer forward, 1 backward, both get concatenated

    def forward(self, x):
        h0 = torch.zeros(self.num_layers*2, x.size(0), self.hidden_size).to(device) #1 layer forward, 1 backward, both get concatenated 
        c0 = torch.zeros(self.num_layers*2, x.size(0), self.hidden_size).to(device)

        out, _ = self.lstm(x, (h0,c0)) #output here is (hidden state, cell state)
        out = self.fc(out[:,-1,:])

        return out
    
#Load data
train_dataset = datasets.MNIST(root = 'dataset/', train = True, transform = transforms.ToTensor(), download = True)
test_dataset = datasets.MNIST(root = 'dataset/', train = False, transform = transforms.ToTensor(), download = True)

train_loader = DataLoader(dataset = train_dataset, batch_size=batch_size, shuffle = True)
test_loader = DataLoader(dataset = test_dataset, batch_size=batch_size, shuffle = True)

#Initialize Network
model = BRNN(input_size, hidden_size, num_layers, num_classes).to(device)

#Loss and Optimizer
criterion = nn.CrossEntropyLoss()
optimizer = optim.Adam(model.parameters(), lr = learning_rate)

#Train Network
for epoch in range(num_epochs): #an epoch means the network has seen all the images in a dataset
    for batch_idx, (data, targets) in enumerate(train_loader): #go through each batch in training loader, see what batch index we have. data = images and correct digit per image
        #Get data to cuda if possible
        data = data.to(device=device).squeeze(1) #will remove the 1 from the Nx1x28x28
        targets = targets.to(device=device)

        #forward
        scores = model(data)
        loss = criterion(scores, targets)

        #backwards
        optimizer.zero_grad()
        loss.backward()

        #gradient descent or adam step
        optimizer.step()

#Check accuracy on training and test model fit
def check_accuracy(loader, model):
    if loader.dataset.train:
        print('Checking accuracy on training data')
    else:
        print("Checking accuracy on test data")
    num_correct=0
    num_samples=0
    model.eval()

    with torch.no_grad():
        for x,y in loader: 
            x= x.to(device= device).squeeze(1)
            y = y.to(device=device)

            scores = model(x) #64x10 
            _, predictions = scores.max(1)
            num_correct += (predictions == y).sum()
            num_samples += predictions.size(0)
        
        print(f'Got {num_correct} / {num_samples} with accuracy {float(num_correct)/float(num_samples)*100:.2f}')

    model.train()
    #return accuracy
check_accuracy(train_loader, model)
check_accuracy(test_loader, model)
