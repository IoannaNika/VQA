import torch.nn as nn
import torch

class LSTM(nn.Module):
    def __init__(self, input_size: int, hidden_size: int):
        super().__init__()

        self._lstm = nn.LSTM(input_size, hidden_size, bidirectional=False, batch_first=True)

    def forward(self, x_t: torch.Tensor):
        x_t = x_t.transpose(1,2)
        _, (h_t, _) = self._lstm(x_t)
        h_t = h_t.squeeze(0)
        # print("Hidden dimension", h_t.shape)
        # concatenate the forward and backward hidden states
        # h_t = torch.cat((h_t[0], h_t[1]), dim=1)
        # print("Hidden dimension", h_t.shape)
        return h_t