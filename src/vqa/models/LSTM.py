import torch.nn as nn
import torch

class LSTM(nn.Module):
    def __init__(self, input_size: int, hidden_size: int):
        super().__init__()

        self._lstm = nn.LSTM(input_size, hidden_size, bidirectional=False, batch_first=True)

    def forward(self, x_t: torch.Tensor):
        _, (h_t, _) = self._lstm(x_t)
        return h_t