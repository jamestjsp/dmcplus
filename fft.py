import numpy as np
import pandas as pd
from scipy.fftpack import fft, ifft, fftfreq
from matplotlib import pyplot as plt
# Number of sample points

df = pd.read_csv('hp.csv', index_col='Time', parse_dates=True)
N = 120
start = 0
for stop in range(N, len(df), N):
    sig1 = df['WELL1.OP'].to_numpy()[start:stop]
    sig2 = df['KOPC_OP'].to_numpy()[start:stop]
    sig1_fft = fft(sig1)
    sig2_fft = fft(sig2)
    xf = fftfreq(N)[:N//2]
    amp_ratio = np.abs(sig2_fft[1:N//2]) / np.abs(sig1_fft[1:N//2])
    plt.plot(xf[1:N//2], amp_ratio)
    start = stop
    print(xf)
plt.show()
