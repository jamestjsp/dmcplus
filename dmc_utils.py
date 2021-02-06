import numpy as np
from scipy import signal
from mdl import get_dmc_model
import matplotlib.pyplot as plt


def interpolate_curve(SteadyStateTime, NumberOfCoefficients, curve, SampleInterval=1):
    xp = np.linspace(1, SteadyStateTime, NumberOfCoefficients)
    upsample = np.arange(1, SteadyStateTime+1, SampleInterval)
    return np.interp(upsample, xp, curve)


def rotate(v, dGain, nGain):
    '''
    This operation rotates the entire response curve around (0,0) to achieve a specified steady-state gain. 
    The operation is performed as follows. The program calculates the difference between the current 
    steady-state gain and the desired value. It then constructs a line that begins at zero and ends 
    at the calculated difference. Finally, this line is added to the response curve to achieve the rotation.

    No special steps are taken to force the final slope to zero. This operation can be used with either ramp 
    or steady-state variables.
    '''
    if isinstance(v, np.ndarray) == False:
        try:
            v = np.array(v)
        except:
            print(f'canot type cast {type(v)} to a numpy array')

    refrence_line = np.linspace(0, nGain - dGain, len(v))
    return v + refrence_line


def gScale(v, dGain, nGain):
    '''
    This operation will scale the coefficients so that the curve will end up with the user-specified 
    steady-state gain at the next to last coefficient. Each coefficient is multiplied by the ratio 
    of the user-specified gain / existing steady-state gain.
    This cannot be applied on Ramp Variable.
    '''
    if isinstance(v, np.ndarray) == False:
        try:
            v = np.array(v)
        except:
            print(f'canot type cast {type(v)} to a numpy array')
    gRatio = nGain / dGain
    return v * gRatio

def get_impulse_responce(SteadyStateTime, NumberOfCoefficients, curve):
    step_responce = interpolate_curve(SteadyStateTime, NumberOfCoefficients, curve)
    impulse_responce = np.diff(step_responce, n=1)
    t = t = np.arange(0, len(impulse_responce))
    return t, impulse_responce

def get_freq_responce(curve):
    _, impulse_responce = get_impulse_responce(SteadyStateTime, NumberOfCoefficients, curve)
    w, h = signal.freqz(impulse_responce)
    return w, np.abs(h)
with open('mdl/Stabi.mdl', 'r') as f:
    model = get_dmc_model(f)
    SteadyStateTime = model['SteadyStateTime']
    NumberOfCoefficients = model['NumberOfCoefficients']
    curve = model['Coefficients']['CV1-iPen-TOP']['MV2-RFX-SP']
    # 'CV2-nBut-BOT': {'MV1-TEMP-SP': -0.08109767701483682}
    v = model['Coefficients']['CV2-nBut-BOT']['MV1-TEMP-SP']
    dGain = model['dGain']['CV2-nBut-BOT']['MV1-TEMP-SP']
    rv = rotate(v, dGain, -0.095)
    scale_v = gScale(v, dGain, -0.095)
    x = np.arange(0, len(v))
    plt.plot(x, v, '-', x, rv, '--', x, scale_v, '-.',)
    plt.legend(['Orginal = -0.0810', 'Rotate = -0.095', 'GScale  = -0.095'])
    plt.grid(color='r', linestyle='--', linewidth=0.5)
    plt.show()
    # plt.plot(t, impulse_responce, '--')
    # plt.show()
    # plt.semilogx(w, np.abs(h), 'r')
    # plt.ylabel('Amplitude (db)', color='b')
    # plt.xlabel('Frequency (rad/sample)', color='b')
    # plt.show()
