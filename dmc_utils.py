import numpy as np
from scipy import signal
from mdl import get_dmc_model
import matplotlib.pyplot as plt


def interpolate_curve(SteadyStateTime, NumberOfCoefficients, curve, SampleInterval=1):
    '''
    This function will interpolate a curve of it is compressed.
    '''
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


def get_impulse_rsponse(SteadyStateTime, NumberOfCoefficients, curve):
    '''
    This function generates the impulse rsponse of step rsponse curve.
    '''
    step_rsponse = interpolate_curve(
        SteadyStateTime, NumberOfCoefficients, curve)
    impulse_rsponse = np.diff(step_rsponse, n=1)
    t = t = np.arange(0, len(impulse_rsponse))
    return t, impulse_rsponse


def get_freq_rsponse(SteadyStateTime, NumberOfCoefficients, curve):
    '''
    This function generates the frequency rsponse of step rsponse curve.
    '''
    _, impulse_rsponse = get_impulse_rsponse(SteadyStateTime, NumberOfCoefficients, curve)
    w, h = signal.freqz(impulse_rsponse)
    return w, np.abs(h)
def calculate_deadtime(curve):
    '''
    This function calulate dead time from the FIR curve.
    It is recomented to pass interpolated curve as cuves can a compressed one. 
    '''
    theta = 0
    for coef in curve:
        if coef ==0:
            theta += 1
        else:
            break
    return theta
if __name__ == '__main__':
    
    with open('mdl/Stabi.mdl', 'r') as f:
        model = get_dmc_model(f)
        
    SteadyStateTime = model['SteadyStateTime']
    NumberOfCoefficients = model['NumberOfCoefficients']
    v = model['Coefficients']['CV2-nBut-BOT']['MV1-TEMP-SP']
    dGain = model['dGain']['CV2-nBut-BOT']['MV1-TEMP-SP']
    ng = - 0.095 #new gain
    rv = rotate(v, dGain, ng)
    scale_v = gScale(v, dGain, ng)
    x = np.arange(0, len(v))
    t, impulse_rsponse = get_impulse_rsponse(SteadyStateTime, NumberOfCoefficients, v)
    w, h = get_freq_rsponse(SteadyStateTime, NumberOfCoefficients, v)
    #Plot Gain correction
    plt.plot(x, v, '-', x, rv, '--', x, scale_v, '-.',)
    plt.legend(['Orginal = -0.0810', f'Rotate = {ng}', f'GScale  = {ng}'])
    plt.grid(color='r', linestyle='--', linewidth=0.5)
    plt.title('Edited Gain using rotate and gScale')
    plt.show()
    # Plot Impulse rsponse
    plt.plot(t, impulse_rsponse, '--')
    plt.grid(color='r', linestyle='--', linewidth=0.5)
    plt.title('Impulse rsponse')
    plt.show()
    # Plot frequency reponce
    plt.semilogx(w, h, 'g')
    plt.ylabel('Amplitude (db)', color='b')
    plt.xlabel('Frequency (rad/sample)', color='b')
    plt.title('Frequency rsponse')
    plt.grid(color='r', linestyle='--', linewidth=0.5)
    plt.show()
