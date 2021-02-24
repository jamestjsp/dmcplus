import numpy as np
from scipy import signal
from scipy.optimize import curve_fit
from robustcontrol.utils import InternalDelay, tf
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
    _, impulse_rsponse = get_impulse_rsponse(
        SteadyStateTime, NumberOfCoefficients, curve)
    w, h = signal.freqz(impulse_rsponse)
    return w, np.abs(h)


def calculate_deadtime(curve):
    '''
    This function calulate dead time from the FIR curve.
    It is recomented to pass interpolated curve as cuves can a compressed one. 
    '''
    theta = 0
    for coef in curve:
        if coef == 0:
            theta += 1
        else:
            break
    return theta


def detect_steady_state(curve, deadtime, SteadyStateTime, threshold=0.0005):
    window_size = int(SteadyStateTime/40) if int(SteadyStateTime/40) > 10 else 10
    start = 0 if deadtime == 0 else deadtime
    end = start + window_size
    while True:
        stdev = np.std(curve[start:end])
        start = end
        end = start + window_size
        if stdev <= threshold:
            break
        elif end > len(curve):
            end = len(curve)
            break
        else:
            continue
    return end


def tf_step(ts, c, theta, K, a, b, isRamp=False):
    s = tf([1, 0])
    tf_model = InternalDelay(
        K * (a*s + 1) / (b*s**2 + c*s + (not isRamp)) * np.exp(-theta*s))

    def uf(t): return np.array([1])
    return tf_model.simulate(uf, ts).flatten()


def fit_tf(curve, K, isRamp, SteadyStateTime, NumberOfCoefficients, fopdt=False):
    extended_curve = interpolate_curve(
        SteadyStateTime, NumberOfCoefficients, curve)
    if np.count_nonzero(extended_curve):
        theta = calculate_deadtime(extended_curve)
        SteadyStateTime = detect_steady_state(
            extended_curve, theta, SteadyStateTime, 0.0001)
        ts = np.arange(0, SteadyStateTime)
        if not fopdt and not isRamp:
            popt, pcov = curve_fit(lambda ts, a, b, c: tf_step(ts, c, theta, K, a, b, isRamp),
                                ts,
                                extended_curve[:SteadyStateTime],
                                bounds=(
                [-100, 0, 0], [100, ((SteadyStateTime-theta)/4)**2, ((SteadyStateTime-theta)/4)])
            )
            a, b, c = popt
            num = [a*K, K]
            den = [b, c, 1]
            deadtime = theta
            return [num, den, deadtime]
        elif fopdt and not isRamp:
            popt, pcov = curve_fit(lambda ts, c, theta: tf_step(ts, c, theta, K, 0, 0, isRamp),
                                ts,
                                extended_curve[:SteadyStateTime],
                                p0 = [1,theta],
                                bounds=([0, theta], [SteadyStateTime/4, theta+SteadyStateTime/4])
            )
            c, theta = popt
            num = [K]
            den = [c, 1]
            deadtime = theta
            return [num, den, deadtime]
        else:
            popt, pcov = curve_fit(lambda ts, theta: tf_step(ts, 1, theta, K, 0, 0, isRamp),
                                ts,
                                extended_curve,
                                p0 = [theta],
                                bounds=(theta, theta+(SteadyStateTime-theta)/4)
            )
            theta = popt
            num = [K]
            den = [1, 0]
            deadtime = theta
            return [num, den, deadtime]
    else:
        return None


def fir2tf(model, fopdt=False):
    SteadyStateTime = model['SteadyStateTime']
    NumberOfCoefficients = model['NumberOfCoefficients']
    deps = model['Dependents']
    inds = model['Independents']
    isRamp = model['isRamp']
    dGain = model['dGain']
    curves = model['Coefficients']
    tf_dict = dict()
    for dep in deps:
        tf_dict[dep] = dict()
        for ind in inds:
            curve = curves[dep][ind]
            K = dGain[dep][ind]
            ramp_flag = isRamp[dep]
            tf = fit_tf(curve, K, ramp_flag, SteadyStateTime,
                        NumberOfCoefficients, fopdt)
            tf_dict[dep][ind] = tf
    return tf_dict


def compare_curve(dep, ind,fir_model, tf_model, SteadyStateTime, NumberOfCoefficients):
    fir_curve = fir_model[dep][ind]
    tf_curve = tf_model[dep][ind]
    ts = np.arange(0, SteadyStateTime)
    fir_step_response = interpolate_curve(
        SteadyStateTime, NumberOfCoefficients, fir_curve)
    tf_model = InternalDelay(
        tf(tf_curve[0], tf_curve[1], deadtime=tf_curve[2]))

    uf = lambda t: np.array([1])
    tf_step_response = tf_model.simulate(uf, ts).flatten()
    plt.title(dep + ' vs ' + ind)
    plt.plot(ts, fir_step_response, '-', color='b', label='dmc_fir_curve')
    plt.plot(ts, tf_step_response, '-', color='m', label='tf_step_response')
    plt.grid(color='r', linestyle='--', linewidth=0.5)
    plt.legend()
    plt.show()

if __name__ == '__main__':

    with open('mdl/upstram_HP.mdl', 'r') as f:
        model = get_dmc_model(f)
    tf_model = fir2tf(model, fopdt=True)
    SteadyStateTime = model['SteadyStateTime']
    NumberOfCoefficients = model['NumberOfCoefficients']

    dep = 'FIC-2102'
    ind = 'FI-2005'
    fir_model = model['Coefficients']
    compare_curve(dep, ind, fir_model, tf_model, SteadyStateTime, NumberOfCoefficients)
    dGain = model['dGain'][dep][ind]
    v = model['Coefficients'][dep][ind]
    ng = -6  # new gain
    rv = rotate(v, dGain, ng)
    scale_v = gScale(v, dGain, ng)
    x = np.arange(0, len(v))
    t, impulse_rsponse = get_impulse_rsponse(
        SteadyStateTime, NumberOfCoefficients, v)
    w, h = get_freq_rsponse(SteadyStateTime, NumberOfCoefficients, v)
    # Plot Gain correction
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
