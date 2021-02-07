import numpy as np
from scipy.optimize import curve_fit
from dmc_utils import interpolate_curve
from robustcontrol.utils import InternalDelay, tf
from matplotlib import pyplot as plt
SteadyStateTime = 360
NumberOfCoefficients = 120
dGain = 0.8535463
isRamp = False
curve = np.array([
    0.00002221823, 0.0009924865, 0.004621449, 0.01261975, 0.03178242, 0.05862877, 0.09139621, 0.1283221, 0.1641151,
    0.1963617, 0.2248391, 0.2493245, 0.268948, 0.286579, 0.3027854, 0.3181353, 0.3331322, 0.3480265, 0.362754, 0.3772505,
    0.3911198, 0.4044743, 0.4173378, 0.4297337, 0.4416337, 0.4533324, 0.4648598, 0.4762461, 0.487401, 0.498352, 0.5090767,
    0.5195526, 0.529628, 0.5394433, 0.5490087, 0.5583341, 0.5673309, 0.5761245, 0.5847198, 0.5931217, 0.6012423, 0.6091794,
    0.6169351, 0.6245115, 0.6318258, 0.6389706, 0.6459485, 0.6527622, 0.6593369, 0.6657568, 0.6720241, 0.6781407, 0.6840372,
    0.6897907, 0.6954033, 0.7008768, 0.7061481, 0.7112887, 0.7163008, 0.7211865, 0.7258887, 0.7304726, 0.7349401, 0.7392929,
    0.7434795, 0.7475588, 0.7515326, 0.7554028, 0.7591226, 0.7627456, 0.7662735, 0.7697079, 0.7730069, 0.7762188, 0.7793452,
    0.7823876, 0.7853084, 0.7881511, 0.7909172, 0.7936081, 0.7961902, 0.7987025, 0.8011464, 0.8035231, 0.8058028, 0.8080204,
    0.8101771, 0.8122739, 0.8142846, 0.81624, 0.8181413, 0.8199896, 0.8217614, 0.8234843, 0.8251593, 0.8267874, 0.8283478,
    0.8298651, 0.83134, 0.8327734, 0.8341473, 0.835483, 0.8367815, 0.8380434, 0.8392529, 0.8404289, 0.8415721, 0.8426832,
    0.8437482, 0.8447838, 0.8457907, 0.8467694, 0.8477077, 0.8486203, 0.8495076, 0.8503703, 0.8511976, 0.8520024, 0.8527851, 0.8535463
])
extended_curve = interpolate_curve(
    SteadyStateTime, NumberOfCoefficients, curve)
s = tf([1, 0])
K = dGain

def tf_step(ts, a, b, c, theta, K, isRamp=False):
    if isRamp:
        tf_model = InternalDelay(
            K * (a*s + 1) / (b*s**2 + c*s + 0) * np.exp(-theta*s))
    else:
        tf_model = InternalDelay(
            K * (a*s + 1) / (b*s**2 + c*s + 1) * np.exp(-theta*s))

    def uf(t): return np.array([1])
    return tf_model.simulate(uf, ts).flatten()

ts = np.arange(0, SteadyStateTime)
popt, pcov = curve_fit(lambda ts, a, b, c,theta:tf_step(ts, a, b, c, theta, K, isRamp=False),
                        ts,
                        extended_curve,
                        bounds=([-50,0, -50, 0],[50,50, 50, 20])
                        )
a, b, c, theta = popt
sim_out = tf_step(ts, a, b, c, theta, K, isRamp=False)
print(popt)
plt.plot(extended_curve)
plt.plot(sim_out)
plt.show()