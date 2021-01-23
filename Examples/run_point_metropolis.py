import ctypes
import time
import math
import numpy as np
from numpy.ctypeslib import ndpointer
import matplotlib.pyplot as plt
from ctypes_classes import CCC, LC_irs

def random_normal_factor(sigma):
    return  np.random.normal(1.0, sigma, None)

class Obs:
    def __init__(self, mag, err):
         self.mag = mag
         self.err = err

class Param:
    def __init__(self, name, value, vary = False, min_value = 0.0, max_value = 1.0):
        self.name = name
        self.value = value
        self.vary = vary
        self.max_value = max_value
        self.min_value = min_value


# simplified metropolis algorithm
class Metropolis:
    def __init__(self, obs, errs, model):
        self.obs = obs
        self.errs = errs
        self.model = model
        self.chi_sequence = []
        self.a_sequence = []
        self.m2_sequence = []
        self.a = model["a"].value
        self.b = model["b"].value
        self.theta = model["theta"].value
        self.m2 = model["m2"].value
        self.m3 = model["m3"].value
        self.pos_ini_x = model["pos_ini_x"].value
        self.pos_ini_y = model["pos_ini_y"].value
        self.pos_fin_x = model["pos_fin_x"].value
        self.pos_fin_y = model["pos_fin_y"].value

    def get_chi_sq(self, ests):
        chi = 0.0
        for ob, err, est in zip(self.obs, self.errs, ests):
            chi += ((ob - est)/err)**2
        return chi

    def iterate(self, iters):
        sigma = 0.2

        lc_list       = np.zeros(len(self.obs), np.double)
        lc_list_trial = np.zeros(len(self.obs), np.double)

        lc_obj = LC_irs(self.a,
                        self.b,
                        self.theta,
                        self.m2,
                        self.m3,
                        source_size,
                        len(self.obs),
                        5)
        lc_obj.get_lc(pos_ini_x, pos_ini_y, pos_fin_x, pos_fin_y)
        lc_obj.copy_lc(lc_list)

        chi_sq = self.get_chi_sq(lc_list)
        number_of_accepted = 0

        for ite in range(0, iters):
            a_trial = self.a
            m2_trial = self.m2
            if model["a"].vary:
                a_sigma = sigma*abs(model["a"].max_value - model["a"].min_value)
                a_trial += np.random.normal(0.0, a_sigma, None)
                if a_trial > model["a"].max_value or a_trial < model["a"].min_value:
                    continue
            if model["m2"].vary:
                m2_sigma = sigma * abs(model["m2"].max_value - model["m2"].min_value)
                m2_trial += np.random.normal(0.0, m2_sigma, None)
                if m2_trial > model["m2"].max_value or m2_trial < model["m2"].min_value:
                    continue
            lc_obj = LC_irs(a_trial, b, theta, m2_trial, m3, source_size, len(self.obs), 5)
            lc_obj.get_lc(pos_ini_x, pos_ini_y, pos_fin_x, pos_fin_y)
            lc_obj.copy_lc(lc_list_trial)
            chi_sq_trial = self.get_chi_sq(lc_list_trial)
            if chi_sq_trial < chi_sq:
                # accept the jump
                self.a = a_trial
                self.m2 = m2_trial
                # ad-hoc of tuning of sigma
                sigma *= max(math.exp((chi_sq_trial-chi_sq)), 0.85)
                chi_sq = chi_sq_trial
                lc_list = lc_list_trial
                self.chi_sequence.append(chi_sq)
                self.a_sequence.append(self.a)
                self.m2_sequence.append(self.m2)
                number_of_accepted += 1
            else:
                # conditionally accept the jump with probability equal to ratio
                # trial likelyhood over current-step likelyhood
                if np.random.uniform() < math.exp(chi_sq-chi_sq_trial):
                    self.a = a_trial
                    self.m2 = m2_trial
                    chi_sq = chi_sq_trial
                    lc_list = lc_list_trial
                    self.chi_sequence.append(chi_sq)
                    self.a_sequence.append(self.a)
                    self.m2_sequence.append(self.m2)
                    number_of_accepted += 1


# Define lens parameters.
a = 1.0
b = 1.0
theta = 1.047197551
m2 = 1/5
m3 = 0.0
length = 500

start_time = time.time()
# Initialise the CriticalCurveCaustic object.
ccc = CCC(a,b,theta,m2,m3,length)
ccc.get_ca()

# Copy the data-points from CCC object to numpy arrays
cc_array = np.zeros(3000, np.cdouble)
ca_array = np.zeros(3000, np.cdouble)
ccc.copy_cc_ca(cc_array, ca_array)

cc_real = []
cc_imag = []
ca_real = []
ca_imag = []

for cc in cc_array:
  cc_real.append(cc.real)
  cc_imag.append(cc.imag)

for ca in ca_array:
  ca_real.append(ca.real)
  ca_imag.append(ca.imag)

# copy lens positions
lens_pos = np.zeros(3, np.cdouble)
ccc.copy_lenses(lens_pos)
lenses_real = []
lenses_imag = []
for lens in lens_pos:
  lenses_real.append(lens.real)
  lenses_imag.append(lens.imag)

# copy the bounding box
cc_min = np.zeros(1, np.cdouble)
cc_max = np.zeros(1, np.cdouble)
ca_min = np.zeros(1, np.cdouble)
ca_max = np.zeros(1, np.cdouble)
ccc.get_bounding_box(cc_min, cc_max, ca_min, ca_max, 1.5)
print("Bounding box cc: ", cc_min, cc_max)
print("Bounding box ca: ", ca_min, ca_max)

# imgs
pos_ini_x = -0.5
pos_ini_y = -0.51
pos_fin_x = 1.0
pos_fin_y = 0.577

# number of steps
lc_steps = 100
source_size = 1e-5
points_per_radius = 300

lc_point_array = np.zeros(lc_steps, np.double)
lc_irs_array = np.zeros(lc_steps, np.double)

lc_irs = LC_irs(a, b, theta, m2, m3, source_size, lc_steps, points_per_radius)

lc_irs.get_lc(pos_ini_x, pos_ini_y, pos_fin_x, pos_fin_y)
lc_irs.copy_lc(lc_point_array)

print("Copied LC")

# Chi squared will be determined while making random changes to modelled curve
chi_sq = 0.0
lc_obs_array = []
lc_err_array = []
for mag in lc_point_array:
    err = 0.1 * np.random.uniform(0.0,1.0)
    mag_obs = random_normal_factor(err) * mag
    lc_obs_array.append(mag_obs)
    lc_err_array.append(err*mag_obs)
    chi_sq += (mag-mag_obs)**2/(err*mag_obs)**2


# Plotting original model
#fig, (ax1, (ax2, ax3)) = plt.subplots(3, 2)

fig = plt.figure()

ax1 = plt.subplot(122)
ax2 = plt.subplot(221)
ax3 = plt.subplot(223)

ax1.set_title("Source Trajectory")

ax1.axis(xmin=cc_min.real, xmax=cc_max.real, ymin=cc_min.imag, ymax=cc_max.imag)
ax1.scatter(ca_real, ca_imag, s = 0.1)
ax1.scatter(lenses_real, lenses_imag, s=200.0, marker = 'o')
ax1.plot([pos_ini_x,pos_fin_x],[pos_ini_y,pos_fin_y])

ax2.set_title("Light Curve")

ax2.errorbar(x = np.linspace(0.0, 1.0, 100) , y = lc_obs_array, yerr = lc_err_array)

chi_label = r'$\chi ^2$'+"="+str("{:.2e}".format(chi_sq)) 

ax3.scatter(np.linspace(0.0, 1.0, 100), lc_obs_array-lc_point_array, label = chi_label)

legend = ax3.legend(loc='upper left', fontsize='x-small')

plt.tight_layout()
fig.savefig("LightCurveGenerated.png", dpi=200)

plt.close('all')

# Metropolis magic

model = {
            "a": Param("a", 0.7, True, 0.2, 1.7),
            "b": Param("b", 1.0),
            "theta": Param("theta", 1.047197551),
            "m2": Param("m2", 0.1, True, 0.001, 0.5),
            "m3": Param("m3", 0.0),
            "pos_ini_x": Param("pos_ini_x", pos_ini_x),
            "pos_ini_y": Param("pos_ini_y", pos_ini_y),
            "pos_fin_x": Param("pos_fin_x", pos_fin_x),
            "pos_fin_y": Param("pos_fin_y", pos_fin_y)
        }

metropolis_obj = Metropolis(lc_obs_array, lc_err_array, model)
metropolis_obj.iterate(1000)
print(metropolis_obj.a_sequence)
print(metropolis_obj.chi_sequence)
print("The chi_sq under metropolis equals = " + str(metropolis_obj.get_chi_sq(lc_point_array)))
print("The final value of a = " + str(metropolis_obj.a))
print("The final value of a = " + str(metropolis_obj.m2))
# Plotting a result of fit

ccc = CCC(metropolis_obj.a,b,theta,m2,m3,length)
ccc.get_ca()
ccc.copy_cc_ca(cc_array, ca_array)

cc_real = []
cc_imag = []
ca_real = []
ca_imag = []

for cc in cc_array:
  cc_real.append(cc.real)
  cc_imag.append(cc.imag)

for ca in ca_array:
  ca_real.append(ca.real)
  ca_imag.append(ca.imag)

# copy lens positions
ccc.copy_lenses(lens_pos)
lenses_real = []
lenses_imag = []
for lens in lens_pos:
  lenses_real.append(lens.real)
  lenses_imag.append(lens.imag)

# copy the bounding box
cc_min = np.zeros(1, np.cdouble)
cc_max = np.zeros(1, np.cdouble)
ca_min = np.zeros(1, np.cdouble)
ca_max = np.zeros(1, np.cdouble)
ccc.get_bounding_box(cc_min, cc_max, ca_min, ca_max, 1.5)

lc_point_array = np.zeros(lc_steps, np.double)
lc_irs_array   = np.zeros(lc_steps, np.double)

lc_irs = LC_irs(a, b, theta, m2, m3, source_size, lc_steps, points_per_radius)

lc_irs.get_lc(pos_ini_x,pos_ini_y,pos_fin_x,pos_fin_y)
lc_irs.copy_lc(lc_point_array)

# Plotting fitted model

fig = plt.figure()

ax1 = plt.subplot(122)
ax2 = plt.subplot(321)
ax3 = plt.subplot(323)
ax4 = plt.subplot(325)

ax1.set_title("Source Trajectory")

ax1.axis(xmin=cc_min.real, xmax=cc_max.real, ymin=cc_min.imag, ymax=cc_max.imag)
ax1.scatter(ca_real, ca_imag, s = 0.1)
ax1.scatter(lenses_real, lenses_imag, s=200.0, marker = 'o')
ax1.plot([pos_ini_x,pos_fin_x],[pos_ini_y, pos_fin_y])

ax2.set_title("Observed Light Curve")
ax2.errorbar(x = np.linspace(0.0, 1.0, 100), y=lc_obs_array, yerr=lc_err_array)
ax2.plot(np.linspace(0.0, 1.0, 100), lc_point_array)


chi_label = r'$\chi ^2$'+"="+str("{:.2e}".format(chi_sq))
ax3.scatter(np.linspace(0.0, 1.0, 100), lc_obs_array-lc_point_array, label=chi_label)
legend = ax3.legend(loc='upper left', fontsize='x-small')

ax4.set_title(r'$\chi ^2$')
ax4.scatter(np.linspace(1, len(metropolis_obj.chi_sequence), num=len(metropolis_obj.chi_sequence)), metropolis_obj.chi_sequence)

plt.tight_layout()
fig.savefig("LightCurveFitted.png", dpi=200)

plt.close('')
