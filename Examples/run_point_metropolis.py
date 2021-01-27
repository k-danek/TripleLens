import ctypes
import time
import math
import copy
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
        self.trial_model = copy.deepcopy(model)  # TODO: think of better way how to represent trial data!
        self.chi_sequence = []
        self.chi_trial_sequence = []
        self.chi_re_sequence = []
        self.a = model["a"].value
        self.b = model["b"].value
        self.theta = model["theta"].value
        self.m2 = model["m2"].value
        self.m3 = model["m3"].value
        self.pos_ini_x = model["pos_ini_x"].value
        self.pos_ini_y = model["pos_ini_y"].value
        self.pos_fin_x = model["pos_fin_x"].value
        self.pos_fin_y = model["pos_fin_y"].value
        self.lc_list = []
        self.chi_sq = 0.0
        self.a_sequence = []
        self.m2_sequence = []

    def get_chi_sq(self, ests):
        chi = 0.0
        for ob, err, est in zip(self.obs, self.errs, ests):
            chi += ((ob - est)/err)**2
        return chi

    def get_trial(self, model_name, sigma):
        if model_name in self.trial_model:
            trial_sigma = sigma * abs(self.trial_model[model_name].max_value - self.trial_model[model_name].min_value)
            trial_value = self.model[model_name].value + np.random.normal(0.0, trial_sigma, None)
            if trial_value > self.trial_model[model_name].max_value or\
               trial_value < self.trial_model[model_name].min_value:
                return copy.copy(self.model[model_name].value)
            else:
                return trial_value
        else:
            print("WARNING: did not find variable: "+model_name)

    def iterate(self, iters):
        sigma = 0.2

        self.lc_list = np.zeros(len(self.obs), np.double)
        lc_list_trial = np.zeros(len(self.obs), np.double)

        lc_obj_init = LC_irs(self.a,
                             self.b,
                             self.theta,
                             self.m2,
                             self.m3,
                             source_size,
                             len(self.obs),
                             5)

        lc_obj_init.get_lc(self.model["pos_ini_x"].value,
                           self.model["pos_ini_y"].value,
                           self.model["pos_fin_x"].value,
                           self.model["pos_fin_y"].value)

        lc_obj_init.copy_lc(self.lc_list)

        self.chi_sq = self.get_chi_sq(self.lc_list)
        number_of_accepted = 0

        # picks a list of keys that have varied parameters
        varied_params = [key for key in model if model[key].vary]

        for ite in range(0, iters):
            for param_name in varied_params:
                self.trial_model[param_name].value = self.get_trial(param_name, sigma)
                # checking the type before inserting into ctypes functions
                if not isinstance(self.trial_model[param_name].value, float):
                    print("WARNING: wrong type for a model parameter: " + param_name + " is "
                          + str(type(self.trial_model[param_name].value)))

            lc_obj_trial = LC_irs(self.trial_model["a"].value,
                                  self.trial_model["b"].value,
                                  self.trial_model["theta"].value,
                                  self.trial_model["m2"].value,
                                  self.trial_model["m3"].value,
                                  source_size,
                                  len(self.obs),
                                  5)

            lc_obj_trial.get_lc(self.trial_model["pos_ini_x"].value,
                                self.trial_model["pos_ini_y"].value,
                                self.trial_model["pos_fin_x"].value,
                                self.trial_model["pos_fin_y"].value)

            lc_obj_trial.copy_lc(lc_list_trial)
            chi_sq_trial = self.get_chi_sq(lc_list_trial)
            if chi_sq_trial < self.chi_sq:
                # accept the jump
                for param_name in varied_params:
                    self.model[param_name].value = copy.copy(self.trial_model[param_name].value)

                # ad-hoc of tuning of sigma
                sigma *= max(math.exp((chi_sq_trial-self.chi_sq)), 0.95)
                self.chi_sq = copy.copy(chi_sq_trial)
                self.chi_sequence.append(self.chi_sq)
                self.lc_list = copy.copy(lc_list_trial)
                number_of_accepted += 1

            else:
                # conditionally accept the jump with probability equal to ratio
                # trial likelihood over current-step likelihood
                # for most of the use this never happens as differences in chi_sq tend to be huge
                if np.random.uniform() < math.exp(self.chi_sq-chi_sq_trial):
                    # accept the jump
                    for param_name in varied_params:
                        self.model[param_name].value = copy.copy(self.trial_model[param_name].value)

                    self.chi_sq = copy.copy(chi_sq_trial)
                    self.lc_list = copy.copy(lc_list_trial)
                    self.chi_sequence.append(self.chi_sq)
                    number_of_accepted += 1
                    self.chi_re_sequence.append(self.get_chi_sq(self.lc_list))
                    self.chi_trial_sequence.append(self.get_chi_sq(lc_list_trial))

        for param in varied_params:
            print("Varied params: " + param)
        for key in self.model:
            print(str(key)+"="+str(self.model[key].value))

        print("final chi sq = "+str(self.chi_sq))
        print(self.chi_sequence)


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
            "a": Param("a", 0.9, True, 0.8, 1.1),
            "b": Param("b", 1.0),
            "theta": Param("theta", 1.047197551),
            "m2": Param("m2", 0.22, True, 0.15, 0.25),
            "m3": Param("m3", 0.0),
            "pos_ini_x": Param("pos_ini_x", -0.45, True, -0.6, -0.2),
            "pos_ini_y": Param("pos_ini_y", -0.4, True, -0.6, -0.3),
            "pos_fin_x": Param("pos_fin_x", 0.7, True, 0.6, 1.05),
            "pos_fin_y": Param("pos_fin_y", 0.4, True, 0.3, 0.7)
        }

metropolis_obj = Metropolis(lc_obs_array, lc_err_array, model)
metropolis_obj.iterate(100000)
print(metropolis_obj.chi_sequence)

a_result = metropolis_obj.model["a"].value
b_result = metropolis_obj.model["b"].value
theta_result = metropolis_obj.model["theta"].value
m2_result = metropolis_obj.model["m2"].value
m3_result = metropolis_obj.model["m3"].value
lc_result = metropolis_obj.lc_list

# Plotting a result of fit

ccc = CCC(a_result,
          b_result,
          theta_result,
          m2_result,
          m3_result,
          length)
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

# Plotting fitted model

fig = plt.figure()

ax1 = plt.subplot(122)
ax2 = plt.subplot(321)
ax3 = plt.subplot(323)
ax4 = plt.subplot(325)

ax1.set_title("Source Trajectory")

ax1.axis(xmin=cc_min.real, xmax=cc_max.real, ymin=cc_min.imag, ymax=cc_max.imag)
ax1.scatter(ca_real, ca_imag, s=0.1)
ax1.scatter(lenses_real, lenses_imag, s=100.0, marker='o')
ax1.plot([pos_ini_x, pos_fin_x], [pos_ini_y, pos_fin_y], label="true", color='#eeeeee')
ax1.plot([metropolis_obj.model["pos_ini_x"].value,
          metropolis_obj.model["pos_fin_x"].value],
         [metropolis_obj.model["pos_ini_y"].value,
          metropolis_obj.model["pos_fin_y"].value],
         label="fitted positions")


ax2.set_title("Observed Light Curve")
ax2.errorbar(x=np.linspace(0.0, 1.0, 100), y=lc_obs_array, yerr=lc_err_array)
ax2.plot(np.linspace(0.0, 1.0, 100), metropolis_obj.lc_list)

chi_label = r'$\chi ^2$'+"="+str("{:.2e}".format(metropolis_obj.chi_sq))
ax3.scatter(np.linspace(0.0, 1.0, 100), lc_obs_array-lc_point_array, label=chi_label)
legend = ax3.legend(loc='upper left', fontsize='x-small')

ax4.set_title(r'$\chi ^2$')
ax4.scatter(np.linspace(1, len(metropolis_obj.chi_sequence),
            num=len(metropolis_obj.chi_sequence)),
            metropolis_obj.chi_sequence)

plt.tight_layout()
fig.savefig("LightCurveFitted.png", dpi=200)

plt.close('')
