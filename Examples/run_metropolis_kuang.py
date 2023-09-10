import ctypes
from ctypes import *
import time
import math
import copy
import numpy as np
from numpy.ctypeslib import ndpointer
import matplotlib.pyplot as plt
from ctypes_classes import CCC, LC_cuda


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
        self.a = copy.deepcopy(model["a"].value)
        self.b = copy.deepcopy(model["b"].value)
        self.theta = copy.deepcopy(model["theta"].value)
        self.m2 = copy.deepcopy(model["m2"].value)
        self.m3 = copy.deepcopy(model["m3"].value)
        self.pos_ini_x = copy.deepcopy(model["pos_ini_x"].value)
        self.pos_ini_y = copy.deepcopy(model["pos_ini_y"].value)
        self.pos_fin_x = copy.deepcopy(model["pos_fin_x"].value)
        self.pos_fin_y = copy.deepcopy(model["pos_fin_y"].value)
        self.limb_dark = copy.deepcopy(model["limb_dark"].value)
        self.source_size = copy.deepcopy(model["source_size"].value)
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
        sigma = 0.03
        points_per_radius = 30
        stringBuffer = create_string_buffer(b"linear")

        self.lc_list = np.zeros(len(self.obs), np.double)
        lc_list_trial = np.zeros(len(self.obs), np.double)

        lc_obj_init = LC_cuda(copy.deepcopy(self.a),
                              copy.deepcopy(self.b),
                              copy.deepcopy(self.theta),
                              copy.deepcopy(self.m2),
                              copy.deepcopy(self.m3),
                              copy.deepcopy(self.source_size),
                              len(self.obs),
                              points_per_radius)

        lc_obj_init.set_limb_darkening(stringBuffer,copy.deepcopy(self.model["limb_dark"].value))

        lc_obj_init.get_lc_cuda(copy.deepcopy(self.model["pos_ini_x"].value),
                                copy.deepcopy(self.model["pos_ini_y"].value),
                                copy.deepcopy(self.model["pos_fin_x"].value),
                                copy.deepcopy(self.model["pos_fin_y"].value))

        lc_obj_init.copy_lc(self.lc_list)

        self.chi_sq = self.get_chi_sq(self.lc_list)
        number_of_accepted = 0

        # picks a list of keys that have varied parameters
        varied_params = [key for key in model if model[key].vary]

        for ite in range(0, iters):
            print("Iteration started: " + str(ite))
            for param_name in varied_params:
                self.trial_model[param_name].value = self.get_trial(param_name, sigma)
                # checking the type before inserting into ctypes functions
                if not isinstance(self.trial_model[param_name].value, float):
                    print("WARNING: wrong type for a model parameter: " + param_name + " is "
                          + str(type(self.trial_model[param_name].value)))


            print("Iteration initialization about to start: " + str(ite))

            #lc_obj_trial = LC_cuda(self.trial_model["a"].value,
            #                       self.trial_model["b"].value,
            #                       self.trial_model["theta"].value,
            #                       self.trial_model["m2"].value,
            #                       self.trial_model["m3"].value,
            #                       source_size,
            #                       len(self.obs),
            #                       points_per_radius)

            print("Iteration initialized: " + str(ite))
            lc_obj_init.set_limb_darkening(stringBuffer,copy.deepcopy(self.model["limb_dark"].value))

            print("Iteration ld set: " + str(ite))
            lc_obj_init.get_lc_cuda(copy.deepcopy(self.trial_model["pos_ini_x"].value),
                                    copy.deepcopy(self.trial_model["pos_ini_y"].value),
                                    copy.deepcopy(self.trial_model["pos_fin_x"].value),
                                    copy.deepcopy(self.trial_model["pos_fin_y"].value))

            print("Iteration calculation finished: " + str(ite) + " chi^2=" + str(self.chi_sq))
            lc_obj_init.copy_lc(lc_list_trial)
            #print("Iteration calculation copied: " + str(ite))

            chi_sq_trial = self.get_chi_sq(lc_list_trial)
            self.chi_trial_sequence.append(chi_sq_trial)

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



# Chi squared will be determined while making random changes to modelled curve
chi_sq = 0.0
lc_obs_array =[15.965657	,
16.020279	,
16.076065	,
16.133014	,
16.191125	,
16.2504	,
16.310838	,
16.372442	,
16.435213	,
16.499154	,
16.564268	,
16.630559	,
16.698032	,
16.76669	,
16.836539	,
16.907586	,
16.979837	,
17.053298	,
17.127977	,
17.203882	,
17.281023	,
17.359406	,
17.439044	,
17.519945	,
17.60212	,
17.68558	,
17.770339	,
17.856407	,
17.943798	,
18.032525	,
18.122604	,
18.214048	,
18.306873	,
18.401096	,
18.496734	,
18.593804	,
18.692325	,
18.792316	,
18.893797	,
18.996789	,
19.101315	,
19.207396	,
19.315058	,
19.424324	,
19.535222	,
19.647779	,
19.762022	,
19.877983	,
19.995693	,
20.115185	,
20.236494	,
20.359656	,
20.484711	,
20.611699	,
20.740664	,
20.871651	,
21.004708	,
21.139889	,
21.277247	,
21.416841	,
21.558735	,
21.702995	,
21.849694	,
21.99891	,
22.150726	,
22.305235	,
22.462536	,
22.622735	,
22.785951	,
22.952313	,
23.121963	,
23.295058	,
23.471773	,
23.652299	,
23.836855	,
24.025683	,
24.219057	,
24.417291	,
24.62074	,
24.829815	,
25.044993	,
25.266831	,
25.495987	,
25.73325	,
25.979572	,
26.236129	,
26.504389	,
28.089206	,
30.185362	,
32.25443	,
34.300544	,
36.328585	,
38.351384	,
40.403571	,
42.577927	,
45.442518	,
47.404358	,
48.559131	,
49.146223	,
49.187147	,
48.548492	,
46.543083	,
45.067104	,
43.85201	,
42.708655	,
41.580693	,
40.443312	,
39.287864	,
38.10926	,
36.904512	,
36.247871	,
36.783625	,
37.348009	,
37.941152	,
38.563584	,
39.216244	,
39.900472	,
40.618014	,
41.371069	,
42.162359	,
42.995222	,
43.873739	,
44.802911	,
45.788907	,
46.839383	,
47.963963	,
49.174903	,
50.488085	,
51.924517	,
53.512788	,
55.29316	,
57.325069	,
59.702241	,
62.588298	,
66.324975	,
72.041504	,
80.382127	,
83.336028	,
84.094025	,
83.4666	,
81.612392	,
78.92873	,
75.544687	,
71.413643	,
66.70181	,
61.459996	,
55.735862	,
49.597716	,
43.091449	,
36.264875	,
29.135985	,
21.740411	,
16.595914	,
16.60562	,
16.606903	,
16.599913	,
16.584838	,
16.5619	,
16.53135	,
16.493468	,
16.448557	,
16.396937	,
16.338947	,
16.274937	,
16.205263	,
16.13029	,
16.050381	,
15.9659	,
15.877208	,
15.784658	,
15.688595	,
15.589357	,
15.487266	,
15.382637	,
15.275768	,
15.166943	,
15.056434	,
14.944495	,
14.831368	,
14.717279	,
14.602439	,
14.487044	,
14.371279	,
14.255311	,
14.139297	,
14.023382	,
13.907696	,
13.79236	,
13.677484	,
13.563167	,
13.449499	,
13.336561	,
13.224425	,
13.113157	,
13.002814	,
12.893447	,
12.7851	,
12.677811	,
12.571615]
lc_err_array = [0.01] * 199


# kkuang parameters

#a =	1.396
#b =	1.168
#theta =	5.332
#m1 =	0.968738799
#m2 =	0.02809342517
#m3 =	0.003167775873
#length = 500
#pos_ini_x = 0.1162461156
#pos_ini_y = -0.006323071317
#pos_fin_x = -0.02947981448
#pos_fin_y = -0.03489234325
#source_size = 0.0044



# Metropolis magic
model = {
            "a": Param("a", 1.396, False, 1.394, 1.398),
            "b": Param("b", 1.168, False, 1.166, 1.170),
            "theta": Param("theta", 5.332),
            "m2": Param("m2", 0.02809342517, False, 0.15, 0.25),
            "m3": Param("m3", 0.003167775873),
            "pos_ini_x": Param("pos_ini_x", 0.11618131837256242, True, 0.10, 0.13),
            "pos_ini_y": Param("pos_ini_y", -0.005533691430735394, True, -0.01, 0.00001),
            "pos_fin_x": Param("pos_fin_x", -0.0303097378671771, True, -0.039, -0.025),
            "pos_fin_y": Param("pos_fin_y", -0.03562362921880297, True, -0.045, -0.032),
            "limb_dark": Param("limb_dark", 0.24123490581531616, True, 0.0, 0.5),
            "source_size": Param("source_size", 0.005, False, 0.0, 0.7)
        }

metropolis_obj = Metropolis(lc_obs_array, lc_err_array, model)

#print(metropolis_obj.chi_sequence)
metropolis_obj.iterate(500)
#print(metropolis_obj.chi_sequence)

a_result = metropolis_obj.model["a"].value
b_result = metropolis_obj.model["b"].value
theta_result = metropolis_obj.model["theta"].value
m2_result = metropolis_obj.model["m2"].value
m3_result = metropolis_obj.model["m3"].value
x1_result = metropolis_obj.model["pos_ini_x"].value
y1_result = metropolis_obj.model["pos_ini_y"].value
x2_result = metropolis_obj.model["pos_fin_x"].value
y2_result = metropolis_obj.model["pos_fin_y"].value
ld_result = metropolis_obj.model["limb_dark"].value
lc_result = metropolis_obj.lc_list

length = 500
cc_array = np.zeros(3000, np.cdouble)
ca_array = np.zeros(3000, np.cdouble)

# Plotting a result of fit
ccc = CCC(copy.copy(a_result),
          copy.copy(b_result),
          copy.copy(theta_result),
          copy.copy(m2_result),
          copy.copy(m3_result),
          copy.copy(length)
        )
ccc.get_ca()
ccc.copy_cc_ca(cc_array, ca_array)
lens_pos = np.zeros(3, np.cdouble)
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


f = open("fit.dat", "w")
for amp in metropolis_obj.lc_list:
    f.write(str(amp)+"\n")
f.close()

f = open("chis_trial.dat", "w")
for chi in metropolis_obj.chi_trial_sequence:
    f.write(str(chi)+"\n")
f.close()

f = open("chis.dat", "w")
for chi in metropolis_obj.chi_sequence:
    f.write(str(chi)+"\n")
f.close()

f = open("reldiff.dat", "w")
for index in range(0,len(lc_obs_array)):
    diff = metropolis_obj.lc_list[index]-lc_obs_array[index]
    reldiff = (metropolis_obj.lc_list[index]-lc_obs_array[index])/lc_obs_array[index]
    f.write(str(lc_obs_array[index])+" "+str(diff)+" "+str(reldiff)+"\n")
f.close()

# Plotting fitted model

fig = plt.figure()

ax1 = plt.subplot(122)
ax2 = plt.subplot(321)
ax3 = plt.subplot(323)
ax4 = plt.subplot(325)

ax1.set_title("Source Trajectory")

ax1.axis(xmin=cc_min.real, xmax=cc_max.real, ymin=cc_min.imag, ymax=cc_max.imag)
ax1.scatter(ca_real, ca_imag, s=0.1)
ax1.scatter(lenses_real, lenses_imag, s=10.0, marker='o')
ax1.plot([x1_result, x2_result], [y1_result, y2_result], label="true", color='#eeeeee')
ax1.plot([metropolis_obj.model["pos_ini_x"].value,
          metropolis_obj.model["pos_fin_x"].value],
         [metropolis_obj.model["pos_ini_y"].value,
          metropolis_obj.model["pos_fin_y"].value],
          label="fitted positions", color='blue')


ax2.set_title("Observed Light Curve")
#ax2.errorbar(x=np.linspace(0.0, 1.0, len(lc_obs_array)), y=lc_obs_array, yerr=lc_err_array)
ax2.plot(np.linspace(0.0, 1.0, len(lc_obs_array)), lc_obs_array, color="red")
ax2.plot(np.linspace(0.0, 1.0, len(metropolis_obj.lc_list)), metropolis_obj.lc_list, color="blue")

chi_label = r'$\chi ^2$'+"="+str("{:.2e}".format(metropolis_obj.chi_sq))
ax3.scatter(np.linspace(0.0, 1.0, len(lc_obs_array)), (lc_obs_array-metropolis_obj.lc_list)/lc_obs_array, label=chi_label)
legend = ax3.legend(loc='upper left', fontsize='x-small')

ax4.set_title(r'$\chi ^2$')
ax4.scatter(np.linspace(0.0, 1.0, len(metropolis_obj.chi_sequence)),
            metropolis_obj.chi_sequence, s=20.0, marker="+")

plt.tight_layout()
fig.savefig("LightCurveFittedKuang.png", dpi=200)
fig.savefig("LightCurveFittedKuang.eps", format="eps")

plt.close('')

print("printed everything, sleeping now")
time.sleep(3)

print("woke up")
