import numpy as np
from matplotlib import pyplot as plt
from scipy.stats import truncnorm
from missionControl import MissionControl

def get_truncated_normal(mean=0, sd=1, low=0, upp=10):
    return truncnorm((low - mean) / sd, (upp - mean) / sd, loc=mean, scale=sd)

constantLoad = 1.1
maxLoad = 4.0
numMCiters = 1000

missionGen = get_truncated_normal(mean = constantLoad, sd = 0.0001, low = 1, upp = 6)
mission = missionGen.rvs(10000)
for idx in range(250):
	# increase load linearly
	missionGen = get_truncated_normal(mean = constantLoad + idx*((maxLoad - constantLoad)/250.), sd = 0.01, low = 1, upp = 6)
	mission = np.concatenate([mission,missionGen.rvs(2)])
for idx in range(250):
	# increase load linearly
	missionGen = get_truncated_normal(mean = maxLoad, sd = 0.01, low = 1, upp = 6)
	mission = np.concatenate([mission,missionGen.rvs(2)])
for idx in range(250):
	missionGen = get_truncated_normal(mean = maxLoad + idx*((constantLoad - maxLoad)/250.), sd = 0.01, low = 1, upp = 6)
	mission = np.concatenate([mission,missionGen.rvs(2)])


missionGen = get_truncated_normal(mean = constantLoad, sd = 0.0001, low = 1, upp = 6)
mission = np.concatenate([mission,missionGen.rvs(10000)])

states = np.zeros([numMCiters, len(mission)])

strains = np.zeros([numMCiters, len(mission)])

for dmx in range(numMCiters):
	print(dmx)
	MC = MissionControl()

	states[dmx,:], strainsDiscard, flutter_speeds_discard = MC.run(mission)

	strains[dmx,:] = np.transpose(MC.runNoisy(states[dmx,:], mission))

np.savetxt('strains_%d'%maxLoad, strains)
np.savetxt('states_%d'%maxLoad, states)
flatStates = states.flatten('C')
plt.hist(flatStates, normed = True)
plt.xlabel('Damage Parameter')
plt.ylabel('Probability Density')
plt.show()

stateVar = np.var(states, axis = 0)
strainVar = np.var(strains, axis = 0)
states = np.mean(states, axis = 0)
strains = np.mean(strains, axis = 0)

fig, axes = plt.subplots(nrows = 3, ncols = 1)
axes[0].set_xlim([0, len(mission)])
axes[0].set_ylim([0, 1])

axes[1].set_xlim([0, len(mission)])
axes[1].set_ylim([min(100*strains)*0.95, max(100*strains)*1.1])

axes[2].set_xlim([0, len(mission)])
axes[2].set_ylim([0, 6])

#axes[0].scatter(range(len(mission)), states, label = 'Damage State', s = 0.01)
axes[0].errorbar(range(len(mission)), states, yerr = 2*stateVar)
axes[0].set_xlabel('Time')
axes[0].set_ylabel('Damage State')



#axes[1].scatter(range(len(mission)), 100*strains, s = 0.01)

axes[1].errorbar(range(len(mission)), 100*strains, yerr = 2*strainVar)
axes[1].set_xlabel('Time')
axes[1].set_ylabel('Strain')



axes[2].scatter(range(len(mission)), mission, s = 0.01)
axes[2].set_xlabel('Time')
axes[2].set_ylabel('Load Factor')



plt.savefig('noisyStrains.png', figsize = [8, 8/1.628])
plt.close()
