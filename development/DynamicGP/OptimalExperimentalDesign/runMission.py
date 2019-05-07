import numpy as np
from matplotlib import pyplot as plt
from scipy.stats import truncnorm
from missionControl import MissionControl

def get_truncated_normal(mean=0, sd=1, low=0, upp=10):
    return truncnorm((low - mean) / sd, (upp - mean) / sd, loc=mean, scale=sd)

missionGen = get_truncated_normal(mean = 1.1, sd = 0.001, low = 1, upp = 6)
mission = missionGen.rvs(10000)
missionGen = get_truncated_normal(mean = 2, sd = 1, low = 1, upp = 6)
mission = np.concatenate([mission,missionGen.rvs(1000)])
missionGen = get_truncated_normal(mean = 1.1, sd = 0.001, low = 1, upp = 6)
mission = np.concatenate([mission,missionGen.rvs(10000)])

damage = np.linspace(0, 1, 1000)
E = 1e9 - 0.5*1e9*damage**(0.9 - damage)
plt.plot(damage, E/1e9, color = 'black')
plt.plot([0, 1], [0.49, 0.49], color = 'black', linestyle = '--')
plt.xlabel('Damage Parameter')
plt.ylabel('$E/E_0$')
plt.xlim([0, 1])
plt.ylim([0, 1.25])
plt.savefig('Modulus degradation.pdf', figsize = [7, 7/1.618])


MC = MissionControl()

MC.demo_damage_probability_2()

mission = 3.5*np.ones(20000)


states, strains, flutter_speeds = MC.run(mission)
plt.plot(states, label = 'Damage State', color = 'red')
plt.scatter(range(len(mission)), mission, label = 'Load Factor', color = 'black', s = 0.01)
plt.plot(flutter_speeds/300, label = 'Normalized Flutter Speed', color = 'blue')
plt.ylim([0,6])
plt.xlabel('Time (seconds)')
plt.ylabel('Flutter speed (mph), Load Factor, (dimless) \n Damage State (dimless)')
plt.legend()
plt.savefig('constantManeuver_mission.pdf', figsize = [8, 8/1.628])
plt.close()

plt.scatter(range(len(mission)), 1000*strains, label = 'First sensor', s = 0.01)
plt.legend()
plt.xlabel('Time (seconds)')
plt.ylabel('Sensor Readings (1000 * strain)')
plt.savefig('noiselessStrains.pdf', figsize = [8, 8/1.628])
plt.close()

plt.scatter(100*strains, flutter_speeds)
plt.ylabel('Flutter Speed')
plt.xlabel('First Sensor')
plt.savefig('noiselessStrains_vs_flutter.pdf', figsize = [8, 8/1.628])
plt.close()


plt.plot(np.linspace(0,1, 1000), MC.flutter_speed_1d(np.linspace(0,1,1000)), color = 'black')
plt.ylabel('Flutter Speed')
plt.xlabel('Damage')
plt.xlim([0,1])
plt.ylim([0,330])
plt.savefig('damage_vs_flutter.pdf', figsize = [7, 7/1.618])
plt.close()

plt.plot(np.linspace(0,1, 1000), 1000*MC.strainFunc(np.linspace(0,1,1000), [24], 1), color = 'black', label = 'Load Factor = 1.0')
plt.plot(np.linspace(0,1, 1000), 1000*MC.strainFunc(np.linspace(0,1,1000), [24], 3.5), color = 'green', label = 'Load Factor = 3.5')
plt.plot(np.linspace(0,1, 1000), 1000*MC.strainFunc(np.linspace(0,1,1000), [24], 6), color = 'blue', label = 'Load Factor = 6.0')
plt.ylabel('Strain (x 1000)')
plt.xlabel('Damage')
plt.legend()
plt.xlim([0,1])
plt.savefig('damage_vs_strain.pdf', figsize = [7, 7/1.618])
plt.close()

strains = MC.runNoisy(states, mission)



plt.scatter(range(len(mission)), 100*strains, label = 'First sensor', s = 0.01)
plt.legend()
plt.xlabel('Time')
plt.ylabel('Sensor Readings (100 * strain)')
plt.savefig('noisyStrains.pdf', figsize = [8, 8/1.628])
plt.close()

plt.scatter(100*strains, flutter_speeds)
plt.ylabel('Flutter Speed')
plt.xlabel('Sensor Readings (1000 * strain)')
plt.savefig('noisyStrains_vs_flutter.pdf', figsize = [8, 8/1.628])
plt.close()