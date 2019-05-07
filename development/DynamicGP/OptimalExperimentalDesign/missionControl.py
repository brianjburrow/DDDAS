import numpy as np
from matplotlib import pyplot as plt
from scipy.stats import truncnorm

class MissionControl():
    def __init__(self):
        # edit parameters here
        self.gauge_locs = [25]
        pass

    def demo_damage_probability(self):
        import matplotlib.colors as colors
        import matplotlib.cm as cmx
        damages = np.linspace(0, 1, 1000)
        load_factor = np.linspace(1, 6, 5000)
        values = np.zeros([len(load_factor), len(damages)])
        jet = cm = plt.get_cmap('jet') 
        cNorm  = colors.Normalize(vmin=1, vmax=6)
        scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)

        Z = [[0,0],[0,0]]
        Min, Max = (1, 6)
        levels = np.linspace(Min,Max,15)
        CS3 = plt.contourf(Z, levels, cmap='jet')
        plt.clf()

        for idx, loads in enumerate(load_factor):
            values[idx] = self.state_evolution(damages, loads)
        for idx,vector in enumerate(values):
            plt.plot(vector, color = scalarMap.to_rgba(load_factor[idx]))


        plt.ylim([0,1.1])
        plt.xlim([0,1.1])
        plt.ylabel('Probability of Damage Propagation')
        plt.xlabel('Current Damage State')
        plt.colorbar(CS3)
        plt.savefig('damage_prob_demo.pdf', figsize = [8, 8/1.628])
        pass
    def demo_damage_probability_2(self):
        load_factor = np.linspace(1, 6, 5000)
        damage = 0.5
        values = np.zeros_like(load_factor)
        fig, ax = plt.subplots(1)

        for idx, loads in enumerate(load_factor):
            blank, values[idx] = self.state_evolution(0.0, loads, demo2 = True)
        ax.plot(load_factor, values, color = 'blue', label = 'Previous State: 0.0')

        for idx, loads in enumerate(load_factor):
            blank, values[idx] = self.state_evolution(damage, loads, demo2 = True)
        ax.plot(load_factor, values, color = 'black', label = 'Previous State: 0.5')

        for idx, loads in enumerate(load_factor):
            blank, values[idx] = self.state_evolution(0.99, loads, demo2 = True)
        ax.plot(load_factor, values, color = 'green', label = 'Previous State: 0.99')

        ax.set_ylabel('Probability of Damage Propagation')
        ax.set_xlabel('Load Factor')
        ax.legend()
        plt.savefig('damage_prob_demo_2.pdf', figsize = [7, 7/1.618])
        pass

    def get_truncated_normal(self, mean=0, sd=1, low=0, upp=10):
        return truncnorm((low - mean) / sd, (upp - mean) / sd, loc=mean, scale=sd)

    def state_evolution(self,current_state, load_factor, demo2 = False):
        b = 0.99975/5 # no longer used
        a = 27./5.
        c = -127./5.
        value = (((load_factor/10)**2))*(current_state - 0.5)**2 / (1000.0) + 5*(self.sigmoid(a*(load_factor) + c - 2*current_state))*b
        u = np.random.uniform(0,1,1)
        try:
            len(current_state)
            return value
        except:
            if u < value:
                damage_sampler = self.get_truncated_normal(
                    mean = 0.001, 
                    sd = 0.1, 
                    low = 0, 
                    upp = (1 - current_state)
                    )
                current_state += damage_sampler.rvs(1)
        if demo2:
            return current_state, value
        else:
            return current_state
        #return damage + damage_event

    def flutter_speed_1d(self,damage):
        return 320 - 300*self.sigmoid(-7 + 14*damage)

    def strainFunc(self,damage, gauge_locs, load_factor):
        E = 1e9 - 0.5*1e9*damage**(0.9 - damage)
        M = np.zeros_like(gauge_locs)
        try:
            for idx, val in enumerate(gauge_locs):
                M[idx] = load_factor * 1000 * val**2 / 2.
        except:
            M = load_factor * 1000 * gauge_locs**2 / 2.
        strains = M/E
        return strains

    def strainFuncNoisy(self,damage, gauge_locs, load_factor):
        E = 1e9 - 0.2*1e9*damage**2
        M = np.zeros_like(gauge_locs)
        for idx, val in enumerate(gauge_locs):
            M[idx] = load_factor * 1000 * val**2 / 2
        strains = np.zeros_like(gauge_locs)
        for idx, val in enumerate(M):
            strains = val/E + np.random.normal(0, 0.000005, size = 1)
        return strains

    def sigmoid(self, x):
      return 1 / (1 + np.exp(-x))

    def run(self, mission):
        states = np.zeros(len(mission))
        strains = np.zeros([len(mission), len(self.gauge_locs)])
        for idx, maneuver in enumerate(mission):
            strains[idx] = self.strainFunc(states[idx], self.gauge_locs, maneuver)
            if idx == len(mission) - 1:
                break
            if states[idx] == 1:
                states[idx + 1] = 1
            else:
                states[idx + 1] = self.state_evolution(states[idx], maneuver)
        flutter_speeds = self.flutter_speed_1d(states)
        return states, strains[:,0], flutter_speeds

    def runNoisy(self, states, mission):
        strains = np.zeros([len(mission), len(self.gauge_locs)])
        for idx, maneuver in enumerate(mission):
            strains[idx] = self.strainFuncNoisy(states[idx], self.gauge_locs, maneuver)
            if idx == len(mission) - 1:
                break
        return strains

