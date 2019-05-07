import numpy as np
from matplotlib import pyplot as plt
from scipy.stats import truncnorm
from scipy.stats import norm
from missionControl import MissionControl as mc
####### Define Helper Functions
def get_truncated_normal(mean=0, sd=1, low=0, upp=10):
    return truncnorm((low - mean) / sd, (upp - mean) / sd, loc=mean, scale=sd)

def mix_gaussian_distributions(x, means, sigmas):
	# Input: 
	# 	x: location to evaluate mixture.  
	# 	means: vector of Gaussian means
	# 	sigmas: vector of Gaussian standard deviations
	# Return:
	# 	GMM probability density evaluated at ''x''.

	num_dists = len(means)
	value = 0
	for idx in range(num_dists):
		value += norm.pdf(x, loc = means[idx], scale = sigmas[idx])
	return value/num_dists

def convert_GMM_to_Gaussian(support, means, sigmas):
	# Create a gaussian mixture model from a set of gaussian random variables
	# which are represented by vectors of mean and standard deviations.
	# Then, the GMM is converted to a single Gaussian distribution
	# Input:
	#       support: discretized support of the GMM.  Vector of values
	#       means:   vector of Gaussian means
	#		sigmas:  vector of Gaussian standard deviations
	# Return: 
	#		mean, standard deviation of the single Gaussian dist.
	num_dists = len(means)
	value = 0
	for idx in range(num_dists):
		value += norm.pdf(support, loc = means[idx], scale = sigmas[idx])
	value /= num_dists
	mean = 0
	variance = 0
	for idx in range(num_dists):
		mean += means[idx]
	mean /= num_dists
	for idx in range(num_dists):
		variance += (sigmas[idx] + (means[idx] - mean)**2)
	variance /= num_dists
	return mean, variance**0.5


def pass_strains_to_matlab(damages, loadFactor, gauge_loc):
	mission = mc()
	strains = 1000*mission.strainFunc(damages, loadFactor, gauge_loc)
	return strains

def pass_flutter_to_matlab(damages, loadFactor, gauge_loc):
	mission = mc()
	flutter = mission.flutter_speed_1d(damages)
	return flutter