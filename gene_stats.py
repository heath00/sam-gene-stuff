import numpy as np
from scipy import stats

class gene:
	def __init__(self, wt_subgene, ko_subgene, gene_num):
		self.wt = wt_subgene
		self.ko = ko_subgene
		self.gene_num = gene_num

	def populate_subgenes(self, runs):
		# gets all files into subgenes
		for run in runs:
			self.wt.populate_array(run.create_slice(self.gene_num))
			self.ko.populate_array(run.create_slice(self.gene_num))
		# self.wt.populate_array(run1.create_slice(self.gene_num))
		# self.wt.populate_array(run2.create_slice(self.gene_num))
		# self.ko.populate_array(run1.create_slice(self.gene_num))
		# self.ko.populate_array(run2.create_slice(self.gene_num))

		# begin stats calculations
		self.wt.stats_setup()
		self.ko.stats_setup()

		self.wt_n = self.wt

		self.name = self.ko.get_name()

	def print_all(self):
		self.ko.print_arrays()
		self.wt.print_arrays()
		print("ddCT gapdh: ", self.ddct_gapdh)
		print("ddCT hprt: ", self.ddct_hprt)
		print("ddCT sd gapdh: ", self.sd_ddct_gapdh)
		print("ddCT sd hprt: ", self.sd_ddct_hprt)
		print("fc gapdh: ", self.fc_gapdh)
		print("fc hprt: ", self.fc_hprt)
		print("std fc gapdh: ", self.sd_fc_gapdh)
		print("std fc hprt: ", self.sd_fc_hprt)
		print("std err gapdh: ", self.std_err_gapdh)
		print("std err hprt: ", self.std_err_hprt)
		print("gapdh ttest: ", self.gapdh_t_test)
		print("hprt ttest: ", self.hprt_t_test)



	def ddct_calculations(self):
		# calculate ddct via ko_mean - wt_mean
		wt_gapdh_mean, wt_hprt_mean = self.wt.get_means()
		ko_gapdh_mean, ko_hprt_mean = self.ko.get_means()
		self.ddct_gapdh = ko_gapdh_mean - wt_gapdh_mean
		self.ddct_hprt = ko_hprt_mean - wt_hprt_mean

		# # calculate std of ddct via sqrt(ko_std^2 + wt_std^2)
		wt_gapdh_std, wt_hprt_std = self.wt.get_stds()
		ko_gapdh_std, ko_hprt_std = self.ko.get_stds()
		# self.sd_ddct_gapdh = (ko_gapdh_std**2 + wt_gapdh_std**2)**.5
		# self.sd_ddct_hprt = (ko_hprt_std**2 + wt_hprt_std**2)**.5

		# CHANGES
		# now calculate (ko_std^2/num_ko_after_ouliers_removed + wt_std^2/num_wt_after_ouliers_removed)
		# t statistic -> degrees of freedom = min(num_ko_after_outliers_removed, num_wt_after_outliers_removed) - 1
		# 12.71 -> 1.960
		# now just print upper bound fold change and lower bound fold change, stop printing the error

		# make t table
		t_table = {1:  12.71,\
				   2:  4.303,\
				   3:  3.182,\
				   4:  2.776,\
				   5:  2.571,\
				   6:  2.447,\
				   7:  2.365,\
				   8:  2.306,\
				   9:  2.262,\
				   10: 2.228,\
				   11: 2.201,\
				   12: 2.179,\
				   13: 2.160,\
				   14: 2.145,\
				   15: 2.131}
		
		gapdh_deg_freedom = min(self.wt.gapdh_len(), self.ko.gapdh_len()) - 1
		print('gapdh_deg_freedom = ' , gapdh_deg_freedom)
		self.sd_ddct_gapdh = t_table[gapdh_deg_freedom] * ( ((ko_gapdh_std**2) / self.ko.gapdh_len()) + ((wt_gapdh_std**2) / self.wt.gapdh_len()) )**.5 

		hprt_deg_freedom = min(self.wt.hprt_len(), self.ko.hprt_len()) - 1

		# this conditional is to avoid the hprt calculation for geomean. messy, but -1 just acts as a flag here since hprt_len is 0
		if hprt_deg_freedom != -1:
			print('hprt_deg_freedom = ' , hprt_deg_freedom)
			self.sd_ddct_hprt = t_table[hprt_deg_freedom] * ( ((ko_hprt_std**2) / self.ko.hprt_len()) + ((wt_hprt_std**2) / self.wt.hprt_len()) )**.5
		else:
			self.sd_ddct_hprt = 1

		# calculate uppers and lowers
		self.ddct_upper_gapdh = self.ddct_gapdh + self.sd_ddct_gapdh
		self.ddct_lower_gapdh = self.ddct_gapdh - self.sd_ddct_gapdh
		self.ddct_upper_hprt = self.ddct_hprt + self.sd_ddct_hprt
		self.ddct_lower_hprt = self.ddct_hprt - self.sd_ddct_hprt

		# self.debug_printer()

	def debug_printer(self):
		print('gene name: ', self.name)
		print('-------------------------------------------')
		print('uppder ddct gapdh', self.ddct_upper_gapdh)
		print('lower ddct gapdh', self.ddct_lower_gapdh)
		print('uppder ddct hprt', self.ddct_upper_hprt)
		print('lower ddct hprt', self.ddct_lower_hprt)

		print('wt_gapdh: ', self.wt.gapdh.get_array())
		print('ko_gapdh: ', self.ko.gapdh.get_array())
		print('wt_hprt: ', self.wt.hprt.get_array())
		print('ko_hprt: ', self.ko.hprt.get_array())

	def foldchange_helper(self, ddct_val):
		if ddct_val < 0:
			return 2**(-ddct_val)
		else:
			return -1/(2**-ddct_val)

	def foldchange_ul_helper(self, ddct_val, ul):
		if ddct_val < 0:
			return 2**(-ul)
		else:
			return -1/(2**-ul)

	def foldchange_calculations(self):
		# these calculations are confusing but here they are

		#original foldchange
		self.fc_gapdh = self.foldchange_helper(self.ddct_gapdh)
		self.fc_hprt = self.foldchange_helper(self.ddct_hprt)

		# foldchange upper/lower calcs
		self.fc_upper_gapdh = self.foldchange_ul_helper(self.ddct_gapdh, self.ddct_upper_gapdh)
		self.fc_lower_gapdh = self.foldchange_ul_helper(self.ddct_gapdh, self.ddct_lower_gapdh)
		self.fc_upper_hprt = self.foldchange_ul_helper(self.ddct_hprt, self.ddct_upper_hprt)
		self.fc_lower_hprt = self.foldchange_ul_helper(self.ddct_hprt, self.ddct_lower_hprt)

		# standard deviations
		self.sd_fc_gapdh = np.std([self.fc_upper_gapdh, self.fc_gapdh, self.fc_lower_gapdh], ddof=1)
		self.sd_fc_hprt = np.std([self.fc_upper_hprt, self.fc_hprt, self.fc_lower_hprt], ddof=1)

	def std_error(self, len_arr1, len_arr2, sd):
		return sd / (min(len_arr1, len_arr2)**.5)

	def std_error_calculations(self):
		self.std_err_gapdh = self.std_error(self.ko.gapdh_len(), self.wt.gapdh_len(), self.sd_fc_gapdh)
		self.std_err_hprt = self.std_error(self.ko.hprt_len(), self.wt.hprt_len(), self.sd_fc_hprt)

	def t_tests(self):
		self.gapdh_t_test = stats.ttest_ind(self.ko.get_gapdh(), self.wt.get_gapdh(), equal_var=False)[1]
		self.hprt_t_test = stats.ttest_ind(self.ko.get_hprt(), self.wt.get_hprt(), equal_var=False)[1]

	def get_post_outlier_numbers(self,stat_dict):
		stat_dict['gapdh_wt'] = self.wt.gapdh_len()
		stat_dict['gapdh_ko'] = self.ko.gapdh_len()
		stat_dict['hprt_wt'] = self.wt.hprt_len()
		stat_dict['hprt_ko'] = self.ko.hprt_len()


class subgene:
	def __init__(self, genotype):
		self.genotype = genotype
		# general data holders
		self.gene_name = ''
		self.gapdh = subarr([])
		self.hprt = subarr([])

	# Functions for array handling

	# Grabs all the numbers and puts them in the arrays
	def populate_array(self, gene_file_slice):
		self.gene_name = str(gene_file_slice[1][2])
		for row in gene_file_slice[1:]:	
			#print(row)
			if row[1] == self.genotype:
				if row[5]:
					self.gapdh.append(float(row[5]))
					self.hprt.append(float(row[12]))

	def get_name(self):
		return self.gene_name

	# Prints array contents
	def print_arrays(self, array_type=True):
		print(self.gene_name, self.genotype, 'gapdh: ', self.gapdh.get_contents())
		print(self.gene_name, self.genotype, 'hprt: ', self.hprt.get_contents())


	# Functions for statistics
	def stats_setup(self):
		self.gapdh.calc_preoutlier_stats()
		self.hprt.calc_preoutlier_stats()

		self.gapdh.remove_outliers()
		self.hprt.remove_outliers()

	def get_means(self):
		return self.gapdh.get_mean(), self.hprt.get_mean()

	def get_stds(self):
		return self.gapdh.get_std(), self.hprt.get_std()

	def gapdh_len(self):
		return self.gapdh.get_len()

	def hprt_len(self):
		return self.hprt.get_len()

	def get_gapdh(self):
		return self.gapdh.get_array()

	def get_hprt(self):
		return self.hprt.get_array()

class subarr:
	def __init__(self, arr):
		self.arr = arr

	def append(self, num):
		# python is weird...
		return self.arr.append(num)

	def print_contents(self):
		print(self.arr)

	def get_contents(self):
		return self.arr

	def calc_preoutlier_stats(self):
		# get q75, q25, iqr
		self.q75, self.q25 = np.percentile(self.arr, [75,25])
		self.iqr = self.q75 - self.q25

	def remove_outliers(self):
		iqr_mult = 1.5*self.iqr
		pre_removal_len = len(self.arr)
		self.arr = [i for i in self.arr if ((i < (iqr_mult + self.q75)) and (i > (self.q25 - iqr_mult)))]

	def get_mean(self):
		self.mean = np.mean(self.arr)
		return self.mean

	def get_std(self):
		self.std = np.std(self.arr, ddof=1)
		return self.std

	def get_len(self):
		return len(self.arr)

	def get_array(self):
		return self.arr