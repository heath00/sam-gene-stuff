import csv
import numpy as np
from scipy import stats


def not_outlier(n, iqr, q25, q75):
	iqr_mult = 1.5*iqr
	return ((n < (iqr_mult + q75)) and (n > (q25 - iqr_mult)))

def ddct_handler(ddct): 
	if ddct < 0:
		return 2**(-ddct)
	else:
		return -1/(2**-ddct)

def fc_ul_handler(ddct, ul):
	if ddct < 0:
		return 2**(-ul)
	else:
		return -1/(2**-ul)

def std_err(len_arr1, len_arr2, sd):
	return sd / (max(len_arr1, len_arr2)**.5)


fname = input("Input file 1: ")


ko_gapdh = []
ko_hprt = []
wt_gapdh = []
wt_hprt = []

with open(fname) as csv_file:
	csv_reader = csv.reader(csv_file, delimiter=',')
	linect = 0
	for row in csv_reader:
		if linect in range(18, 35):
			if row[1] == "KO":
				if row[5]:
					ko_gapdh.append(float(row[5]))
					ko_hprt.append(float(row[12]))
			elif row[1] == "WT":
				if row[5]:
					wt_gapdh.append(float(row[5]))
					wt_hprt.append(float(row[12]))
		linect += 1

fname2 = input("Input file 2: ")

with open(fname2) as csv_file:
	csv_reader = csv.reader(csv_file, delimiter=',')
	linect = 0
	for row in csv_reader:
		if linect in range(18, 35):
			if row[1] == "KO":
				if row[5]:
					ko_gapdh.append(float(row[5]))
					ko_hprt.append(float(row[12]))
			elif row[1] == "WT":
				if row[5]:
					wt_gapdh.append(float(row[5]))
					wt_hprt.append(float(row[12]))
		linect += 1

print("ko_gapdh: ", ko_gapdh)
print("ko_hprt: ", ko_hprt)
print("wt_gapdh: ", wt_gapdh)
print("wt_hprt: ", wt_hprt)

q75_ko_gapdh, q25_ko_gapdh = np.percentile(ko_gapdh, [75, 25])
iqr_ko_gapdh = q75_ko_gapdh - q25_ko_gapdh

q75_wt_gapdh, q25_wt_gapdh = np.percentile(wt_gapdh, [75, 25])
iqr_wt_gapdh = q75_wt_gapdh - q25_wt_gapdh

q75_ko_hprt, q25_ko_hprt = np.percentile(ko_hprt, [75, 25])
iqr_ko_hprt = q75_ko_hprt - q25_ko_hprt

q75_wt_hprt, q25_wt_hprt = np.percentile(wt_hprt, [75, 25])
iqr_wt_hprt = q75_wt_hprt - q25_wt_hprt

print("KO IQR (gapdh) = ", iqr_ko_gapdh)
print("WT IQR (gapdh) = ", iqr_wt_gapdh)
print("KO IQR (hprt) = ", iqr_ko_hprt)
print("WT IQR (hprt) = ", iqr_wt_hprt)

print("AFTER OUTLIER REMOVAL")
print("---------------------")

ko_gapdh = [i for i in ko_gapdh if not_outlier(i, iqr_ko_gapdh, q25_ko_gapdh, q75_ko_gapdh)]
wt_gapdh = [i for i in wt_gapdh if not_outlier(i, iqr_wt_gapdh, q25_wt_gapdh, q75_wt_gapdh)]
ko_hprt = [i for i in ko_hprt if not_outlier(i, iqr_ko_hprt, q25_ko_hprt, q75_ko_hprt)]
wt_hprt = [i for i in wt_hprt if not_outlier(i, iqr_wt_hprt, q25_wt_hprt, q75_wt_hprt)]

ko_gapdh_mean = np.mean(ko_gapdh)
ko_hprt_mean = np.mean(ko_hprt)
wt_gapdh_mean = np.mean(wt_gapdh)
wt_hprt_mean = np.mean(wt_hprt)

ko_gapdh_std = np.std(ko_gapdh, ddof=1)
ko_hprt_std = np.std(ko_hprt, ddof=1)
wt_gapdh_std = np.std(wt_gapdh, ddof=1)
wt_hprt_std = np.std(wt_hprt, ddof=1)

print("ko_gapdh: ", ko_gapdh)
print("ko_gapdh mean: ", ko_gapdh_mean)
print("ko_gapdh std: ", ko_gapdh_std)
print("ko_hprt: ", ko_hprt)
print("ko_hprt mean: ", ko_hprt_mean)
print("ko_hprt std: ", ko_hprt_std)
print("wt_gapdh: ", wt_gapdh)
print("wt_gapdh mean: ", wt_gapdh_mean)
print("wt_gapdh std: ", wt_gapdh_std)
print("wt_hprt: ", wt_hprt)
print("wt_hprt mean: ", wt_hprt_mean)
print("wt_hprt std: ", wt_hprt_std)

ddct_gapdh = ko_gapdh_mean - wt_gapdh_mean
ddct_hprt = ko_hprt_mean - wt_hprt_mean
print("ddCT gapdh: ", ddct_gapdh)
print("ddCT hprt: ", ddct_hprt)

sd_ddct_gapdh = (ko_gapdh_std**2 + wt_gapdh_std**2)**.5
sd_ddct_hprt = (ko_hprt_std**2 + wt_hprt_std**2)**.5
print("sd ddct gapdh: ", sd_ddct_gapdh)
print("sd ddct hprt: ", sd_ddct_hprt)

fc_gapdh = ddct_handler(ddct_gapdh)
fc_hprt = ddct_handler(ddct_hprt)
print("fc gapdh: ", fc_gapdh)
print("fc hprt: ", fc_hprt)

ddct_upper_gapdh = ddct_gapdh + sd_ddct_gapdh
ddct_lower_gapdh = ddct_gapdh - sd_ddct_gapdh
ddct_upper_hprt = ddct_hprt + sd_ddct_hprt
ddct_lower_hprt = ddct_hprt - sd_ddct_hprt

print("ddct upper gapdh: ", ddct_upper_gapdh)
print("ddct lower gapdh: ", ddct_lower_gapdh)
print("ddct upper hprt: ", ddct_upper_hprt)
print("ddct lower hprt: ", ddct_lower_hprt)

fc_upper_gapdh = fc_ul_handler(ddct_gapdh, ddct_upper_gapdh)
fc_lower_gapdh = fc_ul_handler(ddct_gapdh, ddct_lower_gapdh)
fc_upper_hprt = fc_ul_handler(ddct_hprt, ddct_upper_hprt)
fc_lower_hprt = fc_ul_handler(ddct_hprt, ddct_lower_hprt)


print("fc upper gapdh: ", fc_upper_gapdh)
print("fc lower gapdh: ", fc_lower_gapdh)
print("fc upper hprt: ", fc_upper_hprt)
print("fc lower hprt: ", fc_lower_hprt)


sd_fc_gapdh = np.std([fc_upper_gapdh,fc_gapdh,fc_lower_gapdh], ddof=1)
sd_fc_hprt = np.std([fc_upper_hprt, fc_hprt, fc_lower_hprt], ddof=1)

print("sd fc gapdh: ", sd_fc_gapdh)
print("sd fc hprt: ", sd_fc_hprt)

std_err_gapdh = std_err(len(ko_gapdh), len(wt_gapdh), sd_fc_gapdh)
std_err_hprt = std_err(len(ko_hprt), len(wt_hprt), sd_fc_hprt)

print("std err gapdh: ", std_err_gapdh)
print("std err hprt: ", std_err_hprt)

gapdh_t_test = stats.ttest_ind(ko_gapdh, wt_gapdh, equal_var=False)
hprt_t_test = stats.ttest_ind(ko_hprt, wt_hprt, equal_var=False)
print("gapdh t test: ", gapdh_t_test[1])
print("hprt t test: ", hprt_t_test[1])

