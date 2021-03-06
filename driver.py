from gene_file import gene_files
from gene_stats import gene, subgene
import os
choice = input("Want to enter a directory? [y/n]: ")
if choice == 'n':
	inp = input("Enter run filenames separated by a comma (i.e. run 1.csv,run 2.csv): ")
	inp = inp.split(',')
else:
	inp = input("Enter directory name: ")
	inpfiles = os.listdir(inp)
	inpfiles = [x for x in inpfiles if x[-1] == 'x' or x[-1] == 'v']
	inpfiles = [inp + '/' + x for x in inpfiles]
	inp = inpfiles




strt = input("Enter first gene to be considered (0 indexed, so if you want to skip the first you do 1): ")
strt = int(strt)
end = input("Enter the final gene to be considered (the number of the final gene in the spreadsheet like 17 or whatever): ")
end = int(end)
skip = input("Enter any genes you want to skip, separated by a comma (i.e. 3,5 if I want to skip 3 and 5: ")

if len(skip) != 0:
	skip = skip.split(',')
	skip = [int(i) for i in skip]


print(inp)
runs = []
geomean_runs = []
for run in inp:
	this_run = gene_files(run)
	runs.append(this_run)
	geomean_run = gene_files(run, True)
	geomean_runs.append(geomean_run)


buf = 'gapdh,,,,,,,,,hprt,,,,,,,,,geomean\nname,foldchange,foldchange upper bound,foldchange lower bound,p-value,KO #,WT #,,,,name,foldchange,foldchange upper bound,foldchange lower bound,p-value,KO #,WT #,,,,name,foldchange,foldchange upper bound,foldchange lower bound,p-value,KO #,WT #\n'
for i in range(strt, end + 1):
	if skip:
		if i in skip:
			continue
	this_gene = gene(subgene("WT"), subgene("KO"), i)
	this_gene.populate_subgenes(runs)
	this_gene.ddct_calculations()
	this_gene.foldchange_calculations()
	this_gene.std_error_calculations()
	this_gene.t_tests()

	geomean_gene = gene(subgene("WT"), subgene("KO"), i)
	geomean_gene.populate_subgenes(geomean_runs)
	geomean_gene.ddct_calculations()
	geomean_gene.foldchange_calculations()
	geomean_gene.std_error_calculations()
	geomean_gene.t_tests()
	post_outlier_numbers = {}
	this_gene.get_post_outlier_numbers(post_outlier_numbers)
	post_outlier_numbers_geomean = {}
	geomean_gene.get_post_outlier_numbers(post_outlier_numbers_geomean)
	buf += this_gene.name + ',' + str(this_gene.fc_gapdh) + ',' + str(this_gene.fc_upper_gapdh) + ',' + str(this_gene.fc_lower_gapdh) + ',' + str(this_gene.gapdh_t_test) + ',' + str(post_outlier_numbers['gapdh_ko']) + ',' + str(post_outlier_numbers['gapdh_wt']) + ',,,,'
	buf += this_gene.name + ',' + str(this_gene.fc_hprt) + ',' + str(this_gene.fc_upper_hprt) + ',' + str(this_gene.fc_lower_hprt) + ',' + str(this_gene.hprt_t_test) + ',' + str(post_outlier_numbers['hprt_ko']) + ',' + str(post_outlier_numbers['hprt_wt']) + ',,,,'
	buf += geomean_gene.name + ',' + str(geomean_gene.fc_gapdh) + ',' + str(geomean_gene.fc_upper_gapdh) + ',' + str(geomean_gene.fc_lower_gapdh) + ',' +  str(geomean_gene.gapdh_t_test) + ',' + str(post_outlier_numbers_geomean['gapdh_ko']) + ',' + str(post_outlier_numbers_geomean['gapdh_wt']) + '\n'
	print(str(post_outlier_numbers_geomean['hprt_ko']) + "and" + str(post_outlier_numbers_geomean['hprt_wt']))

outname = input("Enter the name of the output file (must be .csv): ")

f = open(outname,'w')
f.write(buf)
f.close()