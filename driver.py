from gene_file import gene_files
from gene_stats import gene, subgene
inp = input("Enter run filenames separated by a comma (i.e. run 1.csv,run 2.csv): ")
strt = input("Enter first gene to be considered (0 indexed, so if you want to skip the first you do 1): ")
strt = int(strt)
end = input("Enter the final gene to be considered (the number of the final gene in the spreadsheet like 17 or whatever): ")
end = int(end)
skip = input("Enter any genes you want to skip, separated by a comma (i.e. 3,5 if I want to skip 3 and 5: ")

if len(skip) != 0:
	skip = skip.split(',')
	skip = [int(i) for i in skip]

inp = inp.split(',')
print(inp)
runs = []
geomean_runs = []
for run in inp:
	this_run = gene_files(run)
	runs.append(this_run)
	geomean_run = gene_files(run, True)
	geomean_runs.append(geomean_run)


buf = 'gapdh,,,,,,,hprt,,,,,,,geomean\nname,foldchange,error,p-value,,,,name,foldchange,error,p-value,,,,name,foldchange,error,p-value\n'
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
	buf += this_gene.name + ',' + str(this_gene.fc_gapdh) + ',' + str(this_gene.std_err_gapdh) + ',' + str(this_gene.gapdh_t_test) + ',,,,'
	buf += this_gene.name + ',' + str(this_gene.fc_hprt) + ',' + str(this_gene.std_err_hprt) + ',' + str(this_gene.hprt_t_test) + ',,,,'
	buf += geomean_gene.name + ',' + str(geomean_gene.fc_gapdh) + ',' + str(geomean_gene.std_err_gapdh) + ',' + str(geomean_gene.gapdh_t_test) + '\n'

outname = input("Enter the name of the output file (must be .csv): ")

f = open(outname,'w')
f.write(buf)
f.close()