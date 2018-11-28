from gene_file import gene_files
from gene_stats import gene, subgene
inp1 = input("Enter run 1: ")
inp2 = input("Enter run 2: ")
strt = input("Enter first gene to be considered (0 indexed, so if you want to skip the first you do 1): ")
strt = int(strt)
end = input("Enter the final gene to be considered (the number of the final gene in the spreadsheet like 17 or whatever): ")
end = int(end)
skip = input("Enter any genes you want to skip, separated by a comma (i.e. 3,5 if I want to skip 3 and 5: ")

if len(skip) != 0:
	skip = skip.split(',')
	skip = [int(i) for i in skip]

run1 = gene_files(inp1)
run2 = gene_files(inp2)

# herc_three = gene(subgene("WT"), subgene("KO"),1)

# herc_three.populate_subgenes(run1,run2)
# herc_three.ddct_calculations()
# herc_three.foldchange_calculations()
# herc_three.std_error_calculations()
buf = 'gapdh,,,,,,,hprt\nname,foldchange,error,p-value,,,,name,foldchange,error,p-value\n'
for i in range(strt, end + 1):
	if skip:
		if i in skip:
			continue
	this_gene = gene(subgene("WT"), subgene("KO"), i)
	this_gene.populate_subgenes(run1, run2)
	this_gene.ddct_calculations()
	this_gene.foldchange_calculations()
	this_gene.std_error_calculations()
	this_gene.t_tests()
	buf += this_gene.name + ',' + str(this_gene.fc_gapdh) + ',' + str(this_gene.std_err_gapdh) + ',' + str(this_gene.gapdh_t_test) + ',,,,'
	buf += this_gene.name + ',' + str(this_gene.fc_hprt) + ',' + str(this_gene.std_err_hprt) + ',' + str(this_gene.hprt_t_test) + '\n'

outname = input("Enter the name of the output file (must be .csv): ")

f = open(outname,'w')
f.write(buf)
f.close()

# g1_wt = subgene("WT")
# g1_ko = subgene("KO")

# my_gene1.populate_array(run1.create_slice(1))
# my_gene1.populate_array(run2.create_slice(1))

# print("--------Pre outlier removal---------")
# my_gene1.print_arrays()
# print("--------Post outlier removal--------")
# my_gene1.calc_stats()
# my_gene1.print_arrays()