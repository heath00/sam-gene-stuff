import csv

class gene_files:
	def __init__(self, fname):
		with open(fname) as csv_file:
			reader = csv.reader(csv_file, delimiter=',')
			self.my_file = list(reader)

	def create_slice(self, gene_num):
		start_line = gene_num*17
		return self.my_file[start_line:start_line + 17]
	
	def get_contents(self):
		return self.my_file
	
	def print_contents(self):
		print(self.my_file)


