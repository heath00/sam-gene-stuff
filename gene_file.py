import csv
import openpyxl
import os

class gene_files:
	def __init__(self, fname, isgeomean=False):

		# if it's an excel file
		if fname[-1] == 'x':
			# wb = openpyxl.load_workbook(fname, data_only=True)

			# print('These sheets are available in ' , fname , ':', wb.get_sheet_names())
			# sheet_ind = input('Enter the index of the sheet you want to access (first one is 0, next is 1, etc.): ')
			# sheet = wb.worksheets[int(sheet_ind)]

			# rows = []
			# for cellObjects in sheet.rows:
			# 	row = []
			# 	for cell in cellObjects[1:]:
			# 		this_val = cell.value
			# 		if this_val:
			# 			row.append(str(this_val))
			# 		else:
			# 			row.append('')
			# 	rows.append(row)
			# self.my_file = rows
			# print(len(self.my_file))
			wb = openpyxl.load_workbook(fname, data_only=True)

			if not isgeomean:
				# print('These sheets are available in ' , fname , ':', wb.get_sheet_names())
				# sheet_ind = input('Enter the index of the sheet you want to access (first one is 0, next is 1, etc.): ')
				sheet_ind = wb.get_sheet_names().index('Analysis')
				print("analysis index is: ", sheet_ind)
			else:
				sheet_ind = wb.get_sheet_names().index('GeoMean')
				print("geomean index is: ", sheet_ind)


			sheet = wb.worksheets[int(sheet_ind)]

			with open(fname[:-4] + 't1.csv', 'w', newline="") as f:
				c = csv.writer(f)
				for r in sheet.rows:
					c.writerow([cell.value for cell in r[1:]])

			with open(fname[:-4] + 't1.csv') as csv_file:
				reader = csv.reader(csv_file, delimiter=',')
				self.my_file=list(reader)

			os.remove(fname[:-4] + 't1.csv')
		# else we treat it as csv	
		else: 
			with open(fname) as csv_file:
				reader = csv.reader(csv_file, delimiter=',')
				self.my_file = list(reader)
			print("length:", len(self.my_file))

	def create_slice(self, gene_num):
		start_line = gene_num*17
		return self.my_file[start_line:start_line + 17]
	
	def get_contents(self):
		return self.my_file
	
	def print_contents(self):
		print(self.my_file)

	def print_len(self):
		print(len(self.my_file))


