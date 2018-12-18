import openpyxl

wb = openpyxl.load_workbook('geo mean\\set 3 run 2.xlsx', data_only=True)
print(wb.get_sheet_names())
sheet = wb.get_sheet_by_name('GeoMean')

rows = []
for cellObjects in sheet.rows:
	row = []
	for cell in cellObjects:
		row.append(str(cell.value))
	rows.append(row)

print(rows)  