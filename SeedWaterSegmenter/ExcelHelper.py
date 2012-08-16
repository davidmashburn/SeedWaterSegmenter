import xlwt
import xlrd

def excelRead(filename,flip=False):
    wb = xlrd.open_workbook(filename)
    sheets = wb.sheets()
    sheetnames = [i.name for i in sheets]
    data=[] # 3d data set... sheets, rows, columns...
    for sheet in sheets:
        #print sheet.name
        data.append([]) #make a data store
        if flip:
            for i in range(sheet.ncols):
                data[-1].append(sheet.col_values(i))
        else:
            for i in range(sheet.nrows):
                data[-1].append(sheet.row_values(i))
    return data,sheetnames
def excelWrite(data,sheetnames,filename,flip=False):
    wb=xlwt.Workbook()
    for i,name in enumerate(sheetnames):
        ws=wb.add_sheet(name)
        for j in range(len(data[i])):
            for k in range(len(data[i][j])):
                if flip:
                    ws.write(k,j,data[i][j][k])
                else:
                    ws.write(j,k,data[i][j][k])
    wb.save(filename)

if __name__=='__main__':
    import numpy as np
    a=np.array([[1,2,3,4],
                [5,6,7,8],
                [9,10,11,12]])
    excelWrite([a],['NumericalListsTo12'],'/home/mashbudn/Desktop/testxlwt.xls')
    print excelRead('/home/mashbudn/Desktop/testxlwt.xls')
