import pandas as pd
import os
import itertools
import math
import xlrd
#xlrd for .xlsx files in jupyter lab

def PandasOutputToMatrix(p, f, df, delim):
    fullpath = os.path.join(p, f)
    print('Output File: %s' % (fullpath))
    df.to_csv(fullpath, sep=delim)

def Log2FoldChange(df, a, b):
    ''' Accepts pandas df, a is index of numerator, b is denominator '''
    ''' Returns pandas series of log2(a/b) fold changes '''
    # divide columns
    div = df.iloc[:, a] / df.iloc[:, b]
    log2FCseries = div.apply(math.log2)
    return(log2FCseries)

def Every3ColumnsForEveryRowCalcMedian(df, x):
    '''finds median of row values for every x columns'''
    '''Assumes first column is row name column'''
    import statistics
    newdf = df.iloc[:, 0]
    newdf = pd.DataFrame(newdf)
    tempdf = df.iloc[:, 1:]  # delete name column, first column # now tempdf is divisible by x
    for i in range(int((tempdf.shape[1])/x)):
        a = i*x
        newcol = []
        for j in range(tempdf.shape[0]):  # for every row, calculate the median of the next three columns from a
            med = statistics.median(tempdf.iloc[j, range(a, a+x)])
            newcol.append(med)
        colname = list(tempdf)[a]
        #print(colname)
        #print(newcol)
        newdf[colname] = newcol
    return(newdf)

def main():
    f = input('Path to .xlsx file please:')
    x = int(input('How many columns do you want to medianize , ex: 3: '))
    #PATH: 'D:\LinuxShare\Projects\Perwin\MASS_SPEC\20190205_MQres_PY.xlsx'

    xl = pd.ExcelFile(f)
    print(xl.sheet_names)

    dfiLFQ = xl.parse('proteinGroups_cleaned', usecols=['Leading_Protein', 'iLFQ EP0_i01', 'iLFQ EP0_i02', 'iLFQ EP0_i03', 'iLFQ EP1_i01', 'iLFQ EP1_i02', 'iLFQ EP1_i03', 'iLFQ EP2_i01', 'iLFQ EP2_i02', 'iLFQ EP2_i03', 'iLFQ EP3_i01', 'iLFQ EP3_i02', 'iLFQ EP3_i03', 'iLFQ EP4_i01', 'iLFQ EP4_i02', 'iLFQ EP4_i03', 'iLFQ EP8_i01', 'iLFQ EP8_i02', 'iLFQ EP8_i03', 'iLFQ EP10_i01', 'iLFQ EP10_i02', 'iLFQ EP10_i03', 'iLFQ EP13_i01', 'iLFQ EP13_i02', 'iLFQ EP13_i03', 'iLFQ EP16_i01', 'iLFQ EP16_i02', 'iLFQ EP16_i03', 'iLFQ EP19_i01', 'iLFQ EP19_i02', 'iLFQ EP19_i03', 'iLFQ EP22_i01', 'iLFQ EP22_i02', 'iLFQ EP22_i03', 'iLFQ EP25_i01', 'iLFQ EP25_i02', 'iLFQ EP25_i03', 'iLFQ EP28_i01', 'iLFQ EP28_i02', 'iLFQ EP28_i03', 'iLFQ EP32_i01', 'iLFQ EP32_i02', 'iLFQ EP32_i03'], skiprows=[1,2,3,4,5,6])
    dfiiTop3 = xl.parse('proteinGroups_cleaned', usecols=['Leading_Protein', 'iiTop3 EP0_i01', 'iiTop3 EP0_i02', 'iiTop3 EP0_i03', 'iiTop3 EP1_i01', 'iiTop3 EP1_i02', 'iiTop3 EP1_i03', 'iiTop3 EP2_i01', 'iiTop3 EP2_i02', 'iiTop3 EP2_i03', 'iiTop3 EP3_i01', 'iiTop3 EP3_i02', 'iiTop3 EP3_i03', 'iiTop3 EP4_i01', 'iiTop3 EP4_i02', 'iiTop3 EP4_i03', 'iiTop3 EP8_i01', 'iiTop3 EP8_i02', 'iiTop3 EP8_i03', 'iiTop3 EP10_i01', 'iiTop3 EP10_i02', 'iiTop3 EP10_i03', 'iiTop3 EP13_i01', 'iiTop3 EP13_i02', 'iiTop3 EP13_i03', 'iiTop3 EP16_i01', 'iiTop3 EP16_i02', 'iiTop3 EP16_i03', 'iiTop3 EP19_i01', 'iiTop3 EP19_i02', 'iiTop3 EP19_i03', 'iiTop3 EP22_i01', 'iiTop3 EP22_i02', 'iiTop3 EP22_i03', 'iiTop3 EP25_i01', 'iiTop3 EP25_i02', 'iiTop3 EP25_i03', 'iiTop3 EP28_i01', 'iiTop3 EP28_i02', 'iiTop3 EP28_i03', 'iiTop3 EP32_i01', 'iiTop3 EP32_i02', 'iiTop3 EP32_i03'], skiprows=[1,2,3,4,5,6])

    print(dfiLFQ.head())
    print(dfiLFQ.shape)
    print(dfiiTop3.head())
    print(dfiLFQ.shape)

    #set first column as row indeces
    #dfiLFQ.index = list(dfiLFQ.iloc[:, 0])
    #dfiiTop3.index = list(dfiiTop3.iloc[:, 0])
    MedDfiLFQ = Every3ColumnsForEveryRowCalcMedian(dfiLFQ, x)
    MedDfiiTop3 = Every3ColumnsForEveryRowCalcMedian(dfiiTop3, x)

    print(type(MedDfiLFQ))
    print(MedDfiLFQ.head())
    print(MedDfiLFQ.shape)
    print(MedDfiiTop3.head())
    print(MedDfiiTop3.shape)

    path, file = os.path.split(f)
    PandasOutputToMatrix(path, 'MedianiLFQ.tsv', MedDfiLFQ, '\t')
    PandasOutputToMatrix(path, 'MedianiiTop3.tsv', MedDfiiTop3, '\t')

    dfoutiLFQ = pd.DataFrame(MedDfiLFQ.iloc[:, 0])
    dfoutiiTop3 = pd.DataFrame(MedDfiiTop3.iloc[:, 0])

    tempiLFQ = MedDfiLFQ.iloc[:, 1:]  # delete name column, first column
    combs = list(itertools.combinations(list(range(tempiLFQ.shape[1])), 2))
    print(combs)
    for t in combs:
        fc = Log2FoldChange(tempiLFQ, t[0], t[1])  # pandas df
        dfoutiLFQ['iLFQ_Log2FC_%dVS%d' % (t[0], t[1])] = fc

    tempiiTop3 = MedDfiiTop3.iloc[:, 1:]  # delete name column, first column
    print(type(tempiiTop3))
    print(tempiiTop3.shape)
    combs = list(itertools.combinations(list(range(tempiiTop3.shape[1])), 2))
    print(combs)
    for t in combs:
        fc = Log2FoldChange(tempiiTop3, t[0], t[1])  # pandas df
        dfoutiiTop3['iiTop3_Log2FC_%dVS%d' % (t[0], t[1])] = fc

    print(dfoutiLFQ.head())
    print(dfoutiLFQ.shape)
    print(dfoutiiTop3.head())
    print(dfoutiiTop3.shape)
    PandasOutputToMatrix(path, 'Log2FCMedianiLFQ.tsv', dfoutiLFQ, '\t')
    PandasOutputToMatrix(path, 'Log2FCMedianiiTop3.tsv', dfoutiiTop3, '\t')



main()

#time points
#index vs DivisionNumber
#0 = 0
#1 = 1
#2 = 2
#3 = 3
#4 = 4
#5 = 8
#6 = 10
#7 = 13
#8 = 16
#9 = 19
#10 = 22
#11 = 25
#12 = 28
#23 = 32