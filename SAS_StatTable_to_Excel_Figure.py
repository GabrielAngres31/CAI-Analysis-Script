import csv
import pyperclip as pc

modeSelect = True

mode = ['SBC','Adj R-Sq'][modeSelect] # True|1 for SBC, False|0 for Adj R-sq

DHInitString = "width height diameter thickness WH WD WT HD HT DT WHD WHT WDT HDT WHDT"
dataHeaders = ["intercept " + DHInitString, DHInitString][modeSelect].split(" ")
dataHeadersMap = dict(zip(dataHeaders,range(len(dataHeaders))))

def headerFormatConvert(x):return [x.lower(),x.replace('x','')]['x'in x]

resultTable = {}
resultPasteArray = []

with open('C:/Users/gjang/Documents/SASUniversityEdition/myfolders/PARL0_STATS.csv', newline='') as STATS_P0, \
     open('C:/Users/gjang/Documents/SASUniversityEdition/myfolders/PARL0_STAT-INTS.csv') as INTS_P0:
    S0reader = csv.reader(STATS_P0, delimiter=',', quotechar='\"')
    I0reader = csv.reader(INTS_P0, delimiter=',', quotechar='\"')
    S_discardRow = discardRow = True
    for Srow in S0reader:
        if S_discardRow:
            S_discardRow = False
            continue
        S_accession, S_whichTable, S_value, S_model = Srow
        S_accession, S_value = int(S_accession), float(S_value)
        if S_whichTable == mode:
            if S_accession not in resultTable.keys():
                resultTable[S_accession] = [0]*len(dataHeaders)
            resultTable[S_accession][dataHeadersMap[headerFormatConvert(S_model)]] = S_value
    if not(modeSelect):
        for row in I0reader:
            if discardRow:
                discardRow = False
                continue
            accession, _, value, _ = row
            accession = int(accession)
            value = float(value)
            resultTable[accession][dataHeadersMap[headerFormatConvert("intercept")]] = value
    
    for _,j in resultTable.items():
        resultPasteArray.append(j)
    resultPaste = "\n".join(["\t".join(map(str,i)) for i in resultPasteArray])
    print(resultPaste)
    pc.copy(resultPaste)
                
        
