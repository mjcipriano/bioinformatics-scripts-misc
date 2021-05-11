#!/usr/bin/python


import pandas as pd
import numpy as np


import glob
import re
import sys
import os
import math

# import easygui
# plotly.tools.set_credentials_file(username='mcipriano', api>>> plotly.tools.set_credentials_file(username='mcipriano', api_key='')
from plotly import tools
import plotly.plotly as py
import plotly.graph_objs as go

colNames=range(1,13)
rowNames = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
colString = "ABCDEFGH"


myFiles = glob.glob("./*p4*.xlsx"); # list of files in current directory. Change this to suit your needs.
myValidFiles = list() # files which are valid will be placed in here
allData = {}
wellData = {}

for thisExcelFile in myFiles :
    if re.search('~', thisExcelFile): # temporary files usually have a tilde in them, this skips them
        print "Skipping ", thisExcelFile, "\n"
        pass
    else:
        print "Reading ", thisExcelFile, "\n"
        thisResultFile = pd.read_excel(thisExcelFile, 'Results') # Read in the file
	thisResultName = os.path.basename(thisExcelFile) # key used for accessing this particular file
	myValidFiles.append(thisResultName) # for iterating through the files in order
	plateData= thisResultFile.iloc[8:17, 5:17] # Get the slice of results from the excel file
	indexedData = pd.DataFrame(plateData.as_matrix(), index=rowNames, columns=colNames) # create an indexed set of this data
	# print indexedData
	allData[thisResultName] = indexedData

valMax = 0
valMin = 99999999999999
allWellNames = []
for tRow in rowNames:
	for tCol in colNames:
		sys.stdout.write(tRow)
		sys.stdout.write("%s" % tCol)
		sys.stdout.write(' ')
		allWellNames.append(tRow + str(tCol))
		wellData[tRow + str(tCol)] = []
		for fn in myValidFiles:
			sys.stdout.write("%s" % allData[fn].get_value(tRow, tCol) )
			sys.stdout.write(',')
			thisWellData = allData[fn].get_value(tRow, tCol)
			wellData[tRow + str(tCol)].append(thisWellData)
			if(valMax < thisWellData):
				valMax = thisWellData
			if(valMin > thisWellData):
				valMin = thisWellData
		print ""
print "Done\n"


# Ploting Well Level Data Start

myFig = tools.make_subplots(rows=8, cols=12, subplot_titles=allWellNames)
myFig['layout'].update(height=1024,width=1280, title='Well Data')
rowNum = 0
numPlots = 0
for tRow in rowNames:
	rowNum += 1
	for tCol in colNames:
		datapoints = len(wellData[tRow + str(tCol)])
		trace = go.Scatter(
			x=range(1, datapoints+1),
			y=wellData[tRow + str(tCol)],
			name=tRow + str(tCol)
		)
		myFig.append_trace(trace, rowNum, tCol)
		numPlots += 1
		print "Loading plot ", tRow, tCol, "plot ", rowNum, tCol
		
			
#print wellData["A1"]
for i in range(1,numPlots+1):
	myFig['layout']['yaxis' + str(i)].update(range=[0,valMax], autorange=False)

# uncomment when ready to update well data TMP
py.iplot(myFig, filename='wellData')


# END Plotting Well Level Data



# Start Replicates data

rc = pd.read_csv("replicates.txt", sep='\t', header=None)
replicatesConfig = pd.DataFrame(rc.as_matrix(), index=rowNames, columns=colNames)
# print replicatesConfig

replicatesData = {}
for tRow in rowNames:
	for tCol in colNames:
		# Get Well Sample Name
		tWellSampleName = replicatesConfig.get_value(tRow, tCol)
		if tWellSampleName not in replicatesData:
			replicatesData[tWellSampleName] = {}
		sys.stdout.write(tRow)
		sys.stdout.write("%s" % tCol)
		sys.stdout.write(' ')
		for fn in myValidFiles:
			sys.stdout.write("%s" % allData[fn].get_value(tRow, tCol) )
			sys.stdout.write(',')
			thisWellData = allData[fn].get_value(tRow, tCol)
			if fn not in replicatesData[tWellSampleName]:
				replicatesData[tWellSampleName][fn] = []
			replicatesData[tWellSampleName][fn].append(thisWellData)
		print ""
# print replicatesData

# Load Extended Config Data

IDnamesOrdered = []
DrugNamesOrdered = []
DrugNamesDict = {}
rex = pd.read_csv("types.txt", sep='\t')

lastD = ""
numrows = 0
numcols = 1
tCol = 1
IDnamesFig = []

for index, row in rex.iterrows():
	# Create single graph
	print row["ID"], row["DrugName"], row["Concentration"]
	IDnamesOrdered.append(row["ID"])
	if row["DrugName"] not in DrugNamesDict:
		DrugNamesDict[row["DrugName"]] = 1
		DrugNamesOrdered.append(row["DrugName"])
	if(lastD == row["DrugName"]):
		tCol += 1
		if(tCol > numcols):
			numcols = tCol
	else:
		numrows += 1
		tCol = 1
	lastD = row["DrugName"]

# redo for titles of graph
lastD = ""
tRow = 0
tCol = 1
for index, row in rex.iterrows():
	if(lastD == row["DrugName"]):
		tCol += 1
	elif(lastD == ""):
		pass
	else:
		numAdd = numcols-tCol
		for i in range(0,numAdd):
			IDnamesFig.append(" ")
		tCol = 1
	IDnamesFig.append(row["ID"])
	lastD = row["DrugName"]


thisWidth = 400*numcols
thisHeight = 500*numrows

print IDnamesFig
myCombFig = tools.make_subplots(rows=numrows, cols=numcols, subplot_titles=IDnamesFig)
myCombFig['layout'].update(height=thisHeight,width=thisWidth, title='Drug Data')

myMultiFig = tools.make_subplots(rows=len(DrugNamesOrdered), cols=1, subplot_titles=DrugNamesOrdered)
myMultiFig['layout'].update(height=thisHeight,width=600, title='Drug Data')


numrow = 0
numcol = 1
lastD = ""
for index, row in rex.iterrows():
	thisGraphAvg = []
	thisGraphStd = []
	datapoints = len(myValidFiles)
	for fn in myValidFiles:
		thisData = replicatesData[row["ID"]][fn]
		print "  ", fn
		print "      ", thisData
		print "      ", "AVG", np.average(thisData), "STDDEV0", np.std(thisData, ddof=0), "STDDEV1", np.std(thisData, ddof=1)
		thisGraphAvg.append(np.average(thisData))
		thisGraphStd.append(np.std(thisData, ddof=1))

	if(lastD == row["DrugName"]):
		numcol += 1
		print "increment col", numcol
	else:
		numrow += 1
		numcol = 1
		print "increment row set col 1", numcol, numrow
	lastD = row["DrugName"]
	trace = go.Scatter(
		x=range(1, datapoints+1),
		y=thisGraphAvg,
		name=row["ID"] + "-" + row["DrugName"] + "-" + str(row["Concentration"]),
		error_y=dict(
			type='data',
			array=thisGraphStd,
			visible=True
		)
	)
	myCombFig.append_trace(trace, numrow, numcol)
	myMultiFig.append_trace(trace, numrow, 1)

totPlots = numrows*numcols
for i in range(1,totPlots+1):
	myCombFig['layout']['yaxis' + str(i)].update(range=[0,valMax], autorange=False)

for i in range(1,len(DrugNamesOrdered)+1):
	myMultiFig['layout']['yaxis' + str(i)].update(range=[0,valMax], autorange=False)

#print myCombFig
# uncomment when ready to update well data TMP
py.iplot(myCombFig, filename='drugData')
py.iplot(myMultiFig, filename='drugMultiData')

print "Done\n"
