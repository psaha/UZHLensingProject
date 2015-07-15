#Get list of candidate names in order of decreasing probability of being a lens
import csv

import operator
data = csv.reader(open('candidates.csv'),delimiter = ',')
sortedlist = sorted(data,key=operator.itemgetter(9),reverse=True)

for row in sortedlist:
    print row[8]