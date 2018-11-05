# -*- coding: utf-8 -*-
import os
import pickle

path = './output/'
listing = os.listdir(path)
globalMp = {};
print(listing)
print(len(listing))

sourceTimeDestination = {} 
for infile in listing:
	# print(infile)
	if infile=='.DS_Store':
		continue
	file1 =  open(path+infile, 'r')
	fileData = file1.readlines()
	index = 0
	for x in range( len( fileData )-1, -1, -1):
		# if ("Stat 600" in fileData[ x ]) or ("Stat 120" in fileData[ x ]) or ("Stat 1080" in fileData[ x ]):
		if ("Stat" in fileData[x]):
			index = x
			break

	# if index == 0:
		# continue
	index=0
	if len(fileData) - index > 30:
		# print("aaaa\n\n\n\n\n")
		continue

	globalMp[ infile ] = {}
	globalMp[ infile ][ "Cabs Dispatched" ] = int( fileData[ index + 1 ].strip("\n").split(" ")[ 3 ] )
	globalMp[ infile ][ "Pickup on way"]  = int( fileData[ index + 2 ].strip("\n").split(" ")[ 4 ] )
	globalMp[ infile ][ "Car pooling"]  = int( fileData[ index + 2 ].strip("\n").split(" ")[ 5 ] )
	globalMp[ infile ][ "Rejection"]  = int( fileData[ index + 3 ].strip("\n").split(" ")[ 2 ] )
	array = {}
	index = index + 3
	array[ 0 ] = []; 	array[ 1 ] = []; 	array[ 2 ] = []
	for y in range( index + 1, len( fileData)-1):		
		if fileData[ y ][0: 14] == "Total distance":
			index = y + 1
			break
		fd = fileData[ y ].strip("\n").split(" ")
		array[ 0 ].append( int( fd[ 3 ] ) )
		array[ 1 ].append( int( fd[ 4 ] ) )
		if( len( fd ) > 5 ):
			array[ 2 ].append( int( fd[ 5 ] ) )
	globalMp[ infile ][ "Waiting" ] = {} 
	globalMp[ infile ][ "Waiting sum" ] = {}
	for k in range(0, 3):
		summ = 0
		for x in range(0, len(array[k]) ):
			summ += x*array[ k ][ x ]
		globalMp[ infile ][ "Waiting" ][ k ] = float(summ)/(sum( array[ k ] ) + 1)
		globalMp[ infile ][ "Waiting sum" ][ k ] = sum( array[ k ] )
		
	mp = {}
	for indexOffset in range( 0, 3):
		passengerK = fileData[ index  + indexOffset ].strip("\n").split(" ")
		mp[ int(passengerK[ 1 ]) ] = float( passengerK[ 3 ])
	index += 3

	if( fileData[ index ][0:9] == "Passenger" ):
		passengerK = fileData[ index ].strip("\n").split(" ")
		mp[ int(passengerK[ 1 ]) ] = float( passengerK[ 3 ] )
		index += 1

	try:
		globalMp[ infile ][ "Distance" ] = mp
	except:
		1
	# print(infile)
	# print(globalMp[infile])
	# print('\n')


	# break

# print(globalMp)
# print(globalMp)
fp1 = open('output.p', 'wb')
pickle.dump(globalMp, fp1)
fp1.close()
