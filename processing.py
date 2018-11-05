import sys

file1 = open(sys.argv[1])

scoreDAG = 0
scoreDAGEX = 0
scoreHuer = 0
timeDAG = 0
timeDAGEX = 0
timeHuer = 0
maxO = 200
counter = 0
for x in range(0, 10000):
	print x
	if( x >= maxO ):
		break
	loopOver = 0
	while(1):	
		while (True):
			temp = file1.readline()
			if( temp == ""):
				loopOver = 1
				break
			temp = temp.strip("\n")
			
			token = temp.split("\t")
			if (token[0] == 'Trip' and len(token) >= 2):
				break
		if( loopOver ):
			break
		score = file1.readline().strip("\n").split("\t")[1:]
		time = file1.readline().strip("\n").split("\t")[1:]
		print score
		if( score[1] != '0'):
			break
	if( score[1] == '0'):
		continue
	if( loopOver ):
		break
	counter += 1
	print counter
	scoreDAG += min( (float( score[1] ) + 1 )/( int( score[0] ) + 1), 15)
	timeDAG += (float( time[1]) + 0.0001)
	scoreDAGEX += min((float( score[2] ) + 1 )/(int( score[0] ) + 1), 15 )
	scoreHuer += min((float( score[3] ) + 1 )/(int( score[0] ) + 1), 15 )
	if( float( time[2]) + 0.0001 > 100 ):
		timeDAGEX += 10
	else:
		timeDAGEX += float( time[2] )
	if( float( time[3]) + 0.0001 > 100 ):
                timeHuer += 20
        else:
                timeHuer += float( time[3] )
maxO = counter 
print scoreDAG/maxO, timeDAG/maxO
print scoreDAGEX/maxO, timeDAGEX/maxO
print scoreHuer/maxO, timeHuer/maxO
