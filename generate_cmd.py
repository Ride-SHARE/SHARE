temp = []

with open('cmd_template.sh','r') as template:
	for line in template:
		temp.append(line)

fname = 1
prefix = '20180910_'
with open('call_param.out','r') as in_file:
	sub_file = open(prefix+'sub.sh','w')
	for line in in_file:
		if len(line.rstrip()) == 0:
			break
		# line = line.rstrip()
		out_file = open(prefix+str(fname)+'.sh','w')
		for t in temp:
			out_file.write(t)
		out_file.write(line)
		out_file.close()
		sub_file.write('qsub '+prefix+str(fname)+'.sh\n')
		fname += 1
			
	sub_file.close()

