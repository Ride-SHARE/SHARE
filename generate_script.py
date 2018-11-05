
import numpy as np

# all
# algos = [(0.1,'dij'),(0.1,'dag'),(0.1,'dex'),(0.3,'dex'),(0.5,'dex')]
algos = [(0.1,'dij'),(0.1,'dag'),(0.1,'dex'),(0.3,'dex')]
assign = 'pxa'
cab_caps = [3]
# variables
alphas = [1.1,1.2,1.3,1.5,1.8,2.1]
alphas_def = 1
loc_dict = {'SG': [2000,4000,6000,8000,10000], 'NY': [2000,4000,6000,8000,10000]}
end_dict = {'SG': 120, 'NY': 120}
cab_def = {'SG': 8000,'NY': 8000}
start_times = [(0,'12am'),(480,'8am'),(960,'4pm')]
start_def = 1

param_set = set()

script_name = 'a.out'
file_prefix = '20180917'

with open('call_param.out','w') as out_file:
	for loc in loc_dict.keys():
		for algo in algos:
			for cab_cap in cab_caps:
				for alpha in alphas:
					num_cab = cab_def[loc]
					start_time = start_times[start_def]
					cmd_str = './' +script_name+ ' '+loc+' '+str(alpha)+' '+str(algo[1])+' '+str(algo[0])+' '+assign+' '+str(num_cab)+' '+str(cab_cap)+' '+str(start_time[0])+' '+str(start_time[0]+end_dict[loc])
					file_str = file_prefix + '/'+loc+'_'+assign+'_'+algo[1]+str(algo[0])+'_a'+str(alpha)+'_cap'+str(cab_cap)+'_n'+str(num_cab)+'_'+start_time[1]
					if not file_str in param_set:
						out_file.write(cmd_str + ' > ' + file_str + '\n')
						param_set.add(file_str)
	for loc in loc_dict.keys():
		for algo in algos:
			for cab_cap in cab_caps:
				for num_cab in loc_dict[loc]:
					alpha = alphas[alphas_def]
					start_time = start_times[start_def]
					cmd_str = './' +script_name+ ' '+loc+' '+str(alpha)+' '+str(algo[1])+' '+str(algo[0])+' '+assign+' '+str(num_cab)+' '+str(cab_cap)+' '+str(start_time[0])+' '+str(start_time[0]+end_dict[loc])
					file_str = file_prefix + '/'+loc+'_'+assign+'_'+algo[1]+str(algo[0])+'_a'+str(alpha)+'_cap'+str(cab_cap)+'_n'+str(num_cab)+'_'+start_time[1]
					if not file_str in param_set:
						out_file.write(cmd_str + ' > ' + file_str + '\n')
						param_set.add(file_str)
	for loc in loc_dict.keys():
		for algo in algos:
			for cab_cap in cab_caps:
				for start_time in start_times:
					alpha = alphas[alphas_def]
					num_cab = cab_def[loc]
					cmd_str = './' +script_name+ ' '+loc+' '+str(alpha)+' '+str(algo[1])+' '+str(algo[0])+' '+assign+' '+str(num_cab)+' '+str(cab_cap)+' '+str(start_time[0])+' '+str(start_time[0]+end_dict[loc])
					file_str = file_prefix + '/'+loc+'_'+assign+'_'+algo[1]+str(algo[0])+'_a'+str(alpha)+'_cap'+str(cab_cap)+'_n'+str(num_cab)+'_'+start_time[1]
					if not file_str in param_set:
						out_file.write(cmd_str + ' > ' + file_str + '\n')
						param_set.add(file_str)
