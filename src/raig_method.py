
import operator
import re
from cliqueinterval import CliqueInterval
		
def infopack(sum_left, sum_right, sum_left_pat, sum_right_pat, delta, pivot, cons_cid, wleft, wright):
	return dict({'both': sum_left.intersection(sum_right),
				'left': sum_left, 
				'right': sum_right, 
				'cons_cid': cons_cid,
				'ucount': len(sum_left.intersection(sum_right)), 
				'rcount':len(sum_right), 
				'lcount':len(sum_left), 
				'ucount_pat': len(sum_left_pat.intersection(sum_right_pat)), 
				'rcount_pat': len(sum_right_pat), 				
				'lcount_pat': len(sum_left_pat), 				
				'p-val' : 1,
				'pivot' : pivot,
				'localc': str(delta)+'_soft'} )

def allowWeight(delta, C, i, j, clist, cnatype):

	cscore = 0
	cscore_ind = -1
	uscore = 0

	return_info = {}
	# for each block, finding a pivot clique with maximum score by excluding noises in all the rest maximal cliques
	for ll in range(j, i+1): # ll : pivot				
		sum_left = set()
		sum_right = set()

		sum_left, sum_right = C.getConsecutiveCliqueInfo_Fast(j, i, ll, cnatype)
		
		right_sample = C.returnSampleLevel(sum_right)
		left_sample = C.returnSampleLevel(sum_left)
			
		WI = 2*min(len(left_sample), len(right_sample))
		
		if WI >= delta:
			if WI > cscore:
				cscore = WI
				cscore_ind = ll
				uscore = len(left_sample.intersection(right_sample))
				return_info = infopack(sum_left, sum_right, left_sample, right_sample, delta, ll, clist[j:i+1], C.getRegion(j,cnatype)[0], C.getRegion(i,cnatype)[1])
			if WI == cscore and uscore <len(left_sample.intersection(right_sample)): # tie, check unique intervals
				cscore = WI
				cscore_ind = ll
				uscore = len(left_sample.intersection(right_sample))
				return_info = infopack(sum_left, sum_right, left_sample, right_sample, delta, ll, clist[j:i+1], C.getRegion(j,cnatype)[0], C.getRegion(i,cnatype)[1])
		
	return cscore, cscore_ind, return_info
	

	
def getWeight( delta, C, i, j, clist, cnatype):
	
	if len(clist) == 1: # for only one clique in genome
		tmp_x, tmp_y, tmp_z = C.getConsecutiveCliqueInfo(j, i - j + 1, cnatype)		
		tmp_x_pat = C.returnSampleLevel(tmp_x)
	
		WI = 2*len(tmp_x_pat)
		if WI >= delta:
			return WI, -1, infopack(tmp_x, tmp_x, tmp_x_pat, tmp_x_pat, delta, -1, clist[j:i+1], C.getRegion(j,cnatype)[0], C.getRegion(i,cnatype)[1])
		else:
			return 0, -1, {}
	else: # multiple cliques in genome
		if i == 0 or i == j: # marginal
			tmp_x, tmp_y, tmp_z = C.getConsecutiveCliqueInfo(j, i - j + 1, cnatype)		
			tmp_x_pat = C.returnSampleLevel(tmp_x)
			tmp_y_pat = C.returnSampleLevel(tmp_y)
			tmp_z_pat = C.returnSampleLevel(tmp_z)

			WI = 2*len(tmp_x_pat)+2*min(len(tmp_y_pat),len(tmp_z_pat))

			if WI >= delta:
				return WI, -1, infopack(tmp_x.union(tmp_z), tmp_x.union(tmp_y), tmp_x_pat.union(tmp_z_pat), tmp_x_pat.union(tmp_y_pat), delta, -1, clist[j:i+1], C.getRegion(j,cnatype)[0], C.getRegion(i,cnatype)[1])
			else:
				return 0, -1, {}
		else:
			soft_WI, soft_i, return_info = allowWeight(delta, C, i, j, clist, cnatype)
			if soft_WI >= delta:
				return soft_WI, soft_i, return_info
			else:
				return 0, -1, {}				

def dynamicProg_RAIG (C, focalRegion, delta, cnatype, edgew):
	opt = dict() # use to store optimal score
	softInd = dict()
	outinfo = dict()
	output = dict()
	unused = dict() # use to store cliques without passing delta (debugging)
	debug_storage = dict()
	lastchange = dict() # use to store dp step
	if cnatype == 'A':
		clist = C.cidOrderListA
	else:
		clist = C.cidOrderListD

	if clist:
		opt[0] = 0			
		starti, endi = 0, len(clist) #getStartEnd(clist)
		for i in range(starti, endi):		
			tmp_max = -1
			tmp_index = i
			tmp_soft = -1
			tmp_info = {}			

			if i == 0:		
				tmp_max, soft, info = getWeight(delta, C, i, i, clist, cnatype)
				tmp_index = i
				tmp_soft = soft
				tmp_info = info

			else:
				for j in range(i, starti-1, -1 ): 
					if C.getRegion(i, cnatype)[0] - C.getRegion(j, cnatype)[0] > focalRegion:				
						break
					else:
						rscore, soft, info = getWeight(delta, C, i, j, clist, cnatype)							
						tmp_score = opt[max(j-1,0)] + rscore
											
						if tmp_max < tmp_score:
							tmp_max = tmp_score
							tmp_index = j	
							tmp_soft = soft
							tmp_info = info

			opt[i] = max(tmp_max, 0)
			lastchange[i] = tmp_index
			softInd[i] = tmp_soft
			outinfo[i] = tmp_info
		
		ls_index = endi - 1			
		while 1:				
			r = ls_index - lastchange[ls_index] + 1
			i = lastchange[ls_index]

			if ls_index == 0:
				WI = opt[ls_index]
			else:
				WI = opt[ls_index] - opt[max(i-1, 0)]					
			if WI > 0:				
				output[clist[i]] = outinfo[ls_index]			
			ls_index = lastchange[ls_index] - 1
			if ls_index not in lastchange:
				break
	
	return output 


