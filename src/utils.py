
 

def formatEdgeId(showppp, outregionA, outregionD):#abbA, abbD):
	edgesetA = dict()
	edgesetD = dict()
	PosA = dict()
	PosD = dict()
	edgeWeightA = dict()
	edgeWeightD = dict()
	edgetoPatient = dict()
	edgeid = 0
	for ppp in showppp:
		if ppp in outregionA:
			for i in range(0, len(outregionA[ppp]), 2):
				edgesetA[edgeid] = list()
				start, end = outregionA[ppp][i], outregionA[ppp][i+1]
				#edgeWeightA[edgeid] = abbA[ppp][i/2]
				edgesetA[edgeid].append(start)
				edgesetA[edgeid].append(end)
				edgetoPatient[edgeid] = ppp				
				if start not in PosA:
					PosA[start] = list()
				PosA[start].append(edgeid)
				if end not in PosA:
					PosA[end] = list()
				PosA[end].append(edgeid)
				edgeid += 1
				
		if ppp in outregionD:
			for i in range(0, len(outregionD[ppp]), 2):
				edgesetD[edgeid] = list()
				#edgeWeightD[edgeid] = abbD[ppp][i/2]
				start, end = outregionD[ppp][i], outregionD[ppp][i+1]
				edgesetD[edgeid].append(start)
				edgesetD[edgeid].append(end)
				edgetoPatient[edgeid] = ppp				
				if start not in PosD:
					PosD[start] = list()
				PosD[start].append(edgeid)
				if end not in PosD:
					PosD[end] = list()
				PosD[end].append(edgeid)
				edgeid += 1	
	return edgesetA, edgesetD, edgetoPatient, edgeWeightA, edgeWeightD, PosA, PosD #, posToEidA, posToEidD


def getCliqueRange(uniqset, edgeset):
	small = edgeset[list(uniqset)[0]][0]
	large = edgeset[list(uniqset)[0]][1]
	tmp = 0
	for eid in uniqset:
		if eid in edgeset:
			if edgeset[eid][0] >= small:
				small = edgeset[eid][0]
			if edgeset[eid][1] <= large:
				large = edgeset[eid][1]
	if large < small:
		tmp = small
		small = large
		large = tmp
	return small, large


def cycle_shift_permutation(regionA, regionD, left, right):
	from collections import deque
	from random import randint

	newA = dict()
	newD = dict()

	for pat in regionA:
		if len(regionA[pat]) > 0:
			#d = deque(regionA[pat])
			rand_i = 2*randint(0, (len(regionA[pat])-1)/2 )
			rand_start = 0
			newA[pat] = list()
			if rand_i == 0:
				#print left, right, regionA[pat][0]-left+1, regionA[pat][-1], right-regionA[pat][-1]+1, left + (regionA[pat][0]-left+1) + (right-regionA[pat][-1]+1)
				rand_start = randint(left, left + (regionA[pat][0]-left+1) + (right-regionA[pat][-1]+1)) 
			else:
				rand_start = randint(left, left + regionA[pat][rand_i] - regionA[pat][rand_i-1] + 1)

			 
			rand_range = regionA[pat][rand_i] - rand_start + 1
			#print rand_i, rand_start, rand_range
			for i in range(rand_i, len(regionA[pat])):
				newA[pat].append(regionA[pat][i] - rand_range + 1)

			for i in range(0, rand_i):
				newA[pat].append(regionA[pat][i] + (right -left) - rand_range + 1)
			
			#print pat, regionA[pat]
			#print pat, newA[pat]

	for pat in regionD:
		if len(regionD[pat]) > 0:		
			rand_i = 2*randint(0, (len(regionD[pat])-1)/2 )
			rand_start = 0
			newD[pat] = list()
			if rand_i == 0:
				rand_start = randint(left, left + (regionD[pat][0]-left+1) + (right-regionD[pat][-1]+1)) 
			else:
				rand_start = randint(left, left + regionD[pat][rand_i] - regionD[pat][rand_i-1] + 1)

			 
			rand_range = regionD[pat][rand_i] - rand_start + 1
			
			for i in range(rand_i, len(regionD[pat])):
				newD[pat].append(regionD[pat][i] - rand_range + 1)

			for i in range(0, rand_i):
				newD[pat].append(regionD[pat][i] + (right -left) - rand_range + 1)
			
			#print pat, regionD[pat]
			#print pat, newD[pat]
	return newA, newD


#testA = dict()
#testA['pat1'] = [200,300,500,700,900,950]
#cycle_shift_permutation(testA, "", 100, 1000)


def collectGenes(metagenebkt, leftbound, rightbound):
	gene_in = set()
	for g in metagenebkt:
		if not (metagenebkt[g][0] > rightbound or metagenebkt[g][1] < leftbound):
			gene_in.add(g)	
	return gene_in


def sort_func_with_same_pos(posList):
	out = list()
	for pos in sorted(posList.keys()):
		out += [(pos, posList[pos])]
	return out
	# from itertools import groupby
	# from operator import itemgetter
	# out = list()
	# seq = posList
	# seq.sort(key = itemgetter(0))
	# groups = groupby(seq, itemgetter(0))
	# for (key, data) in groups:
	# 	out += [(key, [item[1] for item in data])]
	# return out