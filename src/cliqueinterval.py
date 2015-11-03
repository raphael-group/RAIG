#!/usr/bin/env python

class CliqueInterval():
	
	idata = dict() # interval
	cnatype = dict() # type	
	region = dict()
	genes = dict() # genes in the clique
	N_left = dict() # left endpoint
	N_right = dict() # right endpoint
	N_p_left = dict() # left endpoint
	N_p_right = dict() # right endpoint
	N_p = dict()

	output = dict()
	cidOrderListA = list()
	cidOrderListD = list()
	segDegree = dict()
	edgetoPatient = dict()
	avgFreq = dict()
	previous_record = dict()
	
	def __init__(self):
		self.idata = dict() # interval
		self.cnatype = dict() # type
		self.region = dict() # region
		self.genes = dict() # genes in the clique
		self.N_left = dict() # left endpoint
		self.N_right = dict() # right endpoint
		self.N_p_left = dict() # left endpoint
		self.N_p_right = dict() # right endpoint
		self.N_p = dict() # segments in the clique
		self.cidOrderListA = list()
		self.cidOrderListD = list()
		self.segDegree = dict()
		self.edgetoPatient = dict()
		self.avgFreq = dict()
		self.previous_record = dict()
	

	def getAvgFreq(self):
		for i in range(len(self.cidOrderListA)):
			cid = self.cidOrderListA[i]
			self.avgFreq[cid] = [self.region[cid][0], self.region[cid][0], len(self.idata[cid])]

		for i in range(len(self.cidOrderListD)):
			cid = self.cidOrderListA[i]
			self.avgFreq[cid] = [self.region[cid][0], self.region[cid][0], len(self.idata[cid])]

	def orderClique(self):
		for k in sorted(filter(lambda x:self.cnatype[x]=='A', self.region.keys()), key=lambda x:self.region[x][0]):
			self.cidOrderListA.append(k)
		for k in sorted(filter(lambda x:self.cnatype[x]=='D', self.region.keys()), key=lambda x:self.region[x][0]):
			self.cidOrderListD.append(k)		


	def getRegion(self, i, cnatype):
		if cnatype == 'A':
			clist = self.cidOrderListA
		else:
			clist = self.cidOrderListD

		return self.region[clist[i]][0], self.region[clist[i]][1]


	def createBiGraph(self, cliq_results, cutoffA, edgesetA, cutoffD, edgesetD):
		import utils as cna_utils
		cliq_id = 0				
		for allc in cliq_results:			
			for i in allc: # eid
				if i not in self.segDegree:
					self.segDegree[i] = set()
				self.segDegree[i].add(cliq_id)				
				
			self.idata[cliq_id] = set(allc)

			if int(list(allc)[0]) in edgesetA and len(allc) >= round(cutoffA): # remove clique with too few patients
				self.cnatype[cliq_id] = 'A'
				self.region[cliq_id] = cna_utils.getCliqueRange(self.idata[cliq_id], edgesetA)
				
			if int(list(allc)[0]) in edgesetD and len(allc) >= round(cutoffD): # remove clique with too few patients
				self.cnatype[cliq_id] = 'D'				
				self.region[cliq_id] = cna_utils.getCliqueRange(self.idata[cliq_id], edgesetD)
				
			cliq_id += 1

		self.orderClique()	

	def formatOutput(self, i, ls_index, r, clist, delta, cnatype, soft):		

		if soft != -1:			
			sum_left_pat = set()
			sum_left = set()
			sum_right_pat = set()
			sum_right = set()

			for s in range(i, soft+1):
				cross_x, cross_y, cross_z, cross_all = self.getConsecutiveCliqueInfo(s, soft-s + 1, cnatype)
				sum_left.update(cross_z.union(cross_x))				

			for e in range(soft, ls_index+1):
				cross_x, cross_y, cross_z, cross_all = self.getConsecutiveCliqueInfo(soft, e-soft + 1, cnatype)
				sum_right.update(cross_y.union(cross_x))				

			sum_left_pat = self.returnSampleLevel(sum_left)
			sum_right_pat = self.returnSampleLevel(sum_right)

			# for debugging
			return {'both': sum_left.intersection(sum_right),
				'left': sum_left, 
				'right': sum_right, 
				'cons_cid': clist[i:i+1], 
				'wide_region': self.getRegion(i,cnatype),
				'ucount': len(sum_left.intersection(sum_right)), 
				'rcount':len(sum_right), 
				'lcount':len(sum_left), 
				'ucount_pat': len(sum_left_pat.intersection(sum_right_pat)), 
				'rcount_pat': len(sum_right_pat), 				
				'lcount_pat': len(sum_left_pat), 				
				'p-val' : 1,
				'localc': str(delta)+'_soft'} 
		
		else:
			tmp_x, tmp_y, tmp_z, tmp_all = self.getConsecutiveCliqueInfo(i, r, cnatype)		
			tmp_x_pat = self.returnSampleLevel(tmp_x)
			tmp_y_pat = self.returnSampleLevel(tmp_y)
			tmp_z_pat = self.returnSampleLevel(tmp_z)

			return {'both': tmp_x,
				'left': tmp_z.union(tmp_x), 
				'right': tmp_y.union(tmp_x), 
				'cons_cid': clist[i:i+r],
				'wide_region': ( self.getRegion(i,cnatype)[0], self.getRegion(i+r,cnatype)[1] ),
				'ucount': len(tmp_x), 
				'rcount':len(tmp_y.union(tmp_x)), 
				'lcount':len(tmp_z.union(tmp_x)), 
				'ucount_pat': len(tmp_x_pat), 
				'rcount_pat': len(tmp_y_pat.union(tmp_x_pat)), 				
				'lcount_pat': len(tmp_z_pat.union(tmp_x_pat)), 				
				'p-val' : 1,
				'localc': str(delta)} 

	def returnSampleLevel(self, cliq_set):
		return set(self.edgetoPatient[x] for x in cliq_set)

	def getConsecutiveCliqueInfo(self, i, r, cnatype):
		# i: start clique, r: size of block (consequtive cliques)
		if str(i) + "-" + str(r) + "-" + cnatype in self.previous_record:
			return self.previous_record[str(i) + "-" + str(r) + "-" + cnatype]
		else:

			cscSet = set()
			in_intervals = set()
			uni_set = set()		
			upstream_set = set()
			downstream_set = set()			

			if cnatype == 'A':
				clist = self.cidOrderListA
			else:
				clist = self.cidOrderListD


			if len(clist) == 1:						
				for sid in self.idata[clist[0]]:				
					if clist[0] in self.segDegree[sid]:
						uni_set.add(sid)																
			else:		
				cscSet.add(clist[i])
				in_intervals.update(self.idata[clist[i]])
				for j in range(1, r):
					cscSet.add(clist[i+j])
					in_intervals.intersection_update(self.idata[clist[i+j]])
			
				if i == 0: # begin
					if i + r < len(clist):
						cid_r = clist[i + r]
						
						for sid in in_intervals:
							if cid_r not in self.segDegree[sid] :
								uni_set.add(sid)															
							if cid_r in self.segDegree[sid] :
								downstream_set.add(sid)														
						
				elif i + r - 1 == len(clist) - 1 :
					if i > 0:
						cid_l = clist[i - 1]
						
						for sid in in_intervals:							
							if cid_l not in self.segDegree[sid] :
								uni_set.add(sid)																					
							if cid_l in self.segDegree[sid] :
								upstream_set.add(sid)																															
				else:
					cid_l = clist[i - 1]
					cid_r = clist[i + r]
					
					for sid in in_intervals:						
						if cid_l not in self.segDegree[sid] and cid_r not in self.segDegree[sid] :
							uni_set.add(sid)														
						if cid_l in self.segDegree[sid] and cid_r not in self.segDegree[sid]:
							upstream_set.add(sid)													
						if cid_l not in self.segDegree[sid] and cid_r in self.segDegree[sid]:
							downstream_set.add(sid)														
						
			self.previous_record[str(i) + "-" + str(r) + "-" + cnatype] = list([set(uni_set), set(upstream_set), set(downstream_set)])
			return uni_set, upstream_set, downstream_set


	def getConsecutiveCliqueInfo_Fast(self, i, j, pivot, cnatype):
		# i: start clique, j: end clique, pivot, cnatype
		
		cscSet = set()
		in_intervals = set()
		uni_set = set()		
		upstream_set = set()
		downstream_set = set()			

		if cnatype == 'A':
			clist = self.cidOrderListA
		else:
			clist = self.cidOrderListD

		in_intervals = self.idata[clist[pivot]]
		
		if i == 0: # begin
			if j + 1 < len(clist):
				cid_r = clist[j+1]
				for sid in in_intervals:
					if cid_r not in self.segDegree[sid] :
						uni_set.add(sid)															
					if cid_r in self.segDegree[sid] :
						downstream_set.add(sid)														
				
		elif j == len(clist) - 1 :
			if i > 0:
				cid_l = clist[i - 1]				
				for sid in in_intervals:							
					if cid_l not in self.segDegree[sid] :
						uni_set.add(sid)																					
					if cid_l in self.segDegree[sid] :
						upstream_set.add(sid)																															
		else:
			cid_l = clist[i - 1]
			cid_r = clist[j + 1]
			
			for sid in in_intervals:						
				if cid_l not in self.segDegree[sid] and cid_r not in self.segDegree[sid] :
					uni_set.add(sid)														
				if cid_l in self.segDegree[sid] and cid_r not in self.segDegree[sid]:
					upstream_set.add(sid)													
				if cid_l not in self.segDegree[sid] and cid_r in self.segDegree[sid]:
					downstream_set.add(sid)														

		return downstream_set.union(uni_set), upstream_set.union(uni_set)	
		
