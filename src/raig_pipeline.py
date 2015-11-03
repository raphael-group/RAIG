#!/usr/bin/env python

import sys, os, glob
import json, operator, itertools, re, math, time

sys.path.append("../lib/") # point the library path for networkx and statsmodels
sys.path.append("/home/bournewu/pylib/lib/python2.7/site-packages/")
sys.path.append("/home/bournewu/pylib/lib/python2.7/site-packages/")
import networkx as nx
import statsmodels.stats.multitest as smm

from cliqueinterval import CliqueInterval
from genomeinfo import GenomeInfo
import raig_method as rm, utils as cna_utils

class raig_pipeline():
	
	delta_percent=0.05
	clique_cutoff=0.05		
	
	thramp = 0.1
	thrdel = -0.1
	#join_dist = 10
	#dist = 4
	isize = 0
	outfpath = ""
		
	geneLevel = 0
	dlcall = GenomeInfo()
	sample_mutatedGenes = dict()
	gene_mutatedSamples = dict()
	permut_id = 0
	t = 0

	ttype = ""
	log_fout = ""
	vis_fout = ""
	num_permutation = 0
	
	def __init__(self, gi, args):
		
		self.cna_a = 10
		self.cna_d = 10
		self.block_p = args.blocksize
		self.t = args.bpnumber
		#self.dist = args.remove_seg
		#self.join_dist = args.joint_distance
		self.ttype = args.cancer
		self.geneLevel = args.gene_level
		self.thramp = args.amp_cutoff
		self.thrdel = args.del_cutoff
		self.delta_percent = args.delta_p
		#self.clique_cutoff = args.clique_p		
		self.permut_id = args.indvidual_pid

		self.cna_file = args.input	
		#self.dlcall.clearSegInfo()	
		self.dlcall = gi			
		self.vis_foutpath = args.output + 'RAIGOUTPUT_' + args.cancer 
		self.num_permutation = args.num_permutation
		
	def run(self, chrm):

		# output log file					
		leftregion = 0
		rightregion = self.dlcall.p_arm_info[chrm][1]
		self.isize = self.block_p*(rightregion-leftregion)
		if self.permut_id == 0:
			#self.log_fout = open(self.vis_foutpath + "_" + str(chrm) +".log", 'w')
			#self.log_fout.write("gl:" + str(self.geneLevel) + " delta: " + str(self.cna_a) +"-"+str(self.cna_d) + " block_size: " + str(self.isize) + " bpnumber: " + str(self.t) + " thraml: " +  str(self.thramp) + " thrdel: " + str(self.thrdel) + "\n")		
			#self.log_fout.write("delta_p:" + str(self.delta_percent) + " clique_cutoff: " + str(self.clique_cutoff) + "\n")		
			#self.log_fout.write("cnaseg file: " + self.cna_file + "\n")
			#self.log_fout.write("# Chromosome " + chrm + "p... " + "(" + str(leftregion) + "\t" + str(rightregion) + ")\n")		
			self.final_fout = open(self.vis_foutpath + "_" + str(chrm) +".lst", 'w')		
			self.final_fout.write(self.formatHeaderString(self.t))
		else:
			self.final_fout = open(self.vis_foutpath + "_" + str(chrm) + "_" + str(self.permut_id) + ".lst", 'w')		
			self.final_fout.write(self.formatHeaderString(self.t))
				
		self.findmeta(leftregion, rightregion, chrm, "p")		
		leftregion = self.dlcall.q_arm_info[chrm][0]
		rightregion = self.dlcall.q_arm_info[chrm][1]
		self.isize = self.block_p*(rightregion-leftregion)
		
		#if self.permut_id == 0:
			#self.log_fout.write("gl:" + str(self.geneLevel) + " delta: " + str(self.cna_a) +"-"+str(self.cna_d) + " block_size: " + str(self.isize) + " bpnumber: " + str(self.t) + " thraml: " +  str(self.thramp) + " thrdel: " + str(self.thrdel) + "\n")		
			#self.log_fout.write("# Chromosome " + chrm + "q... " + "(" + str(leftregion) + "\t" + str(rightregion) + ")\n")
		
		self.findmeta(leftregion, rightregion, chrm, "q")

		#if self.permut_id == 0:
		#	self.log_fout.close()
		
		self.final_fout.close()
	
	def returnMatrix(self):
		return self.sample_mutatedGenes
			
	
	def formatGraph(self, edgesetA, edgesetD, posA, posD, metagenebkt):
		# 2013/03/20, gene level only applys on Deletion
		# 2015/11/03, update to faster way to construct graph.		
		G = nx.Graph()
		open_things = set()
		for pos_key in sorted(posA.keys()): 
			eidsets = posA[pos_key]				
			# put an edge if eids in the same position
			for eid1, eid2 in itertools.combinations(eidsets, 2): 
				G.add_edge(eid1, eid2)
			# do disclosure procedure belows				
			# step 1: remove eids in open_things if existing in open_things
			exist = open_things.intersection(eidsets)
			new_eids = set(eidsets) - exist
			open_things = open_things - exist
			# step 2: make edges between s_e_pos to open_things
			for eid1, eid2 in itertools.product(new_eids, open_things):
				G.add_edge(eid1, eid2)
			open_things.update(new_eids)

		if not self.geneLevel:			
			open_things = set()
			for pos_key in sorted(posD.keys()): 
				eidsets = posD[pos_key]				
				# put an edge if eids in the same position
				for eid1, eid2 in itertools.combinations(eidsets, 2): 
					G.add_edge(eid1, eid2)
				# do disclosure procedure belows				
				# step 1: remove eids in open_things if existing in open_things
				exist = open_things.intersection(eidsets)
				new_eids = set(eidsets) - exist
				open_things = open_things - exist
				# step 2: make edges between s_e_pos to open_things
				for eid1, eid2 in itertools.product(new_eids, open_things):
					G.add_edge(eid1, eid2)
				open_things.update(new_eids)
		else:						
			for eid1, eid2 in itertools.combinations(edgesetD.keys(), 2):
				if not (edgesetD[eid1][0] > edgesetD[eid2][1] or edgesetD[eid1][1] < edgesetD[eid2][0]) :
					G.add_edge(eid1, eid2)								
				else:					
					for g in metagenebkt: # if two non-overlpping segments are located in the same gene, treat them as overlapping. # 20130509 extend segment to gene boundary
						if (not (edgesetD[eid1][0] > metagenebkt[g][1] or edgesetD[eid1][1] < metagenebkt[g][0])) and (not (edgesetD[eid2][0] > metagenebkt[g][1] or edgesetD[eid2][1] < metagenebkt[g][0])):
							G.add_edge(eid1, eid2)							
										
												
		return G
	

	def RAIG_algo(self, edgesetA, edgesetD, edgetoPatient, edgewA, edgewD, posA, posD, chrm_genebkt, NumAmp, NumDel):
 		
		# format Graph for maximal clique identification				
		G = self.formatGraph(edgesetA, edgesetD, posA, posD, chrm_genebkt)
		# finding all maximal cliques		
		cliq_results = list(nx.find_cliques(G))		
		
		# put all clique information into C, including location, cna type, intervals belonging to the clique
		C = CliqueInterval()
		C.edgetoPatient = dict(edgetoPatient)		
		#C.createBiGraph(cliq_results, round(NumAmp * self.clique_cutoff), edgesetA, round(NumDel * self.clique_cutoff), edgesetD)		
		C.createBiGraph(cliq_results, NumAmp * self.clique_cutoff, edgesetA, NumDel * self.clique_cutoff, edgesetD)		
		
		unusedA = dict()
		unusedD = dict()
		cscPairA = dict()
		cscPairD = dict()

		# start RAIG Algo
		cscPairA = rm.dynamicProg_RAIG(C, self.isize, max(NumAmp*self.delta_percent, self.cna_a), 'A', edgewA)
		cscPairD = rm.dynamicProg_RAIG(C, self.isize, max(NumDel*self.delta_percent, self.cna_d), 'D', edgewD)									
		
		return cscPairA, cscPairD

	def significance_assessment(self, cscPairA, cscPairD, leftregion, rightregion, meta_chrome, arm, AmpPat, DelPat, chrm_genebkt):
		if len(cscPairA.keys()) != 0 or len(cscPairD.keys()) != 0:
			scorelistA, scorelistD = [], []
			for i in range(0, self.num_permutation):				
				permute_regionA, permute_regionD = cna_utils.cycle_shift_permutation(self.dlcall.regionA[meta_chrome][arm], self.dlcall.regionD[meta_chrome][arm], leftregion, rightregion)				
				pedgesetA, pedgesetD, pedgetoPatient, pedgewA, pedgewD, pposA, pposD = cna_utils.formatEdgeId(AmpPat.union(DelPat), permute_regionA, permute_regionD)#, abbA, abbD) 
				pcscPairA, pcscPairD = self.RAIG_algo(pedgesetA, pedgesetD, pedgetoPatient, pedgewA, pedgewD, pposA, pposD, chrm_genebkt, len(AmpPat), len(DelPat))
				if len(pcscPairA.keys()) != 0:
					scorelistA.append(max([2*min(pcscPairA[cid]['lcount'],pcscPairA[cid]['rcount']) for cid in pcscPairA.keys()]))
				else:
					scorelistA.append(0)

				if len(pcscPairD.keys()) != 0:
					scorelistD.append(max([2*min(pcscPairD[cid]['lcount'],pcscPairD[cid]['rcount']) for cid in pcscPairD.keys()]))
				else:
					scorelistD.append(0)

		if len(cscPairA.keys()) != 0:
			pvals = list()
			cidlist = list()
			for cid in cscPairA.keys():
				csc_score = 2*min(cscPairA[cid]['lcount'],cscPairA[cid]['rcount'])
				count = 0
				for s in scorelistA:
					if s > csc_score:
						count += 1
				
				cscPairA[cid]['p-val'] = float(count)/self.num_permutation
				pvals.append(float(count)/self.num_permutation)
				cidlist.append(cid)
			
			corrected_pval = smm.multipletests(pvals, alpha=0.05, method='fdr_bh')[1]
			for i in range(len(cidlist)):
				cscPairA[cidlist[i]]['corrected-p-val'] = corrected_pval[i]
		
		if len(cscPairD.keys()) != 0:
			pvals = list()
			cidlist = list()
			for cid in cscPairD.keys():
				csc_score = 2*min(cscPairD[cid]['lcount'],cscPairD[cid]['rcount'])
				count = 0
				for s in scorelistD:
					if s > csc_score:
						count +=1
				
				cscPairD[cid]['p-val'] = float(count)/self.num_permutation
				pvals.append(float(count)/self.num_permutation)
				cidlist.append(cid)
			
			corrected_pval = smm.multipletests(pvals, alpha=0.05, method='fdr_bh')[1]
			for i in range(len(cidlist)):
				cscPairD[cidlist[i]]['corrected-p-val'] = corrected_pval[i]

	def findmeta(self, leftregion, rightregion, meta_chrome, arm):	

		if self.permut_id != 0:
			
			chrm_genebkt = self.dlcall.gene_container[meta_chrome][arm]
			AmpPat = set(filter(lambda x:len(self.dlcall.regionA[meta_chrome][arm][x])!=0, self.dlcall.regionA[meta_chrome][arm].keys())) if self.dlcall.regionA else set()
			DelPat = set(filter(lambda x:len(self.dlcall.regionD[meta_chrome][arm][x])!=0, self.dlcall.regionD[meta_chrome][arm].keys())) if self.dlcall.regionD else set()
			
			permute_regionA, permute_regionD = cna_utils.cycle_shift_permutation(self.dlcall.regionA[meta_chrome][arm], self.dlcall.regionD[meta_chrome][arm], leftregion, rightregion)		
			pedgesetA, pedgesetD, pedgetoPatient, pedgewA, pedgewD, pposA, pposD = cna_utils.formatEdgeId(AmpPat.union(DelPat), permute_regionA, permute_regionD)#, abbA, abbD) 
			pcscPairA, pcscPairD = self.RAIG_algo(pedgesetA, pedgesetD, pedgetoPatient, pedgewA, pedgewD, pposA, pposD, chrm_genebkt, len(AmpPat), len(DelPat))
			
			for cid in sorted(pcscPairA.keys(), key=lambda x: 2*min(pcscPairA[x]['lcount'],pcscPairA[x]['rcount']), reverse=True):
				self.final_fout.write(str(self.permut_id) + "\t" + self.formatIndividualPermutation(meta_chrome, arm, 'A', pcscPairA[cid]))
			for cid in sorted(pcscPairD.keys(), key=lambda x: 2*min(pcscPairD[x]['lcount'],pcscPairD[x]['rcount']), reverse=True):
				self.final_fout.write(str(self.permut_id) + "\t" + self.formatIndividualPermutation(meta_chrome, arm, 'D', pcscPairD[cid]))
			
			return 0			
			
		else:		
			# retrieve genes located in the region		
			# abbA, abbD = {}, {}
			chrm_genebkt = self.dlcall.gene_container[meta_chrome][arm]
			AmpPat = set(filter(lambda x:len(self.dlcall.regionA[meta_chrome][arm][x])!=0, self.dlcall.regionA[meta_chrome][arm].keys())) if self.dlcall.regionA else set()
			DelPat = set(filter(lambda x:len(self.dlcall.regionD[meta_chrome][arm][x])!=0, self.dlcall.regionD[meta_chrome][arm].keys())) if self.dlcall.regionD else set()
			
			#self.log_fout.write("Amp, Del samples: " + str(len(AmpPat)) + " " + str(len(DelPat)) + "\n")
			#self.log_fout.write("Amp, Del delta: " + str(len(AmpPat)*self.delta_percent) + " " + str(len(DelPat)*self.delta_percent) + "\n")
			
			edgesetA, edgesetD, edgetoPatient, edgewA, edgewD, posA, posD = cna_utils.formatEdgeId(AmpPat.union(DelPat), self.dlcall.regionA[meta_chrome][arm], self.dlcall.regionD[meta_chrome][arm])#, abbA, abbD) 				
			cscPairA, cscPairD = self.RAIG_algo(edgesetA, edgesetD, edgetoPatient, edgewA, edgewD, posA, posD, chrm_genebkt, len(AmpPat), len(DelPat)	)
				
			# assess significance
			if self.num_permutation > 0:
				self.significance_assessment(cscPairA, cscPairD, leftregion, rightregion, meta_chrome, arm, AmpPat, DelPat, chrm_genebkt)
				
			
			# OUTPUT
			temp_outA = dict()
			temp_outD = dict()		
			
			for cid in set(cscPairA.keys()).union(set(cscPairD.keys()))	:  # for each target region
				
				DorA = 'A' if cid in cscPairA else 'D'
				sampleCoverage=set()
				leftmost_t = 0 
				rightmost_t = 0		

				if DorA == "A":							
					cscPairA[cid]['region'] = dict()
					cscPairA[cid]['genes'] = dict()
					for u in range(self.t + 1):
						tmpmeta, leftmost_t, rightmost_t = self.getTargetRegionfromEndPoint(cscPairA[cid], chrm_genebkt, edgesetA, edgetoPatient, u)
						cscPairA[cid]['region'][u] = (leftmost_t, rightmost_t)						
						cscPairA[cid]['genes'][u] = tmpmeta
						
					tmp_key = ",".join(cscPairA[cid]['genes'][self.t]) if cscPairA[cid]['genes'][self.t] else ",".join(cscPairA[cid]['genes'][0])
					
					if tmp_key in temp_outA:
						tmp_key_1 = tmp_key + '-cid-' + str(cid)
						temp_outA[tmp_key_1] = dict()
						temp_outA[tmp_key_1]['output'] = self.formatLogString('A', cid, cscPairA[cid], tmp_key, 0, sampleCoverage)
						temp_outA[tmp_key_1]['print'] = self.formatOutputString(meta_chrome, arm, 'A', cscPairA[cid], tmp_key)
						#temp_outA[tmp_key_1]['cliqnum'] = len(cscPairA[cid]['cons_cid'])
						#temp_outA[tmp_key_1]['uniqnum'] =  cscPairA[cid]['ucount'] 
						temp_outA[tmp_key_1]['cid'] = cid
						temp_outA[tmp_key_1]['pval'] = cscPairA[cid]['p-val']
						
					else:					
						temp_outA[tmp_key] = dict()
						temp_outA[tmp_key]['output'] = self.formatLogString('A', cid, cscPairA[cid], tmp_key, 0, sampleCoverage)
						temp_outA[tmp_key]['print'] = self.formatOutputString(meta_chrome, arm, 'A', cscPairA[cid], tmp_key)
						#temp_outA[tmp_key]['cliqnum'] = len(cscPairA[cid]['cons_cid'])
						#temp_outA[tmp_key]['uniqnum'] =  cscPairA[cid]['ucount'] 
						temp_outA[tmp_key]['cid'] = cid
						temp_outA[tmp_key]['pval'] = cscPairA[cid]['p-val']
				
				if DorA == "D":								
					if self.geneLevel == 1:
						cscPairD[cid]['region'] = dict()
						cscPairD[cid]['genes'] = dict()
						for u in range(self.t + 1):
							tmpmeta, leftmost_t, rightmost_t = self.getTargetRegionfromEndPointAuto(cscPairD[cid], chrm_genebkt, edgesetD, edgetoPatient, u)
							cscPairD[cid]['region'][u] = (leftmost_t, rightmost_t)						
							cscPairD[cid]['genes'][u] = tmpmeta
					else:
						cscPairD[cid]['region'] = dict()
						cscPairD[cid]['genes'] = dict()
						for u in range(self.t + 1):
							tmpmeta, leftmost_t, rightmost_t = self.getTargetRegionfromEndPoint(cscPairD[cid], chrm_genebkt, edgesetD, edgetoPatient, u)
							cscPairD[cid]['region'][u] = (leftmost_t, rightmost_t)						
							cscPairD[cid]['genes'][u] = tmpmeta

					tmp_key = ",".join(cscPairD[cid]['genes'][self.t]) if cscPairD[cid]['genes'][self.t] else ",".join(cscPairD[cid]['genes'][0])

					if tmp_key in temp_outD:
						tmp_key_1 = tmp_key + '-cid-' + str(cid)
						temp_outD[tmp_key_1] = dict()
						temp_outD[tmp_key_1]['output'] = self.formatLogString('D', cid, cscPairD[cid], tmp_key, 0, sampleCoverage)
						temp_outD[tmp_key_1]['print'] = self.formatOutputString(meta_chrome, arm, 'D', cscPairD[cid], tmp_key)
						#temp_outD[tmp_key_1]['cliqnum'] = len(cscPairD[cid]['cons_cid'])
						#temp_outD[tmp_key_1]['uniqnum'] = cscPairD[cid]['ucount']
						temp_outD[tmp_key_1]['cid'] = cid		
						temp_outD[tmp_key_1]['pval'] = cscPairD[cid]['p-val']
						
					else:					
						temp_outD[tmp_key] = dict()
						temp_outD[tmp_key]['output'] = self.formatLogString('D', cid, cscPairD[cid], tmp_key, 0, sampleCoverage)
						temp_outD[tmp_key]['print'] = self.formatOutputString(meta_chrome, arm, 'D', cscPairD[cid], tmp_key)
						#temp_outD[tmp_key]['cliqnum'] = len(cscPairD[cid]['cons_cid'])
						#temp_outD[tmp_key]['uniqnum'] = cscPairD[cid]['ucount']
						temp_outD[tmp_key]['cid'] = cid		
						temp_outD[tmp_key]['pval'] = cscPairD[cid]['p-val']

			for tmp_key in sorted(temp_outA.keys(), key=lambda x:temp_outA[x]['pval']):				
				#self.log_fout.write(temp_outA[tmp_key]['output'])
				self.final_fout.write(temp_outA[tmp_key]['print'])
				
			for tmp_key in sorted(temp_outD.keys(), key=lambda x:temp_outD[x]['pval']):				
				#self.log_fout.write(temp_outD[tmp_key]['output'])
				self.final_fout.write(temp_outD[tmp_key]['print'])
				
			
	# end of findmeta()

	def formatLogString(self, dora, cid, indict, wgenename, ingene_seg, coverage):		
		region_string = "\t".join([str(indict['region'][u][0])+","+str(indict['region'][u][1]) for u in sorted(indict['region'].keys())])
		gene_string = "\t".join([",".join(indict['genes'][u]) for u in sorted(indict['genes'].keys())])
		if self.num_permutation > 0:
			return "TARGET REGION (" + dora + "):" + str(cid) + "\t" + ",".join(str(x) for x in indict['cons_cid']) + "\t" + wgenename + "\tCNA:" + str(len(indict['both'])) + "," +str(len(indict['left'].difference(indict['both'])))+ "," +str(len(indict['right'].difference(indict['both'])))+ ","+str(ingene_seg)+"\t"+ str(len(coverage)) + "\t" + region_string + "\t" + str(indict['lcount']) +"," + str(indict['rcount']) +"\t" + indict['localc'] + "\t" + str(indict['p-val']) +"\t" + str(indict['corrected-p-val']) + "\n"				
		else:
			return "TARGET REGION (" + dora + "):" + str(cid) + "\t" + ",".join(str(x) for x in indict['cons_cid']) + "\t" + wgenename + "\tCNA:" + str(len(indict['both'])) + "," +str(len(indict['left'].difference(indict['both'])))+ "," +str(len(indict['right'].difference(indict['both'])))+ ","+str(ingene_seg)+"\t"+ str(len(coverage)) + "\t" + region_string + "\t" + str(indict['lcount']) +"," + str(indict['rcount']) +"\t" + indict['localc'] + "\t1.0\t1.0\n"				

	def formatHeaderString(self, t):
		region_string = "\t".join([ 'region_start_t='+str(u) + "\t" + 'region_end_t='+str(u) for u in range(t+1)])
		gene_string = "\t".join([ 'gene_t='+str(u) for u in range(t+1)])
		return "\t".join(['chrm+arm', 'cnatype', region_string, gene_string, 'W(B)', 'p-val', 'corrected-p-val']) + "\n"

	def formatOutputString(self, chrm, arm, dora, indict, wgenename):		
		region_string = "\t".join([str(indict['region'][u][0])+"\t"+str(indict['region'][u][1]) for u in sorted(indict['region'].keys())])
		gene_string = "\t".join([",".join(indict['genes'][u]) for u in sorted(indict['genes'].keys())])
		if self.num_permutation > 0:			
			return "\t".join([str(chrm)+arm, dora, region_string, gene_string, str(2*min(indict['lcount'], indict['rcount'])), str(indict['p-val']), str(indict['corrected-p-val'])]) + "\n"
			
		else:			
			return "\t".join([str(chrm)+arm, dora, region_string, gene_string, str(2*min(indict['lcount'], indict['rcount'])), "-", "-"]) + "\n"

	def formatIndividualPermutation(self, chrm, arm, dora, indict):		
			return "\t".join([str(chrm)+arm, dora, str(indict['wide_region'][0]), str(indict['wide_region'][1]), str(2*min(indict['lcount'], indict['rcount']))]) + "\n"		

	def getTargetRegionfromEndPoint(self, indict, metagenebkt, edgeset, edgetoSample, localT):	

		uniqGenes = set()
		leftpoint = set()
		rightpoint = set()		
		# create left and right end point from each left interval and right interval
		for eid in indict['left']:
			leftpoint.add(edgeset[int(eid)][0])
		for eid in indict['right']:
			rightpoint.add(edgeset[int(eid)][1])

		s_leftp = sorted(leftpoint, reverse=True)
		s_rightp = sorted(rightpoint)

		# get the intersect region					
		leftbound, rightbound = self.getIntersectRegionPairing(s_leftp, s_rightp, localT)
		tmpbound = 0
		if leftbound > rightbound:
			tmpbound = leftbound
			leftbound = rightbound
			rightbound = tmpbound
		
		# collect genes in the region
		uniqGenes = cna_utils.collectGenes(metagenebkt, leftbound, rightbound)
										
		if not uniqGenes: # no genes in the clique, output empty set and real clique region
			#tmp1, tmp2 = self.getCliqueRange(indict['left'].union(indict['right']), edgeset)
			#uniqGenes[str(tmp1)] = 0
			#uniqGenes[str(tmp2)] = 0			
			return self.findCloseGenes(metagenebkt, leftbound, rightbound), leftbound, rightbound
		else:
			#sorted_meta = sorted(uniqGenes.keys(), key=lambda g: metagenebkt[g][0])
			return sorted(uniqGenes, key=lambda g: metagenebkt[g][0]), leftbound, rightbound



	def findCloseGenes(self, metagenebkt, left, right):
		minleft = 10**10
		minright = 10**10
		leftg = ""
		rightg= ""
		for g in metagenebkt:
			if metagenebkt[g][1] < left and left - metagenebkt[g][1] < minleft:
				leftg = g
				minleft = left - metagenebkt[g][1]

			if metagenebkt[g][0] > right and metagenebkt[g][0] - right < minright:
				rightg = g
				minright = metagenebkt[g][0] - right

		#if minright < minleft:
		#	return {"["+rightg+"]":set()}
		#else:
		#	return {"["+leftg+"]":set()}
		return set(["["+leftg+","+rightg+"]"])

	def getTargetRegionfromEndPointAuto(self, indict, metagenebkt, edgeset, edgetoSample, localT):	
		# for Gene searching function ON, start counting region outside boudaries
		uniqGenes = set()
		tmpt = 0
		startt = 0 - localT - 1
		leftpoint = set()
		rightpoint = set()		
		# create left and right end point from each left interval and right interval
		for eid in indict['left']:
			leftpoint.add(edgeset[int(eid)][0])
		for eid in indict['right']:
			rightpoint.add(edgeset[int(eid)][1])

		s_leftp = sorted(leftpoint, reverse=True)
		s_rightp = sorted(rightpoint)

		cut_left = 0
		for i in range(len(s_leftp)):
			if all( bp_r >= s_leftp[i] for bp_r in s_rightp ):
				cut_left = i
				break

		cut_right = 0
		for i in range(len(s_rightp)):
			if all(bp_l <= s_rightp[i] for bp_l in s_leftp):
				cut_right = i
				break


		leftbound, rightbound = self.getIntersectRegionPairing(s_leftp[cut_left:], s_rightp[cut_right:], localT)

		tmpbound = 0
		if leftbound > rightbound:
			tmpbound = leftbound
			leftbound = rightbound
			rightbound = tmpbound

		uniqGenes = cna_utils.collectGenes(metagenebkt, leftbound, rightbound)
		"""
		while 1:						
			# get the intersect region					
			leftbound, rightbound = self.getIntersectRegion(s_leftp, s_rightp, tmpt)

			# collect genes in the region
			uniqGenes = cna_utils.collectGenes(metagenebkt, leftbound, rightbound)		
			
			if leftbound == s_leftp[-1] and rightbound == s_rightp[-1]:
				break
			if uniqGenes and startt== 0 - localT - 1: # record the initial t (startt) such that any genes start including
				startt = tmpt
			if uniqGenes and startt + localT == tmpt: # if reach T, stop
				break			
			
			tmpt += 1
			uniqGenes.clear()
		"""					
		if not uniqGenes: # no genes in the clique, output empty set and real clique region			
			return self.findCloseGenes(metagenebkt, leftbound, rightbound), leftbound, rightbound
		else:
			sorted_meta = sorted(uniqGenes, key=lambda g: metagenebkt[g][0])
			#return sorted(uniqGenes, key=lambda g: metagenebkt[g][0]), leftbound, rightbound
			return sorted(uniqGenes, key=lambda g: metagenebkt[g][0]), metagenebkt[sorted_meta[0]][0], metagenebkt[sorted_meta[-1]][1]

	def getIntersectRegion(self, sort_left, sort_right, tmp_t):
		leftbound = 0
		rightbound = 0

		if len(sort_left) + len(sort_right) -2 <= tmp_t:
			leftbound = sort_left[-1]
			rightbound = sort_right[-1]
		else:
			minrange = abs(sort_left[-1] - sort_right[-1])
			for li in range(tmp_t + 1):
				ri = (tmp_t) - li
				if li < len(sort_left) and ri < len(sort_right):
					tmpl = sort_left[li]
					tmpr = sort_right[ri]
					if abs(tmpl - tmpr) < minrange:
						leftbound = tmpl
						rightbound = tmpr

		return leftbound, rightbound

	def getIntersectRegionPairing(self, sort_left, sort_right, tmp_t):
		
		return sort_left[min(tmp_t, len(sort_left)-1)], sort_right[min(tmp_t, len(sort_right)-1)]

	

	
