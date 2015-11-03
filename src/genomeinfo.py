'''
Created on Jul 13, 2012

@author: bournewu
'''

class GenomeInfo():
    '''
    classdocs
    '''
    gene_container = dict()  
    chrmGene = dict() 
    chromregion = dict() 
    regionA = dict() 
    regionD = dict() 
    bregionA = dict() 
    bregionD = dict() 
    weightA = dict()
    weightD = dict()
    broad_cutoff = dict()
    p_arm_info = dict()
    q_arm_info = dict()
    sampleSet = set()
    SNVgeneToSamples = dict()
    probeInfo = dict()
    probeCheck = dict()
    germline_cnv = dict()
    germline_cnv_out = dict()

    marker_file = ""
    
    def __init__(self):
        '''
        Constructor
        '''
        
    def clearSegInfo(self):
        self.regionA = dict()
        self.regionD = dict()
        self.bregionA = dict()
        self.bregionD = dict()
        self.chromregion = dict()

    def loadGermlineCNV(self):
        for line in open(self.germline_cnv_file, 'r'):
            v = line.rstrip().split("\t")
            if v[1] not in self.germline_cnv:
                self.germline_cnv[v[1]] = list()
                self.germline_cnv_out[v[1]] = list()
		
            self.germline_cnv[v[1]].append([int(v[4]), int(v[5])])
            self.germline_cnv_out[v[1]].append([int(v[2]), int(v[3])])

    def passHalfArm(self, left, right, chrm):
        if right <= self.p_arm_info[chrm][1]:
            if right - left + 1 < self.p_arm_info[chrm][2]: # smaller than x% of p arm
                return True, 'p'
            else:
                return False, '-'
        elif left >= self.q_arm_info[chrm][0]:
            if right - left +1 < self.q_arm_info[chrm][2]: # smaller than x% of q arm
                return True, 'q'
            else:
                return False,'-'
        elif left < self.p_arm_info[chrm][1] and right > self.p_arm_info[chrm][1]:
            if right < self.q_arm_info[chrm][0]: # part of segments cross centromere
                if right - left + 1 < self.p_arm_info[chrm][2]: # smaller than x% of p arm
                    return True, 'p'
                else: 
                    return False, '-'
            else:
                # whole segments cross centromere
                return False, '-'
        elif right > self.q_arm_info[chrm][0] and left < self.q_arm_info[chrm][0]: #fixed 20140403 p_arm_info[chrm][0] -> q_arm_info[chrm][0]
            if left > self.p_arm_info[chrm][1]: # part of segments cross centromere
                if right - left + 1 < self.q_arm_info[chrm][2]: # smaller than x% of q arm
                    return True, 'q'
                else:
                    return False, '-'
            else:
                # whole segments cross centromere
                False, '-'        
        else:
            return False, '-'
    
    def readArmInfo(self, arm_json, broad_ratio, pchrm_set):
        import json
        with open(arm_json, 'r') as json_file:
            data= json.load(json_file)
        
        for chrm in data:
            if chrm in pchrm_set:
                # don't consider centromere
                #self.p_arm_info[chrm] = [data[chrm]['p'][0]['start'], data[chrm]['p'][-2]['end'], int(float(data[chrm]['p'][-1]['end'] - data[chrm]['p'][0]['start'])/2)]
                #self.q_arm_info[chrm] = [data[chrm]['q'][1]['start'], data[chrm]['q'][-1]['end'], int(float(data[chrm]['q'][-1]['end'] - data[chrm]['q'][0]['start'])/2)]
                self.p_arm_info[chrm] = [data[chrm]['p'][0]['start'], data[chrm]['p'][-1]['end'], int(float(data[chrm]['p'][-1]['end'] - data[chrm]['p'][0]['start'])*broad_ratio)]
                self.q_arm_info[chrm] = [data[chrm]['q'][0]['start'], data[chrm]['q'][-1]['end'], int(float(data[chrm]['q'][-1]['end'] - data[chrm]['q'][0]['start'])*broad_ratio)]
            
    def readGeneInfo(self, gene_json, pchrm_set):        
        import json
        with open(gene_json, 'r') as json_file:
            data = json.load(json_file)

        for chrm in pchrm_set:
            self.gene_container[chrm] = dict()
            self.gene_container[chrm]['p'] = dict()
            self.gene_container[chrm]['q'] = dict()

            for g, g_range in data[chrm].items():
                if not( g_range[1] < self.p_arm_info[chrm][0] or g_range[0] > self.p_arm_info[chrm][1] ):
                    self.gene_container[chrm]['p'][g] = g_range

                if not( g_range[1] < self.q_arm_info[chrm][0] or g_range[0] > self.q_arm_info[chrm][1] ):
                    self.gene_container[chrm]['q'][g] = g_range

	
    def readSegmentInfo(self, fname, thramp, thrdel, pchrm_set, cond_join_seg):                
        
        for line in open(fname, 'r'):
            line = line.strip("\n")
            v = line.split("\t")
                        
            if v[1] in pchrm_set and v[0] != 'Sample' :
                passHalfTag, arm = self.passHalfArm(int(v[2]), int(v[3]), v[1])

                if v[1] not in self.regionA:
                    self.regionA[v[1]] = dict()
                    self.regionA[v[1]]['p'] = dict()
                    self.regionA[v[1]]['q'] = dict()
                    self.regionD[v[1]] = dict()
                    self.regionD[v[1]]['p'] = dict()
                    self.regionD[v[1]]['q'] = dict()


                if int(v[4]) >= cond_join_seg and passHalfTag:
                    if thramp > 0 and float(v[-1]) >= thramp :#and float(v[-1]) <= capamp:
                        newpat = "-".join(v[0].split('-')[:3])                        
                        self.sampleSet.add(newpat)
                        
                        if newpat not in self.regionA[v[1]][arm]:
                            self.regionA[v[1]][arm][newpat] = list()
                
                        self.regionA[v[1]][arm][newpat].append(int(v[2]))
                        self.regionA[v[1]][arm][newpat].append(int(v[3]))                                            
        
                    if thrdel < 0 and float(v[-1]) <= thrdel :#and float(v[-1]) <= capdel:                        
                        newpat = "-".join(v[0].split('-')[:3])
                        self.sampleSet.add(newpat)
         
                        if newpat not in self.regionD[v[1]][arm]:
                            self.regionD[v[1]][arm][newpat] = list()
                    
                        self.regionD[v[1]][arm][newpat].append(int(v[2]))
                        self.regionD[v[1]][arm][newpat].append(int(v[3]))                        
    
                    if v[1] not in self.chromregion:
                        self.chromregion[v[1]] = [int(v[2]), int(v[3])]
                    else:
                        if int(v[2]) < self.chromregion[v[1]][0]:
                            self.chromregion[v[1]][0] = int(v[2])
                        if int(v[3]) > self.chromregion[v[1]][1]:
                            self.chromregion[v[1]][1] = int(v[3])
                    



    def removeCloseSegments(outregion, meta_chrome, probeInfo, probeCheck, join_dist):
        """
        Remove segments if the number of probes with in the segment is smaller than the cutoff, default: 4
        """
        for ppp in outregion:
            if len(outregion[ppp]) > 0:
                
                i = 0
                while True:
                    if outregion[ppp][i] not in probeCheck[meta_chrome] or outregion[ppp][i + 1] not in probeCheck[meta_chrome]:
                        i += 2
                    else:
                        if probeInfo[meta_chrome].index(outregion[ppp][i + 1]) - probeInfo[meta_chrome].index(outregion[ppp][i]) + 1 < join_dist:
                            rmv1 = outregion[ppp][i]
                            rmv2 = outregion[ppp][i + 1]
                            outregion[ppp].remove(rmv1)
                            outregion[ppp].remove(rmv2)
                        else:
                            i += 2          
                    if i == len(outregion[ppp]) or len(outregion[ppp]) == 0:
                        break

    def mergeCloseSegments(self, pchrm_set, dist):
        """
        merge segments if the number of probes between the segments is smaller than the cutoff, default: 15
        """
        for meta_chrome in pchrm_set:
            if self.regionA:
                for arm, outregion in self.regionA[meta_chrome].items():
                    for ppp in outregion:
                        if len(outregion[ppp]) > 2:
                            i = 0
                            while True:
                                if outregion[ppp][i + 2] not in self.probeCheck[meta_chrome] or outregion[ppp][i + 1] not in self.probeCheck[meta_chrome]:
                                    i += 2
                                else:
                                    if self.probeInfo[meta_chrome].index(outregion[ppp][i + 2]) - self.probeInfo[meta_chrome].index(outregion[ppp][i + 1]) < dist:
                                        rmv1 = outregion[ppp][i + 1]
                                        rmv2 = outregion[ppp][i + 2]
                                        outregion[ppp].remove(rmv1)
                                        outregion[ppp].remove(rmv2)
                                        
                                    else:
                                        i += 2          
                                if i + 2 == len(outregion[ppp]) or len(outregion[ppp]) == 2:
                                    break

            if self.regionD:
                for arm, outregion in self.regionD[meta_chrome].items():
                    for ppp in outregion:
                        if len(outregion[ppp]) > 2:
                            i = 0
                            while True:
                                if outregion[ppp][i + 2] not in self.probeCheck[meta_chrome] or outregion[ppp][i + 1] not in self.probeCheck[meta_chrome]:
                                    i += 2
                                else:
                                    if self.probeInfo[meta_chrome].index(outregion[ppp][i + 2]) - self.probeInfo[meta_chrome].index(outregion[ppp][i + 1]) < dist:
                                        rmv1 = outregion[ppp][i + 1]
                                        rmv2 = outregion[ppp][i + 2]
                                        outregion[ppp].remove(rmv1)
                                        outregion[ppp].remove(rmv2)                                    

                                    else:
                                        i += 2          
                                if i + 2 == len(outregion[ppp]) or len(outregion[ppp]) == 2:
                                    break



    def refineSeg(self, pchrm, pleft, pright):
        newseg = list()
        for i in range(len(self.germline_cnv[pchrm])):
            
            if pleft < self.germline_cnv[pchrm][i][0] and pright > self.germline_cnv[pchrm][i][1]:
                newseg.append([pleft, self.germline_cnv_out[pchrm][i][0]])
                newseg.append([self.germline_cnv_out[pchrm][i][1], pright])
            elif pleft < self.germline_cnv[pchrm][i][0] and pright < self.germline_cnv[pchrm][i][1] and pright > self.germline_cnv[pchrm][i][0]:
                newseg.append([pleft, self.germline_cnv_out[pchrm][i][0]])
            elif pleft > self.germline_cnv[pchrm][i][0] and pright > self.germline_cnv[pchrm][i][1] and pleft < self.germline_cnv[pchrm][i][1]:
                newseg.append([self.germline_cnv_out[pchrm][i][1], pright])

        return newseg

    def readProbeInfo(self):
        fin = open(self.marker_file, 'r')
        for line in fin:
            line = line.strip("\n")
            v = line.split("\t")
            if v[1] == 'X' or v[1] == 'x':
				v[1] = str(23)
            elif v[1] == 'Y' or v[1] == 'y':
                v[1] = str(24)
    
            if v[1] not in self.probeInfo:
                self.probeInfo[v[1]] = list()
                self.probeCheck[v[1]] = set()
            self.probeInfo[v[1]].append(int(v[2]))
            self.probeCheck[v[1]].add(int(v[2]))
        fin.close()
        
        for chrm in self.probeInfo:
            self.probeInfo[chrm].sort()        
