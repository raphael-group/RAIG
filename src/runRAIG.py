#!/usr/bin/python

################################################################################
# Load required modules
import sys, os, time, raig_pipeline as rp
from genomeinfo import GenomeInfo

################################################################################
# Parse args

def parse_args(input_list=None):
    # Parse args
    import argparse
    class Args: pass
    args = Args()
    description = 'Run Clique-based approach to identify recurrent copy number aberrations.'
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-dp', '--delta_p', default=0.05, type=float, help='percentage of samples of the chromosome arm for delta cutoff')
    parser.add_argument('-b', '--blocksize', default=0.01, type=float,
                        help='size of the block to be considered. percent of the length of chromosome arm.')
    parser.add_argument('-t', '--bpnumber', default=2, type=int,
                        help='Number of boundary points to be considered when determing target regions.')
    parser.add_argument('-g', '--gene_level', default=0, type=int, help='Turn on/off the gene level function.')
    parser.add_argument('-brd', '--broad_cutoff', default=0.7, type=float,
                        help='Proportion for considering broad seg.')
    parser.add_argument('-ta', '--amp_cutoff', default=0.1, type=float,
                        help='Amplification cutoff.')
    parser.add_argument('-td', '--del_cutoff', default=-0.1, type=float,
                        help='Deletion cutoff.')    
    parser.add_argument('-p', '--num_permutation', default=0, type=int,
                        help='Number of permutation')
    parser.add_argument('-pid', '--indvidual_pid', default=0, type=int,
                        help='ID for performing permutation parallelly.')    
    parser.add_argument('-s', '--remove_seg', default=0, type=int, help='Minimun probes in the seg.')
    
    parser.add_argument('-n', '--run_chrom', nargs='+', required=True,
                        help='Chromosomes to test.')
    parser.add_argument('-o', '--output',  required=True, help='Output dir.')                    
    parser.add_argument('-i', '--input',  required=True, help='Input segmentation file.')                    
    parser.add_argument('-c', '--cancer',  required=True, help='Cancer name.')                    
    
    if input_list:
        parser.parse_args(input_list, namespace=args)
    else:
        parser.parse_args(namespace=args)

    return args

def run(args):
        
    if os.path.isfile(args.input):
        copyNumberProcessing(args)
    else:
        sys.stderr.write("WARNINGS: File " + args.input + " does not exist! \n") 
    
def copyNumberProcessing(args):    

    # read genome information
    sys.stderr.write("Data loading ...\n")
    gi = GenomeInfo()
    gi.readArmInfo('../required_files/hg19_cytoband.json', args.broad_cutoff, args.run_chrom)
    gi.readSegmentInfo(args.input, args.amp_cutoff, args.del_cutoff, args.run_chrom, args.remove_seg)        
    gi.readGeneInfo('../required_files/hg19_genes.json', args.run_chrom)        
    #gi.readProbeInfo()
    #gi.mergeCloseSegments(args.run_chrom, args.joint_distance)    
    
    # start CNA processing, setup all parameters
    start_time = time.time()
    cna = rp.raig_pipeline(gi, args)
    sys.stderr.write("RAIG starts running ...\n")
            
    for chrm in args.run_chrom:
        sys.stderr.write( "Chromosome " + str(chrm) + "\n")
        cna.run(chrm)  
    print time.time() - start_time          
    
    del cna

if __name__ == "__main__": run( parse_args() )
