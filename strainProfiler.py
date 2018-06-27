#!/usr/bin/env python

import argparse
import pysam
import sys
import os

import numpy as np
import pandas as pd

from Bio import SeqIO
from subprocess import call
from collections import defaultdict
from tqdm import tqdm

__author__ = "Matt Olm"
__version__ = "0.1.0"
__license__ = "MIT"

class Controller():
    '''
    Main controller of program
    '''

    def parseArguments(self, args):
        '''
        Parse user options and call the correct pipeline
        '''
        # Set up .bam file
        bam = self.prepare_bam_file(args)

        # get all the information possible on all scaffolds
        Bdb, Sdb = profile_bam(bam, **args)

        # Parse Sdb a little bit
        Sdb = _parse_Sdb(Sdb)

        # Add the .stb information if requested
        if args.get('stb', None) != None:
            print("stb files are not supported yet- bug Matt about this if you want this feature!")
            pass

        # output
        self.write_output(Bdb, Sdb, args)

    def prepare_bam_file(self, args):
        '''
        From the input, make a sorted .bam file
        '''
        # Set up .bam file
        if args.get('s') != None:
            if args.get('b') != None:
                print('Choose one or the other with the .sam or .bam! Not both')
                sys.exit()

            sam = args.get('s')
            bam = _sam_to_bam(sam)
            bam = _sort_index_bam(bam)

        else:
            bam = args.get('b')
            if os.path.exists(bam + '.bai'):
                pass
            else:
                bam = _sort_index_bam(bam, rm_ori=False)

        if os.stat(bam).st_size == 0:
            print("Failed to generated a sorted .bam file! Make sure you have "+\
                "samtools version 1.6 or greater.")
            sys.exit()
        return bam

    def write_output(self, Bdb, Sdb, args):
        '''
        Write output files
        '''
        out_base = args.get('o')

        # Write a log file
        with open(out_base + '_log', 'w') as o:
            o.write("strainProfiler verion {0}\n".format(__version__))
            for a in args.keys():
                o.write("{0}\t{1}\n".format(a, args[a]))

        # Write the scaffold-level profile
        Bdb.to_csv(out_base + '_scaffoldTable.csv', index=False)

        # Write the SNPs if they exist
        if len(Sdb) == 0:
            pass
        else:
            Sdb.to_pickle(out_base + '_snpLocations.pickle')

def _sam_to_bam(sam):
    '''
    From the location of a .sam file, convert it to a bam file and retun the location
    '''
    if sam[-4:] != '.sam':
        print('Sam file needs to end in .sam')
        sys.exit()

    bam = sam[:-4] + '.bam'
    print("Converting {0} to {1}".format(sam, bam))
    cmd = ['samtools', 'view','-S','-b', sam, '-o', bam]
    call(cmd)

    return bam

def _sort_index_bam(bam, rm_ori=True):
    '''
    From a .bam file, sort and index it. Remove original if rm_ori

    Return path of sorted and indexed bam
    '''
    if bam[-4:] != '.bam':
        print('Bam file needs to end in .sam')
        sys.exit()

    print("sorting {0}".format(bam))
    sorted_bam = bam[:-4] + '.sorted.bam'
    cmd = ['samtools', 'sort', '-o', sorted_bam, bam]
    call(cmd)

    print("Indexing {0}".format(sorted_bam))
    cmd = ['samtools', 'index', sorted_bam, sorted_bam + '.bai']
    call(cmd)

    if rm_ori:
        print("Deleting {0}".format(bam))
        os.remove(bam)

    return sorted_bam

def _parse_Sdb(sdb):
    '''
    Add some information to sdb
    '''
    if len(sdb) == 0:
        return sdb

    sdb['baseCoverage'] = [sum([a,c,t,g]) for a,c,t,g in zip(sdb['A'],sdb['C'],sdb['T'],sdb['G'])]
    sdb['varBase'] = [['A','C','T','G'][[a,c,t,g].index(sorted([a,c,t,g])[-1])]\
                      if ['A','C','T','G'][[a,c,t,g].index(sorted([a,c,t,g])[-1])] != r \
                      else ['A','C','T','G'][[a,c,t,g].index(sorted([a,c,t,g])[-2])] \
                      for a,c,t,g,r in zip(sdb['A'], sdb['C'], sdb['T'], sdb['G'], sdb['refBase'])]
    sdb['varFreq'] = [[a,c,t,g][['A','C','T','G'].index(v)]/s for a,c,t,g,v,s in zip(\
                        sdb['A'], sdb['C'], sdb['T'], sdb['G'], sdb['varBase'], sdb['baseCoverage'])]
    sdb['refFreq'] = [[a,c,t,g][['A','C','T','G'].index(v)]/s for a,c,t,g,v,s in zip(\
                        sdb['A'], sdb['C'], sdb['T'], sdb['G'], sdb['refBase'], sdb['baseCoverage'])]
    return sdb

def profile_bam(bam, **kwargs):
    '''
    Return a dataframe with the complete coverage exhaustion of each scaffold

    If profileANI, also find SNPs and report consensus ANI (slower)

    Return Bdb, Sdb (Sdb will be blank if not profileANI)
    '''
    # get arguments
    fasta = kwargs.get('fasta', None)
    profileANI = kwargs.get('profileANI', True)
    minP = kwargs.get('minP', .8)
    minC = kwargs.get('minC', 5)

    # set up coverage dataframe
    table = defaultdict(list)
    Stable = False

    # set up .sam file
    samfile = pysam.AlignmentFile(bam)

    # set up .fasta file
    if profileANI:
        scaff2sequence = SeqIO.to_dict(SeqIO.parse(fasta, "fasta"))
        Stable = defaultdict(list) # Set up SNP table

    # Get scaffold2length from .bam file
    s2l = {}
    for dic in samfile.header['SQ']:
        s2l[dic['SN']] = dic['LN']

    # Iterate scaffolds
    for scaff in tqdm(s2l, desc='Scaffolds processed'):
        covs = [] # List of coverage values for this scaffold
        if profileANI:
            scaffSeq = scaff2sequence[scaff] # Dictionary of scaffold to seqIO object
            basesCounted = 0 # Count of bases that got through to SNP calling
            snpsCounted = 0 # Count of SNPs

        for pileupcolumn in samfile.pileup(scaff):
            if profileANI:
                # Iterate reads at this position to figure out basecounts
                counts = _get_base_counts(pileupcolumn)
                baseCoverage = sum(counts) # Determine base coverage (this one does NOT include N or indels)
                covs.append(baseCoverage)

                refBase = scaffSeq[pileupcolumn.pos] # Determine reference base
                snp = _call_SNP(counts, refBase, minC, minP) # Call SNP
                if snp == 2: # means base was not counted
                    continue

                elif snp == 1: # means this is a SNP
                    Stable['scaffold'].append(scaff)
                    Stable['position'].append(pileupcolumn.pos)
                    Stable['refBase'].append(refBase)
                    for b, c in zip(['A', 'C', 'T', 'G'], counts):
                        Stable[b].append(c)
                    snpsCounted += 1
                basesCounted += 1

            else: # Just determine the coverage - THIS WILL COUNT READS WITH N or INDEL
                covs.append(pileupcolumn.n)

        # backfill coverage with 0s to account for those with no reads mapping
        zeros = len(covs)
        covs = covs + ([0] * (s2l[scaff] - len(covs)))
        assert len(covs) == s2l[scaff], [covs, s2l[scaff]]

        # fill in all coverage information
        table['scaffold'].append(scaff)
        table['length'].append(s2l[scaff])
        table['breadth'].append(zeros/s2l[scaff])
        table['coverage'].append(np.mean(covs))
        table['median_cov'].append(int(np.median(covs)))
        table['std_cov'].append(np.std(covs))
        table['bases_w_0_coverage'].append(zeros)
        table['max_cov'].append(max(covs))
        table['min_cov'].append(min(covs))

        # fill in all SNP information
        table['SNPs'].append(snpsCounted)
        table['unmaskedBreadth'].append(basesCounted / s2l[scaff])
        if basesCounted == 0:
            table['ANI'].append(0)
        else:
            table['ANI'].append((basesCounted - snpsCounted)/ basesCounted)

    Cdb = pd.DataFrame(table)
    if Stable != False:
        try:
            Sdb = pd.DataFrame(Stable)
            Sdb['scaffold'] = Sdb['scaffold'].astype('category')
            Sdb['refBase'] = Sdb['refBase'].astype('category')
        except TypeError:
            print("No SNPs detected!")
            Sdb = None
    else:
        Sdb = None

    return Cdb, Sdb

P2C = {'A':0, 'C':1, 'T':2, 'G':3}
def _get_base_counts(pileupcolumn):
    '''
    From a pileupcolumn object, return a list with the counts of [A, C, T, G]
    '''
    counts = [0,0,0,0]
    for pileupread in pileupcolumn.pileups:
        if not pileupread.is_del and not pileupread.is_refskip:
            try:
                counts[P2C[pileupread.alignment.query_sequence[pileupread.query_position]]] += 1
            except KeyError: # This would be like an N or something not A/C/T/G
                pass
    return counts

def _call_SNP(counts, base, minC=5, minP=0.8):
    '''
    From a list of counts [A,C,T,G] and the reference base, determine if SNP or not

    Args:
        minC = minimum coverage (if below this, return 2)
        minP = if reference base has less than this fraction, call SNP

    Return:
        0 if not SNP
        1 if SNP
        2 if unCounted
    '''
    baseCoverage = sum(counts)
    if baseCoverage < minC:
        return 2

    try:
        refcount = counts[P2C[base]]
    except KeyError:
        return 2

    if (refcount / baseCoverage) < minP:
        return 1
    return 0

def gen_genome_table(db, stb, min_c = 1):
    gdb = db.copy()
    gdb['bin'] = gdb['scaffold'].map(stb)
    gdb = gdb[gdb['bin'] == gdb['bin']]

    Table = defaultdict(list)
    for binn, d in gdb.groupby('bin'):
        length = d['length'].sum()

        Table['length'].append(length)
        Table['genome'].append(binn)
        Table['breadth'].append(sum(c * l for c,l in zip(d['breadth'], d['length'])) / d['length'].sum())
        Table['coverage'].append(sum(c * l for c,l in zip(d['coverage'], d['length'])) / d['length'].sum())

    Gdb = pd.DataFrame(Table)
    return(Gdb)

def parse_args(args):
    '''
    Parse command line arguemnts
    '''
    parser = argparse.ArgumentParser(description = \
            'Profile the strain in a given sample (v {0})'.format(__version__),\
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    InpArgs = parser.add_argument_group('INPUTS')
    InpArgs.add_argument(\
            '-b', action = 'store', \
            help = 'bam file')
    InpArgs.add_argument(\
            '-s', action = 'store', \
            help = 'sam file (will convert to bam before running)')
    InpArgs.add_argument(\
            '-g','--stb', action = 'store', \
            help = 'scaffold to bin file (to run on whole-genome level)')
    InpArgs.add_argument(\
            '-f', '--fasta', action = 'store', \
            help = 'fasta file (required for ANI profiling)')

    OptArgs = parser.add_argument_group('OPTIONS')
    OptArgs.add_argument(\
            '-p', '--minP', default=.8, type=float,\
            help = 'If reference base has less than this fraction of support, call SNP')
    OptArgs.add_argument(\
            '-c', '--minC', default=5, type=int,\
            help = 'Minumum coverage to call SNPs')

    OutArgs = parser.add_argument_group('OUTPUTS')
    OutArgs.add_argument(\
            '-o', action = 'store', default='strainProfile',\
            help = 'output base name')

    return vars(parser.parse_args(args))

if __name__ == '__main__':
    args = parse_args(sys.argv[1:])
    controller = Controller().parseArguments(args)
