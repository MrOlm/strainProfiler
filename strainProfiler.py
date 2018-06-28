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
__version__ = "0.2.0"
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
            if (os.path.exists(bam + '.bai')) | ((os.path.exists(bam[:-4] + '.bai'))):
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

def _update_cov_table(table, covs, lengt, scaff):
    '''
    Add information to the table

    Args:
        table: table to add to
        covs: list of coverage values
        lengt: length of scaffold
        scaff: name of scaffold
    '''
    nonzeros = len(covs)
    zeros = (lengt - len(covs))
    covs = covs + ([0] * (lengt - len(covs)))

    assert len(covs) == lengt, [covs, lengt]

    # fill in all coverage information
    table['scaffold'].append(scaff)
    table['length'].append(lengt)
    table['breadth'].append(nonzeros/lengt)
    table['coverage'].append(np.mean(covs))
    table['median_cov'].append(int(np.median(covs)))
    table['std_cov'].append(np.std(covs))
    table['bases_w_0_coverage'].append(zeros)
    table['max_cov'].append(max(covs))
    table['min_cov'].append(min(covs))

def _update_covT_table(table, covT, lengt, scaff, debug=False):
    '''
    Add information to the table

    Args:
        table: table to add to
        covs: list of coverage values
        lengt: length of scaffold
        scaff: name of scaffold
    '''
    for mm in sorted(list(covT.keys())):
        covs = _mm_counts_to_counts(covT, mm)
        if covs == [0,0,0,0]:
            covs = [0]*lengt
        nonzeros = np.count_nonzero(covs)
        zeros = lengt - nonzeros

        assert len(covs) == lengt, [covs, lengt, mm]

        # fill in all coverage information
        table['scaffold'].append(scaff)
        table['length'].append(lengt)
        table['breadth'].append(nonzeros/lengt)
        table['coverage'].append(np.mean(covs))
        table['median_cov'].append(int(np.median(covs)))
        table['std_cov'].append(np.std(covs))
        table['bases_w_0_coverage'].append(zeros)
        table['max_cov'].append(max(covs))
        table['min_cov'].append(min(covs))
        table['mm'].append(mm)

    if debug == True:
        covs = _mm_counts_to_counts(covT, max(list(covT.keys())))
        zero_pos = [i+1 for i,x in enumerate(covs) if x==0]
        print(zero_pos)
        print(len(zero_pos))

def _update_snp_table_T(Stable, basesCounted, snpsCounted, refBase, MMcounts,\
        pos, scaff, minC=5, minP=.8):
    '''
    Add information to SNP table
    '''
    x = _mm_counts_to_counts(MMcounts)
    #print("inner type: {0}".format(type(x)))
    #print("mms to check: " + str(sorted(list(MMcounts.keys()))))
    for mm in list(MMcounts.keys()):
        #print('checking {0}'.format(mm))
        counts = _mm_counts_to_counts(MMcounts, mm)
        #print("iii type: {0}".format(type(counts)))
        snp = _call_SNP(counts, refBase, minC, minP) # Call SNP

        #print(type(counts))
        #print("pos {0} has {1} cov at {2}mm and is SNP-type {3}".format(pos, counts.sum(), mm, snp))
        #print(counts)

        if snp == 2: # means base was not counted
           continue

        elif snp == 1: # means this is a SNP
           Stable['scaffold'].append(scaff)
           Stable['position'].append(pos)
           Stable['refBase'].append(refBase)
           for b, c in zip(['A', 'C', 'T', 'G'], counts):
               Stable[b].append(c)
           Stable['mm'].append(mm)

           snpsCounted[mm] += 1

        #print("{0} is True".format(pos))
        basesCounted[mm][pos] = True # count everything that's not skipped

def _calc_counted_bases(basesCounted, maxMM):
    counts = None
    for mm, count in [(mm, count) for mm, count in basesCounted.items() if mm <= maxMM]:
        if counts is None:
            counts = count
        else:
            counts = np.add(counts, count)

    if counts is None:
        return 0

    else:
        return counts.sum()

def _update_snp_covT_table(table, snpsCounted, basesCounted, lengt, scaff, covT,
            minCov):
    # fill in all SNP information
    for mm in sorted(list(covT.keys())):
        covs = _mm_counts_to_counts(covT, mm)
        if covs == [0,0,0,0]:
            counted_basesO = 0
        else:
            zeros = (covs < minCov).sum()
            counted_basesO = lengt - zeros

        counted_snps = sum([snpsCounted[m] for m in list(snpsCounted.keys())\
            if m <= mm])
        counted_bases = _calc_counted_bases(basesCounted, mm)

        # print(basesCounted[mm])
        # print(len(basesCounted[mm]))
        # print(basesCounted[mm].sum())
        # print(counted_bases, counted_basesO)

        table['scaffold'].append(scaff)
        table['mm'].append(mm)
        table['SNPs'].append(counted_snps)
        table['unmaskedBreadth'].append(counted_bases / lengt)
        if counted_bases == 0:
            table['ANI'].append(0)
        else:
            table['ANI'].append((counted_bases - counted_snps)/ counted_bases)

def _make_snp_table(Stable):
    if Stable != False:
        try:
            Sdb = pd.DataFrame(Stable)
            Sdb['scaffold'] = Sdb['scaffold'].astype('category')
            Sdb['refBase'] = Sdb['refBase'].astype('category')
        except KeyError:
            print("No SNPs detected!")
            Sdb = pd.DataFrame()
    else:
        Sdb = pd.DataFrame()

    return Sdb

def _mm_counts_to_counts(MMcounts, maxMM=100):
    '''
    Take mm counts and return just counts
    '''
    counts = None
    for mm, count in [(mm, count) for mm, count in MMcounts.items() if mm <= maxMM]:
        if counts is None:
            counts = count
        else:
            counts = np.add(counts, count)

    if counts is None:
        return np.zeros(4, dtype=int)

    else:
        return counts

def _update_covT(covT, MMcounts, position):
    '''
    Update covT at this position
    '''
    for mm, count in MMcounts.items():
        covT[mm][position] = sum(count)

def run_up_NaN(odb, cols, on='scaffold'):
    '''
    Take NaN values and fill them with the column above.

    For example, if you have mm of [0,1,2,3], and breadth of [.9, .92, NaN, .94],
    change the breadth to [.9, .92, .92, .94]

    If you have [0,1,2,3] and [NaN, NaN, .92, .94], it'll change it to:
    [0, 0, .92, .94]

    NOTE: Must be sorted / indexed in the order that you want run up

    Args:
        odb: original dataframe
        cols: columns to do this filling on
        on: the "key" of the dataframe. If this is all NaN, fill with 0s

    Returns:
        DataFrame: new dataframe
    '''
    Fdb = odb.copy()
    for scaff, db in odb.groupby(on):
        # handle edge case where all are NaN
        if len(db[cols[0]].dropna()) == 0:
            for i, row in db.iterrows():
                for col in cols:
                    Fdb.at[i, col] = 0
            continue

        # do the actual run-up
        top = True
        for i, row in db.iterrows():
            # hangle edge case where top values are NaN
            if top & np.isnan(row[cols[0]]):
                for col in cols:
                    Fdb.at[i, col] = 0
                continue
            else:
                top = False

            # The normal run-up case
            if np.isnan(row['ANI']):
                for col in cols:
                    Fdb.at[i, col] = Fdb.at[i-1, col]

    return Fdb


def _merge_tables_special(Cdb, Adb):
    FIX_COLS = ['ANI', 'SNPs', 'unmaskedBreadth']

    # Handle edge-case where there are no SNPs
    if len(Adb) == 0:
        for c in FIX_COLS:
            Cdb[c] = 0
        return Cdb

    # Make sure anything with a coverage has a SNP
    assert len(set(Adb['mm']) - set(Cdb['mm'])) == 0

    # Do initial merge
    Cdb = pd.merge(Cdb, Adb, on=['scaffold', 'mm'], how='outer')
    Cdb = Cdb.sort_values(['scaffold', 'mm'])
    Cdb = Cdb.reset_index(drop=True)

    # Run up NaN values. For example, in the cases
    Fdb = run_up_NaN(Cdb, FIX_COLS, on='scaffold')

    # Might as well adjust some datatypes
    pass

    return Fdb

def profile_bam(bam, fasta, **kwargs):
    '''
    Return a dataframe with the complete coverage exhaustion of each scaffold

    Bdb = coverage information on all scaffolds
    Sdb = SNP information

    Return Bdb, Sdb
    '''
    # get arguments
    minP = kwargs.get('minP', .8)
    minC = kwargs.get('minC', 5)

    # initialize
    table = defaultdict(list) # set up coverage dataframe
    Atable = defaultdict(list) # set up ANI dataframe
    Stable = defaultdict(list) # Set up SNP table

    samfile = pysam.AlignmentFile(bam) # set up .sam file

    scaff2sequence = SeqIO.to_dict(SeqIO.parse(fasta, "fasta")) # set up .fasta file
    s2l = {s:len(scaff2sequence[s]) for s in list(scaff2sequence.keys())} # Get scaffold2length

    # Iterate scaffolds
    for scaff in tqdm(s2l, desc='Scaffolds processed'):
        covT = defaultdict(lambda:np.zeros(s2l[scaff], dtype=int)) # Dictionary of mm -> positional coverage
        basesCounted = defaultdict(lambda:np.zeros(s2l[scaff], dtype=bool)) # Count of bases that got through to SNP calling
        snpsCounted = defaultdict(int) # Count of SNPs

        for pileupcolumn in samfile.pileup(scaff):
            # Iterate reads at this position to figure out basecounts
            # note: pileupcolumn.pos is 0-based
            MMcounts = _get_base_counts_mm(pileupcolumn)
            _update_covT(covT, MMcounts, pileupcolumn.pos)

            # Call SNPs
            _update_snp_table_T(Stable, basesCounted,\
                    snpsCounted, scaff2sequence[scaff][pileupcolumn.pos], MMcounts,\
                    pileupcolumn.pos, scaff, minC=minC, minP=minP)

        # Update coverage table
        _update_covT_table(table, covT, s2l[scaff], scaff)

        # Update ANI table
        _update_snp_covT_table(Atable, snpsCounted, basesCounted, s2l[scaff], \
                scaff, covT, minC)

    # Make coverage table
    Cdb = _merge_tables_special(pd.DataFrame(table), pd.DataFrame(Atable))

    # Maks SNP table
    Sdb = _make_snp_table(Stable)

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

P2C = {'A':0, 'C':1, 'T':2, 'G':3}
def _get_base_counts_mm(pileupcolumn):
    '''
    From a pileupcolumn object, return a dictionary of readMismatches ->
        list with the counts of [A, C, T, G]
    '''
    table = defaultdict(lambda:np.zeros(4, int))
    for pileupread in pileupcolumn.pileups:
        if not pileupread.is_del and not pileupread.is_refskip:
            try:
                table[pileupread.alignment.get_tag('NM')]\
                [P2C[pileupread.alignment.query_sequence[pileupread.query_position]]] += 1
            except KeyError: # This would be like an N or something not A/C/T/G
                pass
    return table

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
