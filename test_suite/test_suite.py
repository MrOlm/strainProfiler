#!/usr/bin/env python
'''
Run tests of strainProfiler
'''

import os
import glob
import shutil

import warnings
warnings.filterwarnings("ignore")

import pandas as pd
from subprocess import call

def load_data_loc():
    return os.path.join(str(os.getcwd()), \
        'test_files/')

def load_random_test_dir():
    loc = os.path.join(str(os.getcwd()), \
        'test_backend/testdir/')
    return loc

def get_script_loc(script):
    if script == 'strainProfiler.py':
        return os.path.join(str(os.getcwd()), \
            '../strainProfiler.py')

class test_strainProfiler_pile():
    def setUp(self):
        self.script = get_script_loc('strainProfiler.py')
        self.bam = load_data_loc() + \
            'N5_271_010G1_scaffold_min1000.fa-vs-N5_271_010G1.sorted.bam'
        self.fasta = load_data_loc() + \
            'N5_271_010G1_scaffold_min1000.fa'
        self.cc_solution = load_data_loc() + \
            'N5_271_010G1_scaffold_min1000.fa-vs-N5_271_010G1.bam.CB'
        self.pp_solution = load_data_loc() + \
            'N5_271_010G1_scaffold_min1000.fa-vs-N5_271_010G1.scaffoldPP_totals.csv'
        self.test_dir = load_random_test_dir()

        if os.path.isdir(self.test_dir):
            shutil.rmtree(self.test_dir)
        os.mkdir(self.test_dir)

    def tearDown(self):
        if os.path.isdir(self.test_dir):
            shutil.rmtree(self.test_dir)

    def run(self):
        self.setUp()
        self.test1()
        self.tearDown()

        self.setUp()
        self.test2()
        self.tearDown()

        self.setUp()
        self.test3()
        self.tearDown()

    def test1(self):
        '''
        Basic comprehensive test of values
        '''
        # Run program
        out_base = os.path.join(self.test_dir, 'test')
        cmd = [self.script, '-b', self.bam, '-o', out_base, '-f', self.fasta]
        call(cmd)

        # Load output
        Odb = pd.read_csv(out_base + '_scaffoldTable.csv')
        Sdb = pd.read_pickle(out_base + '_snpLocations.pickle')
        assert os.path.isfile(out_base + '_log')

        # Ensure internal consistancy
        assert Odb['ANI'].max() <= 1
        for scaff, db in Sdb.groupby('scaffold'):
            snps = Odb['SNPs'][Odb['scaffold'] == scaff].tolist()[0]
            assert snps == len(db), [snps, len(db)]

        # Compare to calculate_coverage
        Cdb = pd.read_csv(self.cc_solution)
        s2c = Cdb.set_index('scaffold')['coverage'].to_dict()
        s2b = Cdb.set_index('scaffold')['breadth'].to_dict()
        for scaff, db in Odb.groupby('scaffold'):
            assert (db['coverage'].tolist()[0] - s2c[scaff]) < .1, [db['coverage'].tolist()[0], s2c[scaff]]
            assert db['breadth'].tolist()[0] == s2b[scaff]

        # Compare to pileupProfile
        Pdb = pd.read_csv(self.pp_solution)
        Pdb['scaffold'] = [x[1:-3] for x in Pdb['genome']]
        Pdb['consensus_ANI'] = [0 if x != x else x for x in Pdb['consensus_ANI']]
        s2a = Pdb.set_index('scaffold')['consensus_ANI'].to_dict()
        for scaff, db in Odb.groupby('scaffold'):
            if scaff not in s2a:
                continue
            assert (db['ANI'].tolist()[0] - s2a[scaff]) <= 0, \
                [scaff, db['ANI'].tolist()[0], s2a[scaff]]

    def test2(self):
        '''
        Make sure command line arguments are followed
        '''
        # Run program with defaults
        out_base = os.path.join(self.test_dir, 'testD')
        cmd = [self.script, '-b', self.bam, '-o', out_base, '-f', self.fasta]
        call(cmd)

        # Run with adjusted parameters
        out_base = os.path.join(self.test_dir, 'testA')
        cmd = [self.script, '-b', self.bam, '-o', out_base, '-f', self.fasta,\
            '-p', '0.2', '-c', '10']
        call(cmd)

        # Load
        SD = pd.read_pickle(out_base[:-1] + 'D_snpLocations.pickle')
        SA = pd.read_pickle(out_base[:-1] + 'A_snpLocations.pickle')

        # Confirm
        assert abs(SD['refFreq'].max() - 0.8) < .01, SD['refFreq'].max()
        assert abs(SA['refFreq'].max() - 0.2) < .01, SA['refFreq'].max()

        assert SD['baseCoverage'].min() == 5, SD['baseCoverage'].min()
        assert abs(SA['baseCoverage'].min() - 10) <= 1, SA['baseCoverage'].min()

class test_strainProfiler_breadth():
    def setUp(self):
        self.script = get_script_loc('strainProfiler.py')
        self.bam = load_data_loc() + \
            'N5_225.delta.fasta-vs-N5_225_006G1.sub_sorted.bam'
        self.uns_bam = load_data_loc() + \
            'SP_CRL_000G1_concoct_7.fasta-vs-SP_CRL_022G2.bam'
        self.sam = load_data_loc() + \
            'Ig8144_scaffold_2_CLOSED_Huge_Phage.test.shrunk.sam'
        self.cov_sol = load_data_loc() + \
            'N5_225.delta.fasta-vs-N5_225_006G1.cov'
        self.bre_sol = load_data_loc() + \
            'N5_225.delta.fasta-vs-N5_225_006G1.cb.csv'
        self.test_dir = load_random_test_dir()

        if os.path.isdir(self.test_dir):
            shutil.rmtree(self.test_dir)
        os.mkdir(self.test_dir)

    def tearDown(self):
        if os.path.isdir(self.test_dir):
            shutil.rmtree(self.test_dir)

    def run(self):
        self.setUp()
        self.test1()
        self.test2()
        self.test3()
        self.tearDown()

    def test1(self):
        '''
        Basic comprehensive test of breadth calculation
        '''
        # Run program
        out_file = os.path.join(self.test_dir, 'test')
        cmd = [self.script, '-b', self.bam, '-o', out_file]
        print()
        print(' '.join(cmd))
        call(cmd)

        # Get output
        Odb = pd.read_csv(out_file)

        # Make sure of one that I confirmed in geneious
        b = Odb['breadth'][Odb['scaffold'] == "N5_225_000G1_scaffold_10275"].tolist()[0]
        assert b == (33/1002), b

        # Load solutions file with coverage of all scaffolds
        Cdb = pd.read_table(self.cov_sol, header=1, names=['scaffold', 'cov'])
        Cdb['len'] = [int(float(x.split(':')[1])) for x in Cdb['scaffold']]
        Cdb['scaffold'] = [x.split(':')[0] for x in Cdb['scaffold']]

        # Load solutions file with breadth of a subset of scaffolds
        Bdb = pd.read_csv(self.bre_sol)
        Bdb['scaffold'] = [x.split('.')[0] for x in Bdb['genome']]
        del Bdb['genome']
        del Bdb['Unnamed: 0']

        # Make sure coverage matches
        Odb = pd.merge(Cdb, Odb, on='scaffold')
        assert len(Odb[Odb['length'] != Odb['len']]) == 0
        Odb['c_diff'] = [abs(c-b) for c,b in zip(Odb['cov'], Odb['coverage'])]
        assert len(Odb[Odb['c_diff'] > 10]) == 2

        # Make sure breadth matches in all cases where you have it
        g2b = Bdb.set_index('scaffold')['breadth'].to_dict()
        Odb['correct_b'] = Odb['scaffold'].map(g2b)
        bdb = Odb[Odb['correct_b'] == Odb['correct_b']]
        bdb['diff'] = [abs(c-b) for c,b in zip(bdb['breadth'], bdb['correct_b'])]
        assert len(bdb[bdb['diff'] >= .1]) == 1

    def test2(self):
        '''
        Test the ability to make and sort .bam files from .sam
        '''
        # Copy sam to test dir
        new_sam = os.path.join(self.test_dir, os.path.basename(self.sam))
        shutil.copyfile(self.sam, new_sam)

        # Run program
        out_file = os.path.join(self.test_dir, 'test2')
        cmd = [self.script, '-s', new_sam, '-o', out_file]
        call(cmd)

        # Load output
        db = pd.read_csv(out_file)
        assert len(db) == 1

    def test3(self):
        '''
        Test the ability to make and sort .bam files from unsorted .bam
        '''
        # Copy bam to test dir
        uns_bam = os.path.join(self.test_dir, os.path.basename(self.uns_bam))
        shutil.copyfile(self.uns_bam, uns_bam)

        # Run program
        out_file = os.path.join(self.test_dir, 'test3')
        cmd = [self.script, '-b', uns_bam, '-o', out_file]
        call(cmd)

        # Load output
        db = pd.read_csv(out_file)
        assert len(db) == 897

if __name__ == '__main__':
    #test_strainProfiler_breadth().run()
    test_strainProfiler_pile().run()
    print('everything is working swimmingly!')
