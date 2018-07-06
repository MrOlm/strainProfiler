# strainProfiler
Program to analyze strain-level diversity within a population

# Functionality of version 0.3

strainProfiler is able to take a .bam file and determine the breadth, coverage, ANI, ect. at all levels of read-mismatch.
So for 0, 1, 2, 3, ect. mismatches per read, you can figure out these things.

### Output of version 0.3

It produces `_scaffoldTable`, `_snpLocations.pickle`, and `.pickle`. 

`_scaffoldTable` contains the **cumulative** information for each scaffold at each mismatch level.

`_snpLocations.pickle` contains the **cumulative** SNP information for each scaffold.

`.pickle` contains the base information that was used to make those tables. It's size is dictated by the size of the .fasta
file being mapped to, and it can get pretty big. The `--lightRAM` option should make it a good bit smaller


### Other useful python methods

Loading the pickle file
```python
import strainProfiler
Sprofile = strainProfiler.SNVprofile().load(pick)
```

Transforming from scaffold-level to genome-level
```python
def makeGenomeWide_v2(sdb, stb, s2l=None):
    '''
    Make the scaffold table genome-wide from strainProfiler v0.2
    
    Args:
        sdb - the file ending with _scaffoldTable.csv
        stb - dicionary of scaffold -> bin
    '''
    gdb = sdb.copy()
    gdb['genome'] = gdb['scaffold'].map(stb)
    gdb['considered_length'] = [x*y for x,y in zip(gdb['unmaskedBreadth'], gdb['length'])]
    
    # Add blanks for scaffolds that aren't there
    #bdb = pd.DataFrame({s:})
    
    table = defaultdict(list)
    for genome, mm, db in interate_sdb_mm(gdb, s2l=s2l, stb=stb):
        table['genome'].append(genome)
        table['mm'].append(mm)

        # Sum columns
        for col in ['bases_w_0_coverage', 'length']:
            table[col].append(db[col].sum())
            
        # Max columns
        for col in ['SNPs']:
            table[col].append(db[col].max())

        # Weighted average columns
        for col in ['breadth', 'coverage', 'std_cov', 'unmaskedBreadth']:
            table[col].append(sum(x * y for x, y in zip(db[col], db['length'])) / sum(db['length']))

        # Special weighted average
        db['considered_length'] = [x*y for x,y in zip(db['unmaskedBreadth'], db['length'])]
        considered_leng = db['considered_length'].sum()
        
        if considered_leng != 0:
            table['ANI'].append((considered_leng - db['SNPs'].sum()) / considered_leng)
        else:
            table['ANI'].append(0)

        # Special
        table['max_cov'].append(db['max_cov'].max())
        table['min_cov'].append(db['min_cov'].min())
    
    return pd.DataFrame(table)

def interate_sdb_mm(sdb, on='genome', s2l=None, stb=None):
    '''python
    For the dataframe, iterate through each mm and a dataframe with ALL scaffolds at that mm level
    (including blanks)
    '''
    for g, db in sdb.groupby(on):
        if s2l == None:
            gs2l = db.drop_duplicates(subset=['scaffold', 'length'])\
                    .set_index('scaffold')['length'].to_dict()
        else:
            gs2l = {s:s2l[s] for s in [x for x, b in stb.items() if b == g]}
        mms = sorted(db['mm'].unique())
        for mm in mms:
            # get all the ones that you can
            dd = db[db['mm'] <= mm].sort_values('mm').drop_duplicates(subset='scaffold', keep='last')
            #print("mm={0}; len={1}; tl={2}".format(mm, len(dd), len(s2l.keys())))
            
            # backfill with blanks
            dd = _backfill_blanks(dd, gs2l)
            #print("mm={0}; len={1}; tl={2}".format(mm, len(dd), len(s2l.keys())))
            
            yield g, mm, dd
            
def _backfill_blanks(db, s2l):
    scaffs = list(set(s2l.keys()) - set(db['scaffold'].unique()))
    bdb = pd.DataFrame({'scaffold':scaffs, 'length':[s2l[s] for s in scaffs]})
    
    # make some adjustments
    bdb['bases_w_0_coverage'] = bdb['length']
    
    # append
    db = db.append(bdb)
    
    # fill 0
    return db.fillna(0)
```


Plotting the breadth / coverage over all mm levels:
```python
def mm_plot(db, left_val='breadth', right_val='coverage', title='',\
           maxmm=15):
    db = db.sort_values('mm')
    sns.set_style('white')

    # breadth
    fig, ax1 = plt.subplots()
    ax1.plot(db['mm'], db[left_val], ls='-', color='blue')
    if left_val == 'breadth':
        ax1.plot(db['mm'], estimate_breadth(db['coverage']), ls='--', color='lightblue')
    ax1.set_ylabel(left_val, color='blue')
    ax1.set_xlabel('read mismatches')
    ax1.set_ylim(0,1)

    # coverage
    ax2 = ax1.twinx()
    ax2.plot(db['mm'], db[right_val], ls='-', color='red')
    ax2.set_ylabel(right_val, color='red')
    ax2.set_ylim(0,)

    # asthetics
    plt.xlim(0, maxmm)

    plt.title(title)
    
def estimate_breadth(coverage):
    '''
    Estimate breadth based on coverage
    
    Based on the function breadth = -1.000 * e^(0.883 * coverage) + 1.000
    '''
    import numpy as np
    return (-1) * np.exp(-1 * ((0.883) * coverage)) + 1
```

Subsetting to only the highest level of mismatches
```python
tdb = Tdb.sort_values('mm').drop_duplicates(subset=['genome', 'sample'], keep='last')
```

Notebooks where I do this stuff:

*https://biotite.berkeley.edu/j/user/mattolm/notebooks/Infant_Eukaryotes/PublicationQuality/ControlEuks_4_removeContaminants_3.ipynb
*https://biotite.berkeley.edu/j/user/mattolm/notebooks/Infant_Eukaryotes/PublicationQuality/ControlEuks_2_Malassazia_5_strainProfilerv3_circos.ipynb
