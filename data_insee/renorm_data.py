import pandas as pd
import numpy as np
#from plotnine import *
import matplotlib.pyplot as plt

### Custom functions
def classing(row):
    if row['age'] in range(0,5):
        return 1
    elif row['age'] in range(5,10):
        return 2
    elif row['age'] in range(10,15):
        return 3
    elif row['age'] in range(15,20):
        return 4
    elif row['age'] in range(20,25):
        return 5
    elif row['age'] in range(25,30):
        return 6
    elif row['age'] in range(30,35):
        return 7
    elif row['age'] in range(35,40):
        return 8
    elif row['age'] in range(40,45):
        return 9
    elif row['age'] in range(45,50):
        return 10
    elif row['age'] in range(50,55):
        return 11
    elif row['age'] in range(55,60):
        return 12
    elif row['age'] in range(60,65):
        return 13
    elif row['age'] in range(65,70):
        return 14
    elif row['age'] in range(70,75):
        return 15
    elif row['age'] in range(75,80):
        return 16
    elif row['age'] in range(80,85):
        return 17
    elif row['age'] in range(85,90):
        return 18
    elif row['age'] in range(90,95):
        return 19
    elif row['age'] in range(95,100):
        return 20
    else:
        return 21

def classing_norm(row):
    if row['classes'] in range(0,5):
        return 'xpop'
    else:
        return row['classes']-4

def classing_labs(row):
    if row['age'] in range(0,5):
        return '0-4'
    elif row['age'] in range(5,10):
        return '5-9'
    elif row['age'] in range(10,15):
        return '10-14'
    elif row['age'] in range(15,20):
        return '15-19'
    elif row['age'] in range(20,25):
        return '20-24'
    elif row['age'] in range(25,30):
        return '25-29'
    elif row['age'] in range(30,35):
        return '30-34'
    elif row['age'] in range(35,40):
        return '35-39'
    elif row['age'] in range(40,45):
        return '40-44'
    elif row['age'] in range(45,50):
        return '45-49'
    elif row['age'] in range(50,55):
        return '50-54'
    elif row['age'] in range(55,60):
        return '55-59'
    elif row['age'] in range(60,65):
        return '60-64'
    elif row['age'] in range(65,70):
        return '65-69'
    elif row['age'] in range(70,75):
        return '70-74'
    elif row['age'] in range(75,80):
        return '75-79'
    elif row['age'] in range(80,85):
        return '80-84'
    elif row['age'] in range(85,90):
        return '85-89'
    elif row['age'] in range(90,95):
        return '90-94'
    elif row['age'] in range(95,100):
        return '95-99'
    else:
        return '100+'

def cohorter(row):
    time = 1900 + np.arange(start = 0, stop = 201, step = 5)
    for i in time:
        if row['year'] in range(i, i+5):
            return str.join('-', [str(i), str(i+5)])


#### impoprt raw data and minimal handle
proj = pd.read_excel('H:/Il mio Drive/datawork/prolonged_central_demographie_insee_nov21.xlsx', 
                        sheet_name='population', 
                        skiprows=1, 
                        nrows=107, 
                        usecols='A:FE',
                        header=0)
                        
proj.rename(columns={'Âge au 1er janvier':'age'}, inplace=True)

projlong = pd.melt(proj, 
                    id_vars='age', 
                    value_name='pop', 
                    var_name='year') \
                .reset_index(drop=True)


retrodata = pd.read_excel('H:/Il mio Drive/datawork/old_data.xlsx', 
                            sheet_name='midway',
                            usecols='A:CY',
                            header=0,
                            nrows=62) \
                                .rename(columns={'100 et +':100, 
                                            'Population totale': 'tot',
                                            'age':'year'})

retrolong = pd.melt(retrodata, id_vars=['year'], 
                    value_name='pop', 
                    var_name='age') \
                .reset_index(drop=True)

# start data wrangling

### pull out and separate totals
# projections data
projtot = projlong.loc[(projlong['age'] == 'Total') & (projlong['year'] <= 2100)]
projtot.loc[:,'age'] = 'tot'
projlong = projlong.loc[(projlong['age'] != 'Total') & (projlong['year'] <= 2100)]

# historical data
retrotot = retrolong[retrolong['age'] == 'tot']
retrolong = retrolong[retrolong['age'] != 'tot']

totpop = pd.concat([projtot, retrotot], ignore_index=True) \
    .drop(columns=['age']) \
        .sort_values(by=['year']) \
            .drop_duplicates()

# concat asap
popdata = pd.concat([projlong, retrolong], 
                    ignore_index=True,
                    axis=0) \
                        .reset_index(drop=True) \
                            .sort_values(by=['year', 'age']) \
                                .drop_duplicates()

# create classes
popdata['classes'] = popdata.apply(classing, axis=1)

# create normed classes
popdata['classes_norm'] = popdata.apply(classing_norm, axis=1)

# create labels
popdata['label'] = popdata.apply(classing_labs, axis=1)

# create time cohorts
popdata['cohort'] = popdata.apply(cohorter, axis=1)

# pull apart multi-dict
transler = popdata[['classes', 'classes_norm', 'label']] \
                .drop_duplicates() \
                    .reset_index(drop=True)

# sum up by normed classes
popdata_yearsum = popdata.drop(['classes', 'label', 'age'], axis=1) \
                            .groupby(['classes_norm', 'year', 'cohort']) \
                                .sum() \
                                    .reset_index(drop=False)

# avgs over the cohorts
popdata_cohortmean = popdata_yearsum.drop(['year'], axis=1) \
                                    .groupby(['classes_norm', 'cohort']) \
                                        .mean() \
                                            .reset_index(drop=False)

# pull xpop value in 2000 to normalise
xpop2000 = popdata_cohortmean.loc[(popdata_cohortmean['cohort'] == '2000-2005') & 
                                    (popdata_cohortmean['classes_norm'] == 'xpop')] \
                                        .loc[:, 'pop'].values[0]

# normalise data
popdata_norm = popdata_cohortmean
popdata_norm['pop_norm'] = popdata_norm['pop'] / xpop2000

popdata_norm_wide = popdata_norm\
    .drop(['pop'], axis=1) \
    .pivot(index='classes_norm', 
           columns='cohort', 
           values='pop_norm') \
    .reset_index()

popdata_norm_wide.to_excel('D:/emanu/OneDrive/Matlab/cdc_model/data_insee/popdata_norm_wide.xlsx', index = False)

# recompute totals
popdata_total = popdata_norm \
                    .groupby(['cohort']) \
                        .sum() \
                            .reset_index(drop=False) \
                                .rename(columns={'pop_norm': 'pop_norm_tot', 
                                                 'pop': 'pop_tot'})

popdata_total.to_excel('D:/emanu/OneDrive/Matlab/cdc_model/data_insee/popdata_total.xlsx', index=False)

#### test area for plots ######################################################
fig, ax = plt.subplots(figsize=(15,10))
for i in popdata_norm['classes_norm'].unique():
    ss = popdata_norm.loc[popdata_norm['classes_norm'] == i]
    plt.plot(ss['cohort'], ss['pop_norm'], label=i)

plt.legend()
plt.show()

import seaborn as sns
sns.lineplot(x='classes_norm', y='pop_norm', hue='cohort', data=popdata_norm.loc[popdata_norm['classes_norm'] != 'xpop'])
plt.show()

def plotlab(row):
    if row['classes_norm'] in range(1, 8):
        return 'A'
    elif row['classes_norm'] in range(8, 18):
        return 'I'
    else:
        return 'J'

popdata_norm['lab'] = popdata_norm.apply(plotlab, axis=1)

poplot = popdata_norm \
    .drop('pop_norm', axis = 1) \
        .groupby(['cohort', 'lab']) \
            .sum() \
                .reset_index(drop=False) \
                    .rename(columns={'pop': 'pop_lab'})

poplot2 = poplot.groupby(['cohort']) \
                    .sum() \
                        .reset_index(drop=False) \
                            .rename(columns={'pop_lab': 'pop'})

poplot3 = poplot2.merge(poplot, on='cohort')
poplot3['perc'] = poplot3['pop_lab'] / poplot3['pop']
poplot4 = poplot3.drop(['pop', 'pop_lab'], axis=1) \
                    .pivot(index='cohort',
                            columns='lab',
                            values='perc') \
                        .reset_index()





fig_new_pc = plt.figure(figsize=(16,9))
plt.stackplot(poplot4['cohort'], poplot4['J'], poplot4['A'], poplot4['I'], labels=['Pop. Jeune', 'Pop. Active', 'Pop. Inactive'])
fig_new_pc.autofmt_xdate()
plt.legend(loc='upper left')
plt.title('Population par groupe d\'âge en % - INSEE Nov. 2021')
#plt.show()

plt.savefig('D:/emanu/OneDrive/Matlab/cdc_model/plots/pop_pourc_2021.eps', 
            format='eps', 
            dpi=1000)

poplot5 = poplot3 \
    .drop(['perc', 'pop'], axis=1) \
        .pivot(index='cohort',
                columns='lab',
                values='pop_lab') \
                    .reset_index()

fig_new_lvl = plt.figure(figsize=(16,9))
plt.stackplot(poplot5['cohort'], poplot5['J'], poplot5['A'], poplot5['I'], labels=['Pop. Jeune', 'Pop. Active', 'Pop. Inactive'])
fig_new_lvl.autofmt_xdate()
plt.legend(loc='upper left')
plt.title('Population par groupe d\'âge en niveau - INSEE Nov. 2021')
#plt.show()

plt.savefig('D:/emanu/OneDrive/Matlab/cdc_model/plots/pop_niveau_2021.eps', 
            format='eps', 
            dpi=1000)

#### plots from old data ######################################################
oldpop = pd.read_excel('H:/Il mio Drive/datawork/Data_Population_old.xlsx',
                       sheet_name='ww',
                       skiprows=47,
                       usecols='A:W',
                       nrows=41) \
                           .rename(columns={'Unnamed: 0': 'year', 'Population totale': 'tot'}) \
                               .melt(id_vars=['year'],
                                     value_name='pop_norm',
                                     var_name='age')

backnorm = 3462550

oldpop['pop'] = oldpop['pop_norm']*backnorm
oldpoptot = oldpop.loc[oldpop['age'] == 'tot']
oldpop = oldpop.loc[oldpop['age'] != 'tot']

oldpop['cohort'] = oldpop.apply(cohorter, axis=1)
oldpop['classes'] = oldpop.apply(classing, axis=1)
oldpop['classes_norm'] = oldpop.apply(classing_norm, axis=1)
oldpop.drop('classes', axis=1, inplace=True)
oldpop['lab'] = oldpop.apply(plotlab, axis=1)

plo1 = oldpop \
    .drop(['pop_norm', 'age', 'year'], axis=1) \
        .groupby(['cohort', 'lab']) \
            .sum() \
                .reset_index(drop=False) \
                    .rename(columns={'pop': 'pop_lab'})

plo2 = plo1.groupby(['cohort']) \
                    .sum() \
                        .reset_index(drop=False) \
                            .rename(columns={'pop_lab': 'pop'})

plo3 = plo2.merge(plo1, on='cohort')
plo3['perc'] = plo3['pop_lab'] / plo3['pop']
plo4 = plo3.drop(['pop', 'pop_lab'], axis=1) \
                    .pivot(index='cohort',  
                            columns='lab',  
                            values='perc') \
                        .reset_index()

fig_old_pc = plt.figure(figsize=(16,9))
plt.stackplot(plo4['cohort'], plo4['J'], plo4['A'], plo4['I'], labels=['Pop. Jeune', 'Pop. Active', 'Pop. Inactive'])
fig_old_pc.autofmt_xdate()
plt.legend(loc='upper left')
plt.title('Population par groupe d\'âge en % - INSEE Nov. 2020')
#plt.show()

plt.savefig('D:/emanu/OneDrive/Matlab/cdc_model/plots/pop_pourc_2020.eps',
            format='eps',
            dpi=1000)

plo5 = plo3 \
    .drop(['perc', 'pop'], axis=1) \
        .pivot(index='cohort',
                columns='lab',
                values='pop_lab') \
                    .reset_index()

fig_old_lvl = plt.figure(figsize=(16,9))
plt.stackplot(plo5['cohort'], plo5['J'], plo5['A'], plo5['I'], labels=['Pop. Jeune', 'Pop. Active', 'Pop. Inactive'])
fig_old_lvl.autofmt_xdate()
plt.legend(loc='upper left')
plt.title('Population par groupe d\'âge en niveau - INSEE Nov. 2020')
#plt.show()

plt.savefig('D:/emanu/OneDrive/Matlab/cdc_model/plots/pop_niveau_2020.eps',
            format='eps',
            dpi=1000)