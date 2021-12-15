import pandas as pd
import numpy as np
from plotnine import *
import matplotlib.pyplot as plt

def classing(row):
    if row['age'] in range(0,5):
        return 'pop1'
    elif row['age'] in range(5,10):
        return 'pop2'
    elif row['age'] in range(10,15):
        return 'pop3'
    elif row['age'] in range(15,20):
        return 'pop4'
    elif row['age'] in range(20,25):
        return 'pop5'
    elif row['age'] in range(25,30):
        return 'pop6'
    elif row['age'] in range(30,35):
        return 'pop7'
    elif row['age'] in range(35,40):
        return 'pop8'
    elif row['age'] in range(40,45):
        return 'pop9'
    elif row['age'] in range(45,50):
        return 'pop10'
    elif row['age'] in range(50,55):
        return 'pop11'
    elif row['age'] in range(55,60):
        return 'pop12'
    elif row['age'] in range(60,65):
        return 'pop13'
    elif row['age'] in range(65,70):
        return 'pop14'
    elif row['age'] in range(70,75):
        return 'pop15'
    elif row['age'] in range(75,80):
        return 'pop16'
    elif row['age'] in range(80,85):
        return 'pop17'
    elif row['age'] in range(85,90):
        return 'pop18'
    elif row['age'] in range(90,95):
        return 'pop19'
    elif row['age'] in range(95,100):
        return 'pop20'
    else:
        return 'pop21'

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

proj = pd.read_excel('H:/Il mio Drive/datawork/prolonged_central_demographie_insee_nov21.xlsx', 
                        sheet_name='population', 
                        skiprows=1, 
                        nrows=107, 
                        usecols='A:FE',
                        header=0)

proj.rename(columns={'Âge au 1er janvier':'age'}, inplace=True)

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

projlong = pd.melt(proj, id_vars='age', value_name='pop', var_name='year').reset_index(drop=True)

# pull out totals
# projections data
projtot = projlong.loc[(projlong['age'] == 'Total') & (projlong['year'] <= 2100)]
projtot.loc[:, 'age'] = 'tot'
projlong = projlong.loc[(projlong['age'] != 'Total') & (projlong['year'] <= 2100)]
# historical data
retrotot = retrolong[retrolong['age'] == 'tot']
retrolong = retrolong[retrolong['age'] != 'tot']

# create classes
projlong['classes'] = projlong.apply(classing, axis=1)
retrolong['classes'] = retrolong.apply(classing, axis=1)

# create labels
projlong['labels'] = projlong.apply(classing_labs, axis=1)
retrolong['labels'] = retrolong.apply(classing_labs, axis=1)

# sum up by class and year
projlong = projlong.groupby(['year', 'classes']).sum().reset_index()
retrolong = retrolong.groupby(['year', 'classes']).sum().reset_index()

# reclassify years
projlong['cohort'] = projlong.apply(cohorter, axis=1)
retrolong['cohort'] = retrolong.apply(cohorter, axis=1)

projtot['cohort'] = projtot.apply(cohorter, axis=1)
retrotot['cohort'] = retrotot.apply(cohorter, axis=1)

total_pop = pd.concat([retrotot, projtot], axis=0).sort_values(by=['cohort', 'year']).drop(columns = 'age', axis = 1)


# sort data
projlong.sort_values(by=['year', 'classes'], inplace=True)
retrolong.sort_values(by=['year', 'classes'], inplace=True)

# stack them up
popdata = pd.concat([projlong.drop(columns = 'year', axis = 1),
                     retrolong.drop(columns = 'year', axis = 1)], 
                     ignore_index=True) \
            .groupby(['cohort', 'classes']) \
            .mean('pop') \
            .reset_index() \
            .pivot(index='classes', columns='cohort', values='pop') \
            .reset_index() 

popdata.to_csv('H:/Il mio Drive/datawork/popdata.csv', index=False)

total_pop = total_pop \
    .groupby(['cohort']) \
        .mean() \
            .reset_index()

total_pop.to_csv('H:/Il mio Drive/datawork/total_pop.csv', index=False)

# popdata['class_label'] = popdata.apply(classing_labs, axis=1)


(
    ggplot(popdata) 
    + aes(x='cohort', y='pop', color='classes', group='classes') 
    + geom_line()
)

(
    ggplot(total_pop)
    + aes(x = 'year', y = 'pop', color = 'cohort')
    + geom_line()
)

total_pop.groupby(['cohort']).mean(['pop']).reset_index().plot(x='cohort', y='pop', kind='line')
plt.show()

popdata[popdata["classes"] == 'pop1'].melt(id_vars = 'classes', value_name='ii', var_name='rr').plot(x = 'rr', y = 'ii', kind = 'line')
plt.show()