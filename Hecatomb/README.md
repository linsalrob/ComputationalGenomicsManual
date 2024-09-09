# Hecatomb output analysis

# We have provided the data for Hecatomb that we are going to use:

- [bigtable](https://github.com/linsalrob/CF_Data_Analysis/raw/main/hecatomb/bigtable.tsv.gz?download=)
- [VMR table of viral families and hosts](https://github.com/linsalrob/CF_Data_Analysis/raw/main/hecatomb/VMR_MSL39_v1.ascii.tsv.gz?download=)
- [CF Metadata table](https://github.com/linsalrob/CF_Data_Analysis/raw/main/hecatomb/CF_Metadata_Table-2023-03-23.tsv.gz?download=)


# Jupyter analyais

## Load the data into pandas dataframes:

```
import os
import sys
import matplotlib.pyplot as plt
import scipy.stats as stats
import statsmodels.stats.api as sms
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
```

## Copy the files into your Jupyter notebooks and then find the files

```
os.chdir('Data')
os.listdir()
```

## Read the data into pandas data frames

```
data = pd.read_csv('bigtable.tsv.gz',compression='gzip',header=0,sep='\t')
metadata = pd.read_csv('CF_Metadata_Table-2023-03-23.tsv.gz',compression='gzip',header=0,sep='\t')
vmr = pd.read_csv('VMR_MSL39_v1.ascii.tsv.gz', compression='gzip',header=0,sep='\t')
```

Take some time to look at the tables. 

How many columns are there in each?
How many rows?
What are the data types in the rows and columns?

## First, filter the bigtable for _just_ the Viruses

```
viruses = data[(data.kingdom == "Viruses")]
```

## Group them by alignment type, length, and family

```
virusesGroup = viruses.groupby(by=['family','alnType','alnlen','pident'], as_index=False).count()
```

## Now use seaborn to make some pretty plots!

First, we start by creating a style that we like.

```
#styling
sizeScatter = 10 * virusesGroup['count']
sns.set_style("darkgrid")
sns.set_palette("colorblind")
sns.set(rc={'figure.figsize':(12,8)})

# and now plot the data
g = sns.FacetGrid(virusesGroup, col="family", col_wrap=10)
g.map_dataframe(sns.scatterplot, "alnlen", "pident", alpha=.1, hue="alnType", sizes=(100,500), size=sizeScatter)
plt.legend(bbox_to_anchor=(5.0,1), loc=0, borderaxespad=2,ncol=6, shadow=True, labelspacing=1.5, borderpad=1.5)
plt.show()
```


That's a lot of data!

### Restrict this plot to just E value < 1 x 10<sup>-20</sup>

```
virusesFiltered = viruses[viruses.evalue<1e-20]
virusesGroup = virusesFiltered.groupby(by=['family','alnType','alnlen','pident'], as_index=False).count()
g = sns.FacetGrid(virusesGroup, col="family", col_wrap=10)
g.map_dataframe(sns.scatterplot, "alnlen", "pident", alpha=.1, hue="alnType", sizes=(100,500), size=sizeScatter)
for ax in g.axes.flat:
    ax.tick_params(axis='both', labelleft=True, labelbottom=True)
    ax.axhline(y=75, c='red', linestyle='dashed', label="_horizontal")
    ax.axvline(x=150, c='red', linestyle='dashed', label="_vertical")
    
plt.legend(bbox_to_anchor=(5.0,1), loc=0, borderaxespad=2,ncol=6, shadow=True, labelspacing=1.5, borderpad=1.5)
plt.show()
```


## More detailed analysis

There are two proteins that are bogus. We don't know why, so we just ignore. For the next steps, we also _only_ look at the amino acid hits, not the nucleotide hits.

```
data = pd.read_csv('bigtable.tsv.gz',compression='gzip',header=0,sep='\t')

blacklist = ['A0A097ZRK1', 'G0W2I5']
virusesFiltered = data[(data.alnType == "aa") & (data.kingdom == "Viruses") & (~data.targetID.isin(blacklist) & (data.evalue < 1e-20))]
```

## Patient and Sample Date

Now we create new columns for the patient and the sample date

```
virusesFiltered[['patient', 'date', 'Sputum or BAL']] = virusesFiltered['sampleID'].str.split('_', expand=True)
```

# Merge the real data and the host data

```
virusesFiltHost = pd.merge(virusesFiltered, vmr[['Family', 'Host source']], left_on="family", right_on="Family", how='left')
```

Take a look at this new table. What have we accomplished? What was done? How did it work?

# Correct more errors

There are two families that are incorrectly annotated in this table. Let's just manually correct them:

```
really_bacterial = ['unclassified Caudoviricetes family', 'unclassified Crassvirales family']
virusesFiltHost.loc[(virusesFiltHost['family'].isin(really_bacterial)), 'Host source'] = 'bacteria'
```

# Group-wise data

Now lets look at some data, and put that into a group. We'll also remember what we've looked 

```
to_remove = []
bacterial_viruses = virusesFiltHost[(virusesFiltHost['Host source'] == 'bacteria')]
to_remove.append('bacteria')
```

Can you make a plot of just the bacterial viruses like we did?

## How many reads map per group, and what do we know about them

```
host_source = "archaea"
rds = virusesFiltHost[(virusesFiltHost['Host source'] == host_source)].shape[0]
sps = virusesFiltHost[(virusesFiltHost['Host source'] == host_source)].species.unique()
fams = virusesFiltHost[(virusesFiltHost['Host source'] == host_source)].family.unique()
print(f"There are {rds} reads that map to {host_source} viruses, and they belong to {len(sps)} species ", end="")
if len(sps) < 5:
    spsstr = "; ".join(sps)
    print(f"({spsstr}) ", end="")
print(f"and {len(fams)} families: {fams}")
to_remove.append(host_source)
```

Can you repeat this for each of (`plants` OR `plants (S)`), `algae`, `protists`, `invertebrates`, `fungi`, 

# Filter out things we don't want

Once we remove the above, what's left?

```
virusesFiltHost = virusesFiltHost[~virusesFiltHost['Host source'].isin(to_remove)]
virusesFiltHost['Host source'].unique()
```

# Now look at what we have

```
#filter
virusesGroup = virusesFiltHost.groupby(by=['family','alnlen','pident'], as_index=False).count()

#styling
sizeScatter = 10 * virusesGroup['count']
sns.set_style("darkgrid")
sns.set_palette("colorblind")
sns.set(rc={'figure.figsize':(20,10)})

g = sns.FacetGrid(virusesGroup, col="family", col_wrap=6)
g.map_dataframe(sns.scatterplot, "alnlen", "pident", alpha=.8, sizes=(100,500), size=sizeScatter)
for ax in g.axes.flat:
    ax.tick_params(axis='both', labelleft=True, labelbottom=True)
    ax.axhline(y=80, c='red', linestyle='dashed', label="_horizontal")
    ax.axvline(x=150, c='red', linestyle='dashed', label="_vertical")

g.fig.subplots_adjust(hspace=0.4)
g.set_axis_labels("Alignment length", "Percent Identity")
plt.legend(bbox_to_anchor=(5.0,1), loc=0, borderaxespad=2,ncol=6, shadow=True, labelspacing=1.5, borderpad=1.5)
plt.show()
```

# Now look at each family separately

```
virusesFiltHost[virusesFiltHost.family == 'Retroviridae']
```








