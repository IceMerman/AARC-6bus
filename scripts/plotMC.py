#!pip install SciencePlots

import pandas as pd
from matplotlib import pyplot as plt

plt.style.use(['science', 'high-vis', 'grid', 'no-latex']) # 'bright', 'vibrant', 'ieee', 'grid', 'notebook', 'dark_background', 'high-vis'
plt.rcParams['figure.figsize'] = (10, 3)# largo, ancho
#plt.rcParams['legend.loc'] = 'upper left'

df = pd.read_csv('MCreg.csv', sep =',')

df.columns = ['taus' if idx==0 else val for idx, val in enumerate(df.columns)]

df.index = df.taus
df=df.iloc[:,1:]
df=df.transpose()
df=df*301

taus = [i-1 for i in [1,5,10,15,20,24]] #list(range(25))
latexLegends = [fr'$t - {i+1} \leq \tau$' for i in taus]
datalen = len(df.index)
ax = df.iloc[:,taus].plot.line()

fig = plt.gcf()
#fig.set_size_inches(3, 4)
fig.set_label('MonteCarlo simulations')

#ax = plt.gca()
ax.set_xlabel('Relative uncertainty level [%]')
ax.set_ylabel('Number of violated constraints')

customXticks = [idx for idx, val in enumerate(df.index) if (idx+1)%3==0]
ax.set_xticks(customXticks)
customXlabels = ['{:.0f}'.format(float(i)*100) for i in df.index[customXticks]]
ax.set_xticklabels(customXlabels, rotation=90)

legend = ax.legend(labels = latexLegends, loc='upper left', shadow=False)
#plt.savefig('MC_fullScale.svg')
plt.savefig('6busMC_fullScale.svg')
#plt.show()

# Zoom

datalen = len(df.index)

#             center    left values : center right values
df2 = df.iloc[datalen//2 - 4 : datalen//2 + 10, taus]

ax = df2.plot.line()

fig = plt.gcf()
#fig.set_size_inches(3, 4)
fig.set_label('MonteCarlo simulations')

#ax = plt.gca()
ax.set_xlabel('Relative uncertainty level [%]')
ax.set_ylabel('Number of violated constraints')

customXticks = [idx for idx, val in enumerate(df2.index) if (idx+1)%2==0]
ax.set_xticks(customXticks)
customXlabels = ['{:.0f}'.format(float(i)*100) for i in df2.index[customXticks]]
ax.set_xticklabels(customXlabels, rotation=90)

legend = ax.legend(labels = latexLegends, loc='upper left', shadow=False)
#plt.savefig('MC_Zoom.svg')
plt.savefig('6busMC_Zoom.svg')
#plt.show()
