import pandas as pd
from matplotlib import pyplot as plt
from matplotlib import ticker

plt.style.use(['science', 'vibrant', 'grid', 'no-latex']) # 'bright', 'vibrant', 'ieee', 'grid', 'notebook', 'dark_background', 'high-vis'
plt.rcParams['figure.figsize'] = (8, 5)# largo, ancho

df = pd.read_csv('NTimeReg.csv', sep =',')
df.columns = ['Tau', 'FO', 'Time', 'size']
df = df[['Tau', 'FO', 'Time']]
#df.drop(df.index[0], inplace=True)
#df.drop(df.index[-1], inplace=True)

df2 = pd.read_csv('NP_NTimeReg.csv', sep =',')
df2.columns = ['Tau', 'FO', 'TimeNP', 'size']
df2 = df2[['Tau', 'FO', 'TimeNP']]

df['TimeNP'] = df2.TimeNP

fig, ax1 = plt.subplots()

ax1.yaxis.set_major_formatter(ticker.StrMethodFormatter('{x:,.0f}'))

color = '#007FFF'
ax1.set_xlabel(r'Past time periods')
ax1.set_ylabel('OF [USD]') #, color=color
#ax1.plot(df.Tau, df.FO, color=color)
bar1 = ax1.bar(df.Tau, df.FO, width=0.8, color=color, alpha=0.7, label='OF [USD]')
ax1.set_ylim(bottom=126000, top=126600)
ax1.tick_params(axis='y') #, labelcolor=color

# TODO: tick label
ax1.set_xticks(df.Tau)
ax1.set_xticklabels([fr'$t - {i} \leq \tau$' for i in df.Tau], rotation=90)
legend1 = ax1.legend(labels = ['FO [USD]'], loc='lower right', shadow=False)

ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

ax1.yaxis.set_major_formatter(ticker.StrMethodFormatter('{x:,.0f}'))

color = '#FF7F00'
ax2.set_ylabel('Run time [s]') #, color=color  # we already handled the x-label with ax1
lns1 = ax2.plot(df.Tau, df.Time, color=color, lw=1, ls='dashed', marker='v', ms=4, label='Run time with algorithm 1 [s]')
lns2 = ax2.plot(df.Tau, df.TimeNP, color="#FF00AA", lw=1, ls='dashed', marker='*', ms=4, label='Run time without algorithm 1 [s]')
ax2.set_ylim(bottom=0, top=120)
ax2.tick_params(axis='y') #, labelcolor=color

fig.tight_layout()  # otherwise the right y-label is slightly clipped
plt.xticks(rotation=70)

# added these three lines
lns = [bar1, lns1[0], lns2[0]]
labs = [l.get_label() for l in lns]
ax1.legend(lns, labs, loc='lower right')

#legend2 = ax2.legend(labels = ['Time [s]'], loc='center right', shadow=False)

plt.savefig('6bus_FO_Tau.pdf')
plt.close()
#plt.show()
