import pandas
import ggplot
import glob
import os

data = {}

for ff in glob.glob('930*'):
    n = os.path.split(ff)[1].split('_')[3]
    tmp = pandas.read_csv(ff)
    for i in tmp.columns:
        tmp[i] = tmp[i] / float(tmp.ix[0][i])
    tmp = tmp.drop(tmp.index[[0]])
    data[n] = tmp

max_val = []
n_val = []
iterations = map(lambda x: str(x),range(2,13))
for iteration in iterations:
    data_i = data[iteration]
    max_val += data_i.max().tolist()
    n_val += [iteration]*data_i.shape[1]

rr = pandas.DataFrame([n_val,max_val]).T
rr.columns = ['Iteration','Max']

ggplot.ggplot(rr, ggplot.aes(x='Iteration', y='Max')) + ggplot.geom_boxplot()


vals = []
n_val = []
iterations = map(lambda x: str(x),range(2,13))
for iteration in iterations:
    data_i = data[iteration]
    vals += (data_i==0).sum().tolist()
    n_val += [iteration]*data_i.shape[1]

rr = pandas.DataFrame([n_val,vals]).T
rr.columns = ['Iteration','Equal to 0']

ggplot.ggplot(rr, ggplot.aes(x='Iteration', y='Equal to 0')) + ggplot.geom_boxplot()


vals = []
n_val = []
iterations = map(lambda x: str(x),range(2,13))
for iteration in iterations:
    data_i = data[iteration]
    vals += data_i.quantile(0.99).tolist()
    n_val += [iteration]*data_i.shape[1]

rr = pandas.DataFrame([n_val,vals]).T
rr.columns = ['Iteration','Median']

ggplot.ggplot(rr, ggplot.aes(x='Iteration', y='Median')) + ggplot.geom_boxplot() + ggplot.ylim(0,0.00025)
