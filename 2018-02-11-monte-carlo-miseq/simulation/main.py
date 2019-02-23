import pandas
import glob
import os

data = {}

for ff in glob.glob('*csv'):
    filename = os.path.split(ff)[1]
    if filename.find('sampling.csv')!=-1:
        continue
    params = filename.split('_')

    tmp = pandas.read_csv(ff)
    for i in tmp.columns:
        tmp[i] = tmp[i] / float(tmp.ix[0][i])
    tmp = tmp.drop(tmp.index[[0]])
    data[filename] = {
        'params':{
            'initPop':params[0],
            'polymerase_error_truseq':params[1],
            'polymerase_error_replig':params[2],
            'cycles_replig':params[3],
            'cycles_truseq':params[4],
            'pcr_efficiency_truseq':params[5],
            'pcr_efficiency_replig':params[6],
            'replig_sampling':params[7].replace('.csv','')
        },
        'data':tmp
    }

d1 = []
d2 = []
d3 = []
d4 = []
d5 = []
for filename, data_i in data.items():
    params = data_i['params']
    data_i = data_i['data']
    d1 += data_i.max().tolist()
    d2 += [params['initPop']]*data_i.shape[1]
    d3 += [params['replig_sampling']]*data_i.shape[1]
    d4 += [params['pcr_efficiency_truseq']]*data_i.shape[1]
    d5 += [params['polymerase_error_replig']]*data_i.shape[1]

rr = pandas.DataFrame([d1,d2,d3,d4,d5]).T
rr.columns = ['Max','InitialPop','Sampling',"PCREfficiency","PolymeraseError"]
rr.to_csv('by_initPop_sampling.csv')
