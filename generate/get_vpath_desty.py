from tpath import TPath
import json 
import numpy as np
from v_opt import VOpt
import operator


dinx = 50
fname = '/q/storage/yuanye/work/data/AAL_short_%d.csv'%dinx
B = 3

def norms(dicts2, is_norm):
    for e_key in dicts2:
        sums = sum(dicts2[e_key].values())
        for time in dicts2[e_key]:
            dicts2[e_key][time] = round(100.0*dicts2[e_key][time]/sums, 6)
        if is_norm:
            dicts2[e_key] =  dict(sorted(dicts2[e_key].items(), key=operator.itemgetter(1), reverse=True))
    return dicts2

def get_vopt(vopt, dicts, is_path=False, is_norm=False):
    for e_key in dicts:
        time_freq = dicts[e_key] 
        time_cost, time_freq = [float(l) for l in time_freq.keys()], [float(l) for l in time_freq.values()]
        pvalues = {}
        if len(time_cost) > 1:
            pvalues = vopt.v_opt(time_cost, time_freq, B) # pvalues is a dict, contains new key and value
            dicts[e_key] = pvalues
        else:
            if not is_path:
                k1 = str(int(np.mean(time_cost)))
                pvalues[k1] = 1.0
    return norms(dicts, is_norm)


tpath = TPath(fname)
data = tpath.load()
edge_desty = tpath.get_edge_time_freq(data)

vopt = VOpt()

edge_desty_ = get_vopt(vopt, edge_desty, False, True)


keys = list(edge_desty_.keys())
print(keys[0])
print(edge_desty_[keys[0]])
edge_desty_1 = {}
for edge in edge_desty_:
    #ky = list(edge.keys())
    edge_ = {}
    ky = edge_desty_[edge]
    for k in ky.keys():
        edge_[str(np.abs(k))] = float(ky[k])/100
    edge_desty_1[edge] = edge_

#ffname = '../res3/M_vedge_desty_10_2.json'
ffname = './res%d/M_edge_desty.json'%dinx
with open(ffname, 'w') as fw:
    json.dump(edge_desty_1, fw, indent=4)
