import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys, os, gc
import threading
import multiprocessing 

from multiprocessing import Process

from v_opt import VOpt
import json
import operator
import networkx as nx
import pickle

class Overlp():
    def __init__(self, filename, subpath, p_path, subpath_range):
        self.filename = filename
        self.subpath = subpath
        self.p_path = p_path
        self.subpath_range = subpath_range
    
    def cal_overlap(self, path_keys):
        path_keys_ = [set([int(c) for c in key.split('-')]) for key in path_keys]
        #path_keys_ = [set(key.split('-')) for key in path_keys]
        N = len(path_keys_)
        print(N)
        overlap = {}
        for i in range(N-1):
            if i% 10 == 0:
                print(i)
            for j in range(i+1, N):
                intes = path_keys_[i].intersection( path_keys_[j])
                #intes = np.intersect1d(path_keys_[i], path_keys_[j],assume_unique=False)
                if len(intes) > 1 and len(intes) < min(len(path_keys_[i]), len(path_keys_[j])):
                    if i not in overlap:
                        overlap[i] = []
                    if j not in overlap:
                        overlap[j] = []
                    overlap[i].append(j)
                    overlap[j].append(i)
        return overlap

    def norms(self, dicts2, is_norm):
        for e_key in dicts2:
            sums = sum(dicts2[e_key].values())
            for time in dicts2[e_key]:
                dicts2[e_key][time] = round(100.0*dicts2[e_key][time]/sums, 6)
            if is_norm:
                dicts2[e_key] =  dict(sorted(dicts2[e_key].items(), key=operator.itemgetter(1), reverse=True))
        return dicts2

    def get_vopt(self, vopt, dicts, is_path=False, is_norm=False):
        for e_key in dicts:
            time_freq = dicts[e_key] 
            time_cost, time_freq = [float(l) for l in time_freq.keys()], [float(l) for l in time_freq.values()]
            pvalues = {}
            if len(time_cost) > 1:
                pvalues = vopt.v_opt(time_cost, time_freq, self.B) # pvalues is a dict, contains new key and value
                dicts[e_key] = pvalues
            else:
                if not is_path:
                    k1 = str(int(np.mean(time_cost)))
                    pvalues[k1] = 1.0
        return self.norms(dicts, is_norm)

    def get_inters2(self, path_1, path_2, flag=False):
        #print('path 1 2')
        #print(path_1)
        #print(path_2)
        end = path_1.find(';')
        #end = path_1.find('-', end+1)
        lhead = path_1[:end]
        end2 = path_2.find(lhead)
        if end2 == -1:
            if flag: return 0
            return self.get_inters2(path_2, path_1, True)
        else:
            b = min(len(path_2)-end2, len(path_1))
            a = path_2[end2:] == path_1[:b]
            #print('sub path 1 2')
            #print(path_2[end2:])
            #print(path_1[:b])
            if not a: return 0
            else: return len(path_2[end2:].split(';'))

    def special2(self, path_keys, path_count):
        rev_dict = {}
        N = len(path_keys)
        print(N)
        M = 300
        for i in range(N):            
            key = path_keys[i]
            key = set(key.split(';'))
            for k_ in key:
                if k_ not in rev_dict:
                    rev_dict[k_] = []
                rev_dict[k_].append(i)
        overlap, path_count_ = {}, {}
        for i in range(N):
            overlap[i] = [0] * M
        gc.disable()
        for i in range(N):
            if i%100 == 0: print(i)
            #if i%10000 == 0: 
            #    gc.enable()
            #    gc.disable()
            key = path_keys[i]
            key2 = set(key.split(';'))
            A = {}
            for k_ in key2:
                valu = rev_dict[k_]
                for va in valu: # valu is a int list
                    if va not in A:
                        A[va] = 1
            for va in A:
                if va != i:
                    pat_ = path_count[va]
                    path_len = len(pat_.split(';'))
                    intes = self.get_inters2(key, pat_)
                    #print(intes)
                    if intes > 0 and intes<min(len(key2), path_len):
                        overlap[i][0] += 1
                        if overlap[i][0] >= M and overlap[i][0]%M == 0:
                            overlap[i].extend([0]*M)
                        overlap[i][overlap[i][0]] = va
            #sys.exit()
        gc.enable()
        for k in overlap:
            overlap[k] = overlap[k][1:overlap[k][0]+1]
        return overlap

    def get_k(self, key):
        end, j = -1, 0
        key = key.split('-')
        key_ = key[0]
        for kk in key[1:]:
            if j % 2 == 0:
                key_ += '-'
            else:
                key_ += ';'
            key_ += kk
            j += 1
        return key_
    
    def get_subpath(self, ):
        path_desty, path_count = {}, {}
        iters = 0
        for i in self.subpath_range:
            with open(self.subpath+self.p_path+str(i)+'.json') as json_file:
                jsdata = json.load(json_file)
                for key in jsdata.keys():
                    #print(key)
                    key_ = self.get_k(key)
                    skey = key_.split(';')
                    if len(skey) > len(set(skey)):
                        continue
                    #print(key_)
                    path_desty[key_] = jsdata[key]
                    path_count[iters] = key_
                    iters += 1

        overlap = self.special2(list(path_desty.keys()), path_count)
        return path_desty, path_count, overlap

    def write_json(self, js_dict, fname):
        #json.dumps(js_dict, fname)
        with open(fname, 'w') as fw:
            json.dump(js_dict, fw, indent=4)
    def write_file(self, js_dict, fname):
        fn = open(fname, 'w')
        for key in js_dict:
            val=','.join(str(l) for l in js_dict[key])
            fn.write(str(key) +':' + val + '\n')
        fn.close()

    def main(self, ):
        path_desty, path_count, overlap = self.get_subpath()
        print('begin to write file ...')
        self.write_json(path_desty, self.subpath+'path_desty2.json')
        self.write_json(path_count, self.subpath+'path_count2.json')
        #self.write_json(overlap, self.subpath+'overlap.json')
        self.write_file(overlap, self.subpath+'overlap2.txt')
        

if __name__ == '__main__':
    dinx = 50
    filename = '/q/storage/yuanye/work/data/trips_real_short_%d.csv'%dinx
    subpath = './subpath/path_travel_time_'
    subpath = './res%d/'%dinx
    p_path = 'path_travel_time_'
    path_desty_name = 'path_desty2.json'
    path_count_name = 'path_count2.json'
    subpath_range = [l for l in range(2, 39)]
    subpath_range.append(28)
    overlp = Overlp(filename, subpath, p_path, subpath_range)
    overlp.main()


