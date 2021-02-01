import time
import numpy as np
import queue 
import json 
import sys, os
import multiprocessing
from multiprocessing import Process
from pqdict import PQDict
import pickle 
import argparse

import networkx as nx
from tpath import TPath

class Rout:
    #def __init__(self, G, maxsize, node_size, eta, node_list, graph_store_name, time_budget ):
    def __init__(self, subpath, subpath2, graph_store_name, time_budget, filename, fpath_desty, fvedge_desty, maxsize, degree_file, fedge_desty, umatrix_path, upath, process_num, speed_file, eta=12, sigma=3):
        self.maxsize = maxsize
        self.eta = eta
        self.sigma = sigma
        #self.Q = queue.Queue(maxsize=maxsize)
        self.graph_store_name = graph_store_name
        self.time_budget = time_budget 
        self.filename = filename
        self.fpath_desty = fpath_desty
        self.fvedge_desty = fvedge_desty
        self.subpath = subpath
        self.subpath2 = subpath2
        self.degree_file = degree_file
        self.fedge_desty = fedge_desty
        self.umatrix_path = umatrix_path
        self.upath = upath
        self.process_num = process_num
        self.speed_file = speed_file

    def get_speed(self, ):
        speed_dict = {}
        with open(self.speed_file) as fn:
            for line in fn:
                line = line.strip().split('\t')
                speed_dict[line[0]] = {3600*float(line[1])/float(line[2]): 1.0}
        return speed_dict

    def get_graph(self, edge_desty):
        speed_dict = self.get_speed()
        all_nodes, all_edges = set(), set()
        for key in speed_dict:
            line = key.split('-')
            all_nodes.add(line[0])
            all_nodes.add(line[1])
            all_edges.add(key)
        all_nodes, all_edges = list(all_nodes), list(all_edges)
        All_edges = []
        for edge in all_edges:
            edge_ = edge.split('-')
            if edge in edge_desty:
                cost1 = edge_desty[edge].keys()
                cost = min(float(l) for l in cost1)
                All_edges.append((edge_[0], edge_[1], cost))
            elif edge in speed_dict:
                cost1 = speed_dict[edge].keys()
                cost = min(float(l) for l in cost1)
                All_edges.append((edge_[0], edge_[1], cost))
        for edge in speed_dict:
            if edge not in edge_desty:
                edge_desty[edge] = speed_dict[edge]
        G2 = nx.DiGraph()
        G2.add_nodes_from(all_nodes)
        G2.add_weighted_edges_from(All_edges)

        fn = open(self.subpath+self.graph_store_name)
        edges, edges_p, nodes, nodes_p = [], [], set(), set()
        l, l2 = 0, 0
        for line in fn:
            line = line.strip()
            if line not in edge_desty :
                if line in speed_dict:
                    edge_desty[line] = speed_dict[line]
                else:
                    l += 1
                continue
            l2 += 1
            line = line.split('-')

            edges.append((line[0], line[1]))
            nodes.add(line[0])
            nodes.add(line[1])
        fn.close()
        print('%d %d'%(l, l2))
        for key in speed_dict:
            kk = key.strip().split('-')
            nodes.add(kk[0])
            nodes.add(kk[1])
            edges.append((kk[0], kk[1]))
            if key not in edge_desty:
                edge_desty[key] = speed_dict[key]
        #nodes_other = nodes - nodes_p
        nodes = list(nodes)
        nodes_p = list(nodes_p)
        print('len nodes %d'%len(nodes))
        print('len nodes_p %d'%len(nodes_p))
        G = nx.DiGraph()
        G.add_nodes_from(nodes)
        G.add_edges_from(edges)
        return edges, nodes, edge_desty, G2
    
    def get_dict(self, ):
        
        with open(self.subpath+self.fedge_desty) as js_file:
            edge_desty = json.load(js_file)
        return edge_desty

    def get_nodes(self, ):
        gnodes = os.listdir('./res3/u_mul_matrix3/')
        return gnodes


    def compute_one_rowu(self, v_i, v_d, sucpre, edge_desty, nodes_order, U, Q, G, iset, pred):
        st1, vimin, v_s = v_i, 0.0, sucpre[v_i]
        vi_order, vs_order = nodes_order[v_i], nodes_order[v_s]
        while st1 != v_d:
            st2, st1 = st1 , pred[st1]
            tedge = st2 + '-' + st1
            if tedge in edge_desty:
                vimin += min(list([float(l) for l in edge_desty[tedge].keys()]))
            else:
                print('error ... ')
        row_l = vimin
        row_s = self.eta * self.sigma
        for i in range(self.eta):
            if i * self.sigma < row_l:
                U[vi_order][i] = 0.0
        row_x = row_l
        #row_inx = int(row_l / self.sigma)
        while row_x < self.eta * self.sigma:
            i = int(row_x / self.sigma)
            U[vi_order][i] = 0
            for row_c in edge_desty[v_i+'-'+v_s]:
                #inx_xc = int((row_x-float(row_c))/self.sigma) 
                inx_xc =  int((row_x-float(row_c))/self.sigma) 
                if inx_xc < 0 : continue
                U[vi_order][i] = U[vi_order][i] + edge_desty[v_i+'-'+v_s][row_c] * U[vs_order][inx_xc]
            if U[vi_order][i] < 1:
                row_x = row_x + self.sigma
            else:
                break
        row_s = min(row_s, row_x)
        for i in range(int(row_s/self.sigma), self.eta):
            U[vi_order][i] = 1.0 
        return True


    def get_dijkstra3(self, G,  target):
        Gk = G.reverse()
        inf = float('inf')
        D = {target:0}
        Que = PQDict(D)
        P = {}
        nodes = Gk.nodes
        U = set(nodes)
        while U:
            #print('len U %d'%len(U))
            #print('len Q %d'%len(Que))
            if len(Que) == 0: break
            (v, d) = Que.popitem()
            D[v] = d
            U.remove(v)
            #if v == target: break
            neigh = list(Gk.successors(v))
            for u in neigh:
                if u in U:
                    d = D[v] + Gk[v][u]['weight']
                    if d < Que.get(u, inf):
                        Que[u] = d
                        P[u] = v
        return P


    def rout_(self, v_d, edge_desty, edges, nodes, nodes_order, G):
        pred = self.get_dijkstra3(G, v_d)
        N = len(nodes)
        U = np.zeros(N * self.eta).reshape(N, self.eta)
        U.fill(-1)
        Q = []
        iset = set()
        print('pid %d, v d %s' %(os.getpid(), v_d))
        if len(Q) != 0:
            print('error, Q is not empty, exit')
            sys.exit()
        Q.append(v_d)
        iset.add(v_d)
        for j in range(self.eta):
            U[nodes_order[v_d]][j] = 1.0
        lsucess, lorder = '', -1
        flag = True
        sucpre = {}
        while len(Q) != 0:
            v_i = Q.pop(0)
            v_order = nodes_order[v_i]
            #print('v_order %d'%v_order)
            #print('len of predecessors: %d'%(len(list(self.G.predecessors(v_i)))))
            if v_i == v_d:
                U[v_order].fill(1) # = 1
                iset.add(v_i)
            else:
                #self.computeU(v_i, v_order, edge_desty, nodes_order, U, Q, G, iset)
                self.compute_one_rowu(v_i, v_d, sucpre, edge_desty, nodes_order, U, Q, G, iset, pred)
                iset.add(v_i)
        
            #print(list(G.predecessors(v_i)))
            for v in list(G.predecessors(v_i)):
                if U[nodes_order[v]][0] == -1:
                    if v not in Q :
                    #if v not in iset:
                        sucpre[v] = v_i
                        Q.append(v)
        print('pid %d, v_d %s, len iset %d'%(os.getpid(), v_d, len(iset)))
        U = np.round(U, 4)
        #sys.exit()
        return U
    
    def rout(self, sub_nodes, edge_desty, edges, nodes, nodes_order, G):
        print('pid %d'%os.getpid())
        for v_d_ in sub_nodes:
            U_1 = self.rout_(v_d_, edge_desty, edges, nodes, nodes_order, G)
            self.save_matrix(v_d_, U_1, nodes)
    
    def write_json(self, dicts, name):
        with open(self.subpath+name, 'w' ) as fw:
            json.dump(dicts, fw, indent=4)

    def write_file(self, name, ilist):
        strs = '\n'.join(str(ilist_) for ilist_ in ilist)
        fn = open(self.upath+name, 'w')
        fn.write(strs)
        fn.close()

    def save_matrix(self, name, U_2, nodes):
        strs = ''
        for ix, ky in enumerate(nodes):
            U_ = U_2[ix]
            AL = np.argwhere(U_ > 0)
            AH = np.argwhere(U_ > 0.999999)
            if len(AL) == 0:
                al = 0
                strs += ky+';'+str(-1)+';'+str(-1)+'\n'
            else:
                al = AL[0][0]
                if len(AH) == 0:
                    ah = self.eta
                else:
                    ah = AH[0][0]
                strs += ky+';'+str(al-1)+';'+str(ah)+';'+';'.join(str(l) for l in U_[al:ah])+'\n'
            #strs += ky+';'+';'.join(str(l) for l in U_) + '\n'
        fw = open(self.umatrix_path+name, 'w')
        fw.write(strs)
        fw.close()

    def collect_res(self, sub_res):
        self.result.extend(sub_res)

    def main(self,):
        if not os.path.exits(self.subpath2):
            os.mkdir(self.subpath2)
        if not os.path.exits(self.umatrix_path):
            os.mkdir(self.umatrix_path)
        edge_desty = self.get_dict()
        edges, nodes, edge_desty, G = self.get_graph(edge_desty)
        nodes_order, nodes_order_r = {}, {}
        self.result = []
        final_inx = []
        i = 0
        for node in nodes:
            nodes_order[node] = i
            nodes_order_r[i] = node
            i += 1
        N = len(nodes)

        if N % self.process_num == 0:
            t_inx = int(N/self.process_num)
        else:
            t_inx = int(N/self.process_num)+1
        
        pool = multiprocessing.Pool(self.process_num)
        print('begin ... ')
        for len_thr in range(self.process_num):
            mins = min((len_thr+1)*t_inx, N)
            sub_array = nodes[len_thr*t_inx:mins]
            print(len(sub_array))
            pool.apply_async(self.rout, args=(sub_array, edge_desty, edges, nodes,  nodes_order, G))
        pool.close()
        pool.join()
        
        print('len result %d'%len(self.result))
        print('process end ...')

if __name__ == '__main__':
    threads_num = 15
    process_num = 20
    time_budget = 1000
    maxsize = 10000
    #eta, sigma = 333, 30
    #eta, sigma = 111, 90
    #eta, sigma = 170, 60
    #eta, sigma = 800, 10
    parser = argparse.ArgumentParser(description='T-BS')
    parser.add_argument('--sig', default=0, type=int)
    args = parser.parse_args()
    if args.sig == 0:
        sigma, eta = 10, 800
    elif args.sig == 1:
        sigma, eta = 30, 333
    elif args.sig == 2:
        sigma, eta = 60, 170
    elif args.sig == 3:
        sigma, eta = 90, 111
    else:
        print('wrong sig , exit')
        sys.exit()
    print('eta: %d, sigma: %d'%(eta, sigma))

    dinx = 3
    subpath = '/q/storage/yuanye/work/georgi/genvpath/res%d/'%dinx
    subpath2 = './res%d'%dinx
    filename = '/q/storage/yuanye/work/data/AAL_short_%d.csv'%dinx
    fpath_desty = 'KKdesty_num_%d.json'%threads_num
    fvedge_desty = 'KK_vedge_desty2.json'
    fedge_desty = 'M_edge_desty.json'
    graph_store_name = 'KKgraph_%d.txt'%threads_num
    degree_file = 'KKdegree2_%d.json'%threads_num
    umatrix_path = subpath2+'u_mul_matrix_sig%d/'%sigma
    upath = subpath+'u_mul_path_sig%d/'%sigma
    speed_file = '/q/storage/yuanye/work/data/AAL_NGR'
    rout = Rout(subpath, subpath2, graph_store_name, time_budget, filename, fpath_desty, fvedge_desty, maxsize, degree_file, fedge_desty, umatrix_path, upath, process_num, speed_file, eta, sigma)

    rout.main()

