# route-sota

## generate 
This folder includes the generation scripts for virtual path and matrix U 

### virtual path generation
4 scripts needed to generate the virtual path. 

```
python get_tpath.py
python get_overlap.py
python get_vpath.py
python get_vpath_desty.py
```

### matrix U (heuristic table) generation 
1 script needed to generate the matrix U 

```
python get_policy_U.py 
```

## routing
This folder includes 7 routing algorithms implemented in paper. 

```
# Binary heuristic routing using Euclidean distance 
./run_T-B-EU.sh

# Binary heuristic routing using shortest path tree 
./run_T-B-E.sh 

# Binary heuristic routing using shortest path tree and T-paths 
./run_T-B-P.sh 

# Binary heuristic routing using matrix (heuristic table) U
./run_T-BS.sh 

# Routing only with virtual paths, no any heuristic 
./run_V-None.sh 

# Routing with virtual paths and matrix (heuristic table) U 
./run_V-BS.sh  
```

## data
The data used in the experiment. 


