#NDIM
2
#PRIM
1.0 0.0
0.0 1.0
#SUPER
4 0
0 4
#ORB
s0 0.0d0 0.0d0 0.0d0 #0
s1 0.5d0 0.0d0 0.0d0 #1
s2 0.0d0 0.5d0 0.0d0 #2
#HAMILT 
0 1 0.5 0.0 0.0 1.0 1.0 0.0
0 1 -0.5 0.0 0.0 1.0 1.0 0.0
0 2 0.0 0.5 0.0 1.0 1.0 0.0
0 2 0.0 -0.5 0.0 1.0 1.0 0.0
0 0 0.0 0.0 0.0 0.0 0.0 0.0
1 1 0.0 0.0 0.0 0.0 0.0 0.0
2 2 0.0 0.0 0.0 0.0 0.0 0.0
#SYMM
d  0.0d0 0.0d0 0.0d0 1.0d0 0.0d0 0.0d0
d  0.0d0 0.0d0 0.0d0 0.0d0 1.0d0 0.0d0
c4 0.0d0 0.0d0 0.0d0 0.0d0 0.0d0 1.0d0
#PHASE
1 0
0 1
s0 0.0 0.0 0.0  1.0
s1 0.5 0.0 0.0  -1.0
s2 0.0 0.5 0.0  -1.0
#BONDS   
0 0  0.0  0.0  0.0  # 1
1 1  0.0  0.0  0.0  # 2
2 2  0.0  0.0  0.0  # 3
0 1  0.5  0.0  0.0  # 4
0 1  -0.5  0.0  0.0 # 5
0 2  0.0  0.5  0.0  # 6
0 2  0.0  -0.5  0.0 # 7
#PAIR    
            1    2    3    4    5    6    7
s-wave     1.0  1.0  1.0  0.0  0.0  0.0  0.0
#END