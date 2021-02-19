 %mem=200MW                         
 %chk=tmp
 %NProcShared=          12
 # m062x/6-311+g(d,p) irc(reverse,calcall,stepsize=3,maxpoints=10)      
 # int=ultrafine  nosym  iop(7/33=1)                                    

 geom            0

           0           2
   c1                                                        
   c2  c1 rcc1                                               
   h3  c1 rch1 c2 ahcc1                                      
   h4  c1 rch2 c2 ahcc2 h3 b1                                
   h5  c2 rch3 c1 ahcc3 h3 b2                                
   h6  c2 rch4 c1 ahcc4 h5 b3                                
   H1          2 rts    1 aabs1    3 babs1                   

 RCC1                             1.3386000000000000     
 RCH1                             1.0833999999999999     
 RCH2                             1.0833999999999999     
 RCH3                             1.0839000000000001     
 RCH4                             1.0839000000000001     
 AHCC1                            121.45780000000001     
 AHCC2                            121.45750000000000     
 AHCC3                            121.34860000000000     
 AHCC4                            121.35050000000000     
 B1                               181.92560000000000     
 B2                              -5.5297999999999998     
 B3                               189.14410000000001     
 AABS1                            106.58499999999999     
 BABS1                            89.077200000000005     
 RTS                              1.9963000000000000     

