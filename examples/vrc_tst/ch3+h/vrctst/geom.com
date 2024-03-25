 %mem=400MW                         
 %chk=tmp
 %NProcShared=          32
 # uwb97xd/aug-cc-pvtz opt=(noeig,intern,maxcyc=50)                     
 # int=ultrafine nosym freq iop(7/33=1) guess=(mix,always)              

 geom            0

           0           1
 C1                                                          
 H2   C1   R1                                                
 H3   C1   R2   H2   A2                                      
 H4   C1   R3   H2   A3   H3   D3                            
 H21 C1 BB H2 ANG1 H3 DH1                                    

 R1                               1.0821000000000001     
 R2                               1.0821000000000001     
 R3                               1.0820000000000001     
 A2                               117.02900000000000     
 A3                               117.05200000000001     
 D3                               213.56530000000001     
 DH1                              106.71920000000000     
 ANG1                             100.03080000000000     

 BB                               2.0000000000000000     

