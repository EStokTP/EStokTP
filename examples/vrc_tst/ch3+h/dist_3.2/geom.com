 %mem=400MW                         
 %chk=tmp
 %NProcShared=          32
 #  uwb97xd/aug-cc-pvtz opt=(noeig,intern,maxcyc=50)                    
 #  int=ultrafine nosym  iop(7/33=1) guess=(mix,always)                 

 geom            0

           0           1
 C1                                                          
 H2   C1   R1                                                
 H3   C1   R2   H2   A2                                      
 H4   C1   R3   H2   A3   H3   D3                            
  H21         1 rts    2 aabs1    3 babs1                    

 BABS1                            90.947900000000004     
 AABS1                            90.486999999999995     

 R3                               1.0782000000000000     
 A2                               119.97210000000000     
 A3                               120.02070000000001     
 D3                               179.80760000000001     
 R2                               1.0782000000000000     
 R1                               1.0782000000000000     
 RTS                              3.2000000000000002     

