 %mem=400MW                         
 %chk=tmp
 %NProcShared=          16
 #UWB97XD/JUN-CC-PVTZ
 #freq=readfc iop(7/33=1)                                  guess=read geom=check int=ultra

 geom            0

           0           2
                                                                                 
                                                                                 
 H3   C1   R2   O2   A2                                                          
 H4   C1   R3   O2   A3   H3   D3                                                
 H5   O2   R4   C1   A4   H3   D4                                                

 R1                               1.3575999999999999     
 R2                               1.0779000000000001     
 R3                               1.0817000000000001     
 R4                              0.95730000000000004     
 A2                               113.96800000000000     
 A3                               119.04490000000000     
 A4                               109.87820000000001     
 D3                               205.75290000000001     
 D4                               175.17810000000000     

