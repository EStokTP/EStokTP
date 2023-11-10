#!/home/chimica2/abaggioli/bin/python3

import os,subprocess,time

spec = None
log = open( "output/hl_composite.log" , "w" )

def Msg( msg , leave ):
  if leave<2:
    if msg[-1]==".":
      log.write( "%s\n" % ( msg ) )
    else:
      log.write( "%s%s.\n" % ( msg + " at " , time.asctime() ) )
    log.flush()
  if leave>0:
    log.write( "%s%s.\n" % ( "Termination with error at " , time.asctime() ) )
    log.close()
    if spec is not None:
      rename = comm( [ "mv" , "output/hl_composite.log" , "output/" + fname + ".log" ] )
    exit()

def comm( COMM ):
  comm = subprocess.run( COMM , stdout=subprocess.PIPE , stderr=subprocess.PIPE , universal_newlines=True )
  comm
  return ( comm.stdout , comm.stderr )

def EQSplit( K , S ):
  for l in K:
    for s in S:
      EQ[l] = EQ[l].replace( s , " " + s + " " if s!="**" else " $$ " )
    EQ[l] = EQ[l].replace( "$$" , "**" )
    EQ[l] = EQ[l].split()

def EQJoin( K ):
  for l in K:
    EQ[l] = "".join( EQ[l] )

def EQReplace( K , S ):
  EQSplit( K , S )
  for l in K:
    for i,I in enumerate(EQ[l]):
      if I not in S:
        try:
          float( I )
        except:
          if I in K and len(EQ[I])<3:
            try:
              float( "".join(EQ[I]) )
              EQ[l][i] = "".join(EQ[I])
            except: continue
  EQJoin( K )

def EQCheck( K , S ):
  Z = [ True for i in range(len(K)) ]
  for i,(j,J) in enumerate(EQ.items()):
    try:
      float( J )
      Z[i] = False
    except:
      d = 0
      for s in S:
        J = J.replace( s , " " )
      for D in J.split():
        try:
          float( D )
        except:
          d = d + 1
      if d>0: Z[i] = False
  return Z

def EQSolve( a , S ):
  s = 0
  while s<len(S):
    if a[0]=="+": a = a[1:]
    A = a[1:] if a[0]=="-" else a
    if S[s] in A:
      b = a.replace("*"," * ").replace(" *  * "," ** ").replace("/"," / ").replace("+"," + ").replace("-"," - ")
      B = b.split()
      if len(B)==1: break
      for i,I in enumerate(B):
        if i==0: continue
        if I==S[s]:
          if s==0:
            B[i] = "*"
            if B[i+1]=="-":
              B[i+1] = str( 1. / float( "".join( B[i+1:i+3] ) ) )
              del B[i+2]
            else:
              B[i+1] = str( 1. / float( B[i+1] ) )
          else:
            ant = 2 if i>1 and s in (3,4) and B[i-2] in S[3:5] else 1
            des = 3 if B[i+1] in S[3:5] else 2
            old = "".join( B[i-ant:i+des] )
            try: eva = eval( old )
            except: Msg( "  ERROR - could not solve for 'E', possible invalid expression in the $SOLVE input section." , 1 )
            if eva==0: new = ""
            else: new = "+" + "%.24f"%eva if eva>0 and i>1 and s in (3,4) and B[i-2] in S[3:5] else "%.24f"%eva
            a = a.replace( old , new )
            break
      if s==0:
        s = 1
        a = "".join( B )
    else: s = s + 1
  return a

### import files ###

Msg( "Beginning HL Composite calculation" , 0 )

if os.path.isfile( "hl_composite.xyz" ):
  with open( "hl_composite.xyz" ) as f:
    G = [ ln for ln in f ]
  spec = G[1].split()[0]
  chrg = G[1].split()[1]
  mult = G[1].split()[2]
  fname = spec + "_hlcomp"
  Msg( "  'hl_composite.xyz' imported successfully." , 0 )
else:
  Msg( "  ERROR - 'hl_composite.xyz' not found." , 1 )

if os.path.isfile( "output/%s.restart" % ( fname ) ):
  with open( "output/%s.restart" % ( fname ) ) as f:
    L = { ln.split()[0] : ln.split()[1] for ln in f }
  Msg( "  'output/%s.restart' imported successfully." % ( fname ) , 0 )
else:
  L = []

if os.path.isfile( "hl_composite.dat" ):
  with open( "hl_composite.dat" ) as f:
    F = [ ln for ln in f ]
  Msg( "  'hl_composite.dat' imported successfully." , 0 )
else:
  Msg( "  ERROR - 'hl_composite.dat' not found." , 1 )

### parse input and fill EQ and EN ###

n = 0
EN = []
EQ = {}
PA = {}
while n!=len(F):
  if F[n].split() and F[n].split()[0].upper()=="$ENERGY":
    EN.append( F[n].upper().split()[1:3] )
    EN[-1][0] = EN[-1][0].split(",")
    EN[-1].append( [] )
    EQ[ EN[-1][1] ] = "0.0"
    n = n + 1
    while True:
      EN[-1][-1].append( F[n] )
      n = n + 1
      if F[n].split() and F[n].split()[0].upper()=="$END3":
        break
    EN[-1][-1].append( "$END3" )
  elif F[n].split() and F[n].split()[0].upper()=="$PARAM":
    PA[ F[n].split()[1].upper() ] = [ x for x in F[n].upper().split()[2:4] ]
  elif F[n].split() and F[n].split()[0].upper()=="$SOLVE":
    n = n + 1
    while True:
      eq = "".join( F[n].split() ).split( "=" )
      if len(eq)==2:
        EQ[eq[0].upper()] = eq[1].upper()
      else:
        Msg( "    ERROR - invalid expression for '%s', too many '='." % ( eq[0] ) , 1 )
      n = n + 1
      if F[n].split() and F[n].split()[0].upper()=="$END":
        break
  n = n + 1

### check that input is correct ###

pErr = []
PM = { "G09"    : ["KB","MB","GB","TB","KW","MW","GW","TW"] ,\
       "G16"    : ["KB","MB","GB","TB","KW","MW","GW","TW"] ,\
       "MOLPRO" : ["K","M","G"] ,\
       "MRCC"   : ["MB","GB"] }
K = tuple( EQ )
S = ["/","**","*","+","-","(",")"]
EQSplit( K , S )
if "E" not in K:
  Msg( "    ERROR - missing expression for total energy 'E'." , 1 )
for i,en in enumerate(EN):
  if en[0][0] not in ["G09","G16","MOLPRO","MRCC"]:
    Msg( "    ERROR - software '%s' is not supported." % ( en[0][0] ) , 1 )
  try:
    float( en[1] )
  except:
    for nen in en[1]:
      if nen not in "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ1234567890_":
        Msg( "    ERROR - invalid energy name '%s', underscore '_' is the only allowed symbol." % ( en[1] ) , 1 )
  else:
    Msg( "    ERROR - invalid energy name '%s', cannot be a number." % ( en[1] ) , 1 )
  endj = 0
  for ln in en[2]:
    if ln.split() and ln.split()[0].upper()=="$END"+str(endj):
      endj = endj + 1
  if endj<3:
    Msg( "    ERROR - missing or misplaced '$END%d' statement at step %d." % ( endj , i+1 ) , 1 )
  if len(en[0])==1:
    if en[0][0] in PA:
      en[0].append( PA[ en[0][0] ][0] )
      en[0].append( PA[ en[0][0] ][1] )
    else:
      Msg( "    ERROR - no parameters are provided for %s at step %d." % ( en[0][0] , i+1 ) , 1 )
  try:
    int( en[0][1] )
  except:
    pErr.append( [ 0 , en[0][1] , en[0][0] , i+1 ] )
  if en[0][2].upper()[(0-len(PM[en[0][0]][0])):] not in PM[en[0][0]]:
    pErr.append( [ 1 , en[0][2] , en[0][0] , i+1 ] )
  else:
    try:
      int( en[0][2][:(0-len(PM[en[0][0]][0]))] )
    except:
      pErr,append( [ 1 , en[0][2] , en[0][0] , i+1 ] )
if len(pErr)>0:
  for pe in pErr:
    pErrm = [ "    ERROR - invalid number of processors '%s' for software %s at step %d." % tuple( pe[1:] ) ,\
              "    ERROR - invalid memory declaration '%s' for software %s at step %d."   % tuple( pe[1:] ) ]
    Msg( pErrm[ pe[0] ] , 0 )
  Msg( "" , 2 )
dupl = [ en[1] for en in EN ]
if len(dupl)!=len(set(dupl)):
  Msg( "    ERROR - energy names must be unique." , 1 )
for l in K:
  for i,I in enumerate(EQ[l]):
    if I not in S:
      try: float( I )
      except:
        if I not in K:
          Msg( "    ERROR - unknown energy '%s' used in '%s'." % ( I , l ) , 1 )
    if I in S[:-2] and EQ[l][i+1] in S[:-2]:
      Msg( "    ERROR - invalid expression in '%s', consetutive arithmetic operators not allowed, use parentheses." % ( l ) , 1 )
EQJoin( K )
for l in K:
  if "()" in EQ[l]:
    Msg( "    ERROR - invalid use of parentheses in '%s'." % ( l ) , 1 )
  npo = 0
  npc = 0
  for X in EQ[l]:
    if   X=="(": npo = npo + 1
    elif X==")": npc = npc + 1
  if npo!=npc:
    Msg( "    ERROR - unmatched parentheses in '%s'." % ( l ) , 1 )
Msg( "  'hl_composite.dat' parsed successfully." , 0 )

### produce energies in order as they appear in input ###

fail = []
restart = open( "output/" + fname + ".restart" , "w" )
for i,en in enumerate(EN):
  zi = str(i+1).zfill(2)
  Msg( "  Step %s..." % ( str(i+1).zfill(2) ) , 0 )
  # check if results from a previous run exist
  if en[1] in L:
    EQ[ en[1] ] = L[ en[1] ]
    Msg( "    '%s' recovered from 'output/%s.restart'." % ( en[1] , fname ) , 0 )
    Msg( "    Setting '%s' equal to %s." % ( en[1] , L[ en[1] ] ) , 0 )
    restart.write( "%s%30s\n" % ( en[1] , EQ[ en[1] ] ) )
    restart.flush()
  else:
    # concoct input files
    com = []
    n = 0
    T = 0 if int(mult)==1 else 1
    out = open( fname + zi + ".inp" , "w" )
    if en[0][0]=="G09" or en[0][0]=="G16":
      out.write( "%" + "Chk=%s\n" % ( fname + zi ) )
      out.write( "%" + "NProcShared=%s\n" % ( en[0][1] ) )
      out.write( "%" + "Mem=%s\n" % ( en[0][2] ) )
      if T:
        while True:
          if en[2][n].split() and en[2][n].split()[0].upper()=="$END0":
            break
          n = n + 1
        n = n + 1
      while True:
        if en[2][n].split() and en[2][n].split()[0].upper()=="$END"+str(T):
          break
        if "%"+"MEM" not in en[2][n].upper() and "%"+"NPROCSHARED" not in en[2][n].upper() and "%"+"CPU" not in en[2][n].upper() and "%"+"CHK" not in en[2][n].upper():
          out.write( en[2][n] )
        n = n + 1
      out.write( "\ncomposite_hl_step_%s %s\n\n%3s%3s\n" % ( str(i+1).zfill(2) , en[1] , chrg , mult ) )
      for g in G[2:]:
        out.write( g )
      while True:
        if en[2][n].split() and en[2][n].split()[0].upper()=="$END"+str(T+1):
          break
        n = n + 1
      n = n + 1
      while True:
        if en[2][n].split() and en[2][n].split()[0].upper()=="$END"+str(T+2):
          break
        out.write( en[2][n] )
        n = n + 1
      out.write( "\n" )
    elif en[0][0]=="MOLPRO":
      out.write( "memory,%s,%s\n" % ( en[0][2][:-1] , en[0][2][-1:] ) )
      if T:
        while True:
          if en[2][n].split() and en[2][n].split()[0].upper()=="$END0":
            break
          n = n + 1
        n = n + 1
      while True:
        if en[2][n].split() and en[2][n].split()[0].upper()=="$END"+str(T):
          break
        if "MEMORY" not in en[2][n].upper():
          out.write( en[2][n] )
        n = n + 1
      out.write( "geometry={angstrom\n" )
      for g in G[2:-1]:
        for y,Y in enumerate(g.split()):
          out.write( ",%14s" % ( Y ) if y>0 else Y )
        out.write( "\n" )
      out.write( "}\n\nset,spin=%d\n\n! composite_hl_step_%s %s\n\n" % ( int(mult)-1 , str(i+1).zfill(2) , en[1] ) )
      while True:
        if en[2][n].split() and en[2][n].split()[0].upper()=="$END"+str(T+1):
          break
        n = n + 1
      n = n + 1
      while True:
        if en[2][n].split() and en[2][n].split()[0].upper()=="$END"+str(T+2):
          break
        out.write( en[2][n] )
        n = n + 1
      out.write( "\nput,molden,%s.molden\nhl_out=energy\n---\n" % ( fname + zi ) )
    elif en[0][0]=="MRCC":
      out.write( "# composite_hl_step_%s %s\ncharge=%s\nmult=%s\n" % ( str(i+1).zfill(2) , en[1] , chrg , mult ) )
      out.write( "mem=%s\n" % ( en[0][2] ) )
      if T:
        while True:
          if en[2][n].split() and en[2][n].split()[0].upper()=="$END0":
            break
          n = n + 1
        n = n + 1
      while True:
        if en[2][n].split() and en[2][n].split()[0].upper()=="$END"+str(T):
          break
        if "MEM=" not in en[2][n].upper() and "UNIT=" not in en[2][n].upper() and "GEOM=" not in en[2][n].upper():
          out.write( en[2][n] )
        n = n + 1
      out.write( "unit=angstrom\ngeom=xyz\n" )
      for g in G:
        out.write( g )
      while True:
        if en[2][n].split() and en[2][n].split()[0].upper()=="$END"+str(T+1):
          break
        n = n + 1
      n = n + 1
      while True:
        if en[2][n].split() and en[2][n].split()[0].upper()=="$END"+str(T+2):
          break
        out.write( en[2][n] )
        n = n + 1
      out.write( "\n" )
    out.close()
    Msg( "    Input file generated successfully." , 0 )
    # run job
    if en[0][0]=="G09" or en[0][0]=="G16":
      runCheck = [ "which" , en[0][0].lower() ]
      runJob = [ [ en[0][0].lower() , fname + zi + ".inp" ] ]
    elif en[0][0]=="MOLPRO":
      runCheck = [ "which" , "molprop" ]
      runJob = [ [ "cp" , fname + zi + ".inp" , "molpro.inp" ] ,\
                 [ "molprop" , "-n" , en[0][1] , "molpro.inp" ] ,\
                 [ "rm" , "molpro.inp" ] ,\
                 [ "mv" , "molpro.out" , fname + zi + ".out" ] ,\
                 [ "mv" , "molpro.xml" , fname + zi + ".xml" ] ]
    elif en[0][0]=="MRCC":
      runCheck = [ "which" , "mrccp" ]
      runJob = [ [ "mrccp" , "-n" , en[0][1] , fname + zi + ".inp" ] ,\
                 [ "mv" , "mrcc.log" , fname + zi + ".log" ] ]
    if comm( runCheck )[0]=="":
      Msg( "    WARNING - cannot run %s, moving on." % ( en[0][0] ) , 0 )
      fail.append( en[1] )
      continue
    else:
      Msg( "    Starting %s" % ( en[0][0] ) , 0 )
      for rj in runJob: 
        com.append( comm( rj ) )
      Msg( "    Calculation ended" , 0 )
    # read output
    if en[0][0]=="G09" or en[0][0]=="G16":
      if comm( [ "grep" , "Normal termination" , fname + zi + ".log" ] )[0]=="":
        Msg( "    WARNING - no normal termination occurred, moving on." , 0 )
        fail.append( en[1] )
        continue
      if comm( [ "which" , "formchk" ] )[0]!="":
        fchk = "formchk"
      elif comm( [ "which" , "formchk" + en[0][0][1:] ] )[0]!="":
        fchk = "formchk" + en[0][0][1:]
      else:
        Msg( "    WARNING - cannot format checkpoint file, moving on." , 0 )
        fail.append( en[1] )
        continue
      com.append( comm( [ fchk , fname + zi + ".chk" ] ) )
      com.append( comm( [ "grep" , "Total Energy" , fname + zi + ".fchk" ] ) )
      if com[-1][0]=="":
        Msg( "    WARNING - no total energy found in '%s', moving on." % ( fname + zi + ".fchk" ) , 0 )
        fail.append( en[1] )
        continue
      else:
        EQ[ en[1] ] = str( float( com[-1][0].split()[-1] ) )
        Msg( "    Setting '%s' equal to %s." % ( en[1] , EQ[ en[1] ] ) , 0 )
        com.append( comm( [ "cp" , "%s.log" % ( fname + zi ) , "hl_logs/%s.log" % ( fname + zi ) ] ) )
        Msg( "    Output file stored as 'hl_logs/%s.log'." % ( fname + zi ) , 0 )
    elif en[0][0]=="MOLPRO":
      if comm( [ "grep" , "Molpro calculation terminated" , fname + zi + ".out" ] )[0]=="":
        Msg( "    WARNING - no normal termination occurred, moving on." , 0 )
        fail.append( en[1] )
        continue
      com.append( comm( [ "grep" , "SETTING HL_OUT" , fname + zi + ".out" ] ) )
      if com[-1][0]=="" or com[-1][0].split()[-1]!="AU":
        Msg( "    WARNING - no total energy found in '%s', moving on." % ( fname + zi + ".out" ) , 0 )
        fail.append( en[1] )
        continue
      else:
        EQ[ en[1] ] = str( float( com[-1][0].split()[-2] ) )
        Msg( "    Setting '%s' equal to %s." % ( en[1] , EQ[ en[1] ] ) , 0 )
        com.append( comm( [ "cp" , "%s.out" % ( fname + zi ) , "hl_logs/%s.out" % ( fname + zi ) ] ) )
        Msg( "    Output file stored as 'hl_logs/%s.log'." % ( fname + zi ) , 0 )
    elif en[0][0]=="MRCC":
      if comm( [ "grep" , "Normal termination" , fname + zi + ".log" ] )[0]=="":
        Msg( "    WARNING - no normal termination occurred, moving on." , 0 )
        fail.append( en[1] )
        continue
      com.append( comm( [ "grep" , "-B" , "5" , "Normal" , fname + zi + ".log" ] ) )
      EQ[ en[1] ] = str( float( com[-1][0].split("\n")[0].split()[-1] ) )
      Msg( "    Setting '%s' equal to %s." % ( en[1] , EQ[ en[1] ] ) , 0 )
      com.append( comm( [ "cp" , "%s.log" % ( fname + zi ) , "hl_logs/%s.log" % ( fname + zi ) ] ) )
      Msg( "    Output file stored as 'hl_logs/%s.log'." % ( fname + zi ) , 0 )
    restart.write( "%s%30s\n" % ( en[1] , EQ[ en[1] ] ) )
    restart.flush()
    Msg( "    Appending total energy to 'output/%s.restart'." % ( fname ) , 0 )
restart.close()

if len(fail)>0:
  Msg( "  ERROR - %d failed steps." % ( len(fail) ) , 1 )
else:
  Msg( "  All steps completed successfully." , 0 )
  
### solve for E ###

new = 0
K = tuple( EQ )
for l in K:
  while "(" in EQ[l]:
    pc = []
    for x,X in enumerate(EQ[l]):
      if   X=="(": O = x
      elif X==")": pc.append( x )
    for IC in pc:
      if IC>O:
        C = IC
        break
    if EQ[l][O+1] in ("*","/"): Msg( "  ERROR - could not solve for 'E', possible invalid expression in the $SOLVE input section." , 1 )
    NEW = "NEW" + str(new).zfill(3)
    EQ[NEW] = EQ[l][O+1:C]
    EQ[l] = EQ[l].replace( EQ[l][O:C+1] , NEW )
    new = new + 1

K = tuple( EQ )
S = ["/","**","*","+","-"]
nstp = 0
while True:
  nstp = nstp + 1
  try:
    float( EQ["E"] )
    break
  except:
    EQReplace( K , S )
    Z = EQCheck( K , S )
    for z,l in enumerate(K):
      if Z[z]:
        EQ[l] = EQSolve( EQ[l] , S )
  if nstp>1000:
    Msg( "  ERROR - could not solve for 'E', possible invalid expression in the $SOLVE input section." , 1 )

Msg( "  Setting final energy to %.12f." % ( float( EQ["E"] ) ) , 0 )
final = open( "output/hl_composite.en" , "w" )
final.write( "%24.12f" % ( float( EQ["E"] ) ) )
final.close()
Msg( "  Storing final energy in 'output/hl_composite.en'." , 0 )
Msg( "Normal termination" , 0 )
log.close()
ren = comm( [ "mv" , "output/hl_composite.log" , "output/" + fname + ".log" ] )

