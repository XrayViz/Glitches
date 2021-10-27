#!/usr/bin/python3

### Calculate error for a given file with glitches and for given UC, omega, phi ###
  
import sys
import numpy as np
import math as m

conv_fact = 12.3984

# rotate vector Ain around axis B
def Rotating(angle, AinX, AinY, AinZ, BX, BY, BZ):
    AbsB = m.sqrt(BX*BX+BY*BY+BZ*BZ)
    if AbsB==0: 
        return 0, 0, 0
    Skal = (1.-m.cos(angle))*(BX*AinX+BY*AinY+BZ*AinZ)/AbsB
    AoutX = AinX*m.cos(angle)+(BX*Skal+(BY*AinZ-BZ*AinY)*m.sin(angle))/AbsB
    AoutY = AinY*m.cos(angle)+(BY*Skal+(BZ*AinX-BX*AinZ)*m.sin(angle))/AbsB
    AoutZ = AinZ*m.cos(angle)+(BZ*Skal+(BX*AinY-BY*AinX)*m.sin(angle))/AbsB
    return AoutX, AoutY, AoutZ

def MultScal(AX, AY, AZ, BX, BY, BZ):
    return AX*BX+AY*BY+AZ*BZ

def AbsVec(AX, AY, AZ):
    return m.sqrt(AX*AX+AY*AY+AZ*AZ)

def AnglVec(AX, AY, AZ, BX, BY, BZ):
    mmm = MultScal(AX,AY,AZ,BX,BY,BZ)/(AbsVec(AX,AY,AZ)*AbsVec(BX,BY,BZ))
    if (mmm >= 1): 
        return 0.
    else: 
      return m.acos(mmm)

# not needed I think
#def MultVect(AX, AY, AZ, BX, BY, BZ):
#    CX = AY*BZ-AZ*BY
#    CY = AZ*BX-AX*BZ
#    CZ = AX*BY-AY*BX
#    return CX, CY, CZ

def CalcH(aR, bR, cR, _h, _k, _l):
#XViz    Hy = aR*_h
#XViz    Hx = bR*_k
    Hx = aR*_h
    Hy = bR*_k
    Hz = cR*_l
    return Hx, Hy, Hz

def CalcK0(K, omega, phi):
    K0x = -K
    K0y = 0.
    K0z = 0.
    # that's not the best - rotations are not commutative. 
    K0x_, K0y_, K0z_ = Rotating(omega,K0x,K0y,K0z,0,1,0)
    return Rotating(phi,K0x_,K0y_,K0z_,0,0,1)

def CheckSelectionRules(_h, _k, _l):
    return (((_h%2==0 and _k%2==0 and _l%2==0) and ((_h+_k+_l)%4==0) or (_h%2!=0 and _k%2!=0 and _l%2!=0)))

def IsNBeam(_hx,_hy,_hz,k0x,k0y,k0z,K,tolerance):
    if abs((_hx+k0x)**2+(_hy+k0y)**2+(_hz+k0z)**2-K*K) < (tolerance*K)**2:
        return True
    else: 
        return False

if (len(sys.argv) < 5):
    print("Usage: nk_calc_error_rev.py  file_with_exp_E cell_a omega phi")
    sys.exit(0)

# all in A, stupid cubic cell
cell_a = 3.57

cell_a = float(sys.argv[2])
omega = float(sys.argv[3])*m.pi/float(180)
phi = float(sys.argv[4])*m.pi/float(180)

cell_b = cell_a
cell_c = cell_a
#XViz aR = -1/cell_a
aR = 1/cell_a
bR = 1/cell_b
cR = 1/cell_c

Emin = 10.
Emax = 20.
numE = 2001
tolerance = 0.08

# reading experimental spectrum
#fileE = open("energies_exp.txt","r")  
fileE = open(sys.argv[1],"r")  
energ_list=[]
for line in fileE:
  val = float(line.strip('\n'))
  if (val > 0):
    energ_list.append(val)
fileE.close()
expLen = len(energ_list)

#for omega in np.arange(omStart, omEnd+angStep, angStep):
#  for phi in np.arange(phStart, phEnd+angStep, angStep):
donehkls = []            # hkls for which energy is found
donehkls.append((0,0,0))
doneEn = []              # corresponding energies
doneEn.append(0)
energ_dist = [0] * expLen
energ_num = [0] * expLen
# max energy and K=1/lambda
K = Emax/conv_fact
# boundary for calculation - restricted by 2K sphere - this is 2K cube
NReci = int(2.*cell_a*K+0.1);
# Calculate current direction of incident beam for current omega and phi
K0x, K0y, K0z = CalcK0(1, omega, phi)
# loop over all possibly excited reflections
for hi in range(-NReci, NReci+1):
  for ki in range(-NReci, NReci+1):
    for li in range(-NReci, NReci+1):
        #skip already found reflections
        if ((hi,ki,li) in donehkls): 
            continue
        # selection rules for Di unit cell
        if not(CheckSelectionRules(hi,ki,li)):
            continue
        # calculating vector H for each reciprocal point
        hx, hy, hz = CalcH(aR, bR, cR, hi, ki, li)
        # reject outside 2K sphere
        Hm = m.sqrt(hx*hx+hy*hy+hz*hz)
        if Hm > 2.*K: 
            continue
        # MAIN function. Checks how far each recipr. point is from current Ewald sphere. Distance should be within K*tolerance

        ang = AnglVec(hx,hy,hz,K0x,K0y,K0z)
        eTrue = 0.5*conv_fact*Hm/(m.cos(-ang))
        if eTrue >= Emin and eTrue <= Emax:
            #calculating exact energy for this reflection
            donehkls.append((hi,ki,li))
            doneEn.append(eTrue)
            
#debug                    print(hi, ki, li, energ, eTrue) 
#debug                  print(energ, NReci, hi, ki, li, hx, hy, hz, K0x, K0y, K0z)

total_dist = 0
for l in range(0, expLen):
    dist0 = 10000
    for k in range(1, len(doneEn)):
        dist = (doneEn[k]-energ_list[l])**2
        if (dist < dist0):
            dist0 = dist
    total_dist += dist0

total_dist = 1e3*m.sqrt(total_dist/float(expLen))

print('UC: {0:.6f}, Omega: {1:.5f}, Phi: {2:.5f}, Aver.Sq.Dist: {3:.5f}'.format(cell_a, omega*180./m.pi, phi*180./m.pi, total_dist))

fileE = open("best_energies","w")  
for i in range(1, len(doneEn)):
    print('{0:d},{1:d},{2:d}  {3:.6f}'.format(donehkls[i][0],donehkls[i][1],donehkls[i][2],doneEn[i]), file=fileE)
fileE.close()

