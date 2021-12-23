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


if len(sys.argv) < 5:
    print("The angles of a CRL at which no glitches are present in the desired E range.")
    print("Usage:")
    print("  nk_safe_range.py Emin(eV) Emax(eV) cOm(deg) cPhi(deg) a(A) [omx omy omz] [phx phy phz]")
    print("where [omx omy omz] and [phx phy phz] - rotation axes (not implemented yet)")
    sys.exit(0)

Emin = 0.001*float(sys.argv[1])  #keV
Emax = 0.001*float(sys.argv[2])
omega = float(sys.argv[3])*m.pi/float(180)   #rad
phi = float(sys.argv[4])*m.pi/float(180)
cell_a = float(sys.argv[5])

eTol = 0.01 * (Emax - Emin)

omRange = 25.5 * m.pi/float(180)
phRange = 25.5 * m.pi/float(180)

omSteps05 = 10
phSteps05 = 10

omStart = omega - 0.5*omRange
phStart = phi - 0.5*phRange

#cell_a = 3.57
aR = 1/cell_a

goodAng = []  

K = Emax/conv_fact
NReci = int(2.*cell_a*K+0.1);

for omi in range(-omSteps05, omSteps05+1):
  for phi in range(-phSteps05, phSteps05+1):
    com = omStart + omi*omRange/(2.*omSteps05)
    cph = phStart + phi*phRange/(2.*phSteps05)

    glitch = 0
    for hi in range(-NReci, NReci+1):
      for ki in range(-NReci, NReci+1):
        for li in range(-NReci, NReci+1):
    
          # selection rules for Di unit cell
          if not(CheckSelectionRules(hi,ki,li)):
              continue
          # calculating vector H for each reciprocal point
          hx, hy, hz = CalcH(aR, aR, aR, hi, ki, li)
          # reject outside 2K sphere
          Hm = m.sqrt(hx*hx+hy*hy+hz*hz)
          if Hm > 2.*K: 
              continue
    
          if hi==0 and ki==0 and li==0:
              continue
    
          K0x, K0y, K0z = CalcK0(1, com, cph)
          ang = AnglVec(hx,hy,hz,K0x,K0y,K0z)
          eTrue = 0.5*conv_fact*Hm/(m.cos(-ang))


#          print(eTrue)
          if eTrue > Emin-eTol and eTrue < Emax+eTol:
              glitch = 1
              hi = NReci+1
              ki = NReci+1
              li = NReci+1
#              print(eTrue)
              continue

#    sys.exit(0)
    if glitch == 0:
        goodAng.append((com,cph))

if len(goodAng) == 0:
    print("No  configuration is found. Try to decrease the range")
    sys.exit(0)


minSum = 1e10
bestS = 0
for i in range(0, len(goodAng)):
#    print("goodangs ", goodAng[i][0]*180/m.pi, goodAng[i][1]*180/m.pi)
    asum = m.sqrt((omega-goodAng[i][0])**2 + (phi-goodAng[i][1])**2)
    if asum < minSum:
        minSum = asum
        bestS = i

print("The good angles are, omega=", goodAng[bestS][0]*180/m.pi, "phi=", goodAng[bestS][1]*180/m.pi)

sys.exit(0)


#for omega in np.arange(omStart, omEnd+angStep, angStep):
#  for phi in np.arange(phStart, phEnd+angStep, angStep):
#    print('{0:.2f},{1:.2f} '.format(omega*180/m.pi, phi*180/m.pi))

#print('UC: {0:.6f}, Omega: {1:.5f}, Phi: {2:.5f}, Best.Sq.Dist: {3:.5f}'.format(bestUC, bestOmega*180./m.pi, bestPhi*180./m.pi, bestDist))

numang = round((omEnd-omStart+angStep)*(phEnd-phStart+angStep)/angStep/angStep)
print(numang,len(donehkls),len(doneEn))
fileE = open("plot2d","w")  
print('0 0 ', end="", file=fileE)
for i in range(0, len(donehkls)):
  print('{0:d},{1:d},{2:d} '.format(donehkls[i][0],donehkls[i][1],donehkls[i][2]), end="", file=fileE)
print(' ', file=fileE)

count=0
for omega in np.arange(omStart, omEnd+angStep, angStep):
  for phi in np.arange(phStart, phEnd+angStep, angStep):
    print('{0:.2f} {1:.2f} '.format(omega*180/m.pi, phi*180/m.pi), end="", file=fileE)
    for i in range(0, len(donehkls)):
      print('{0:.6f} '.format(doneEn[i*numang+count]), end="", file=fileE)
    count += 1
    print(' ', file=fileE)

#for i in range(0, len(donehkls)):
#  for j in range(i*numang, (i+1)*numang):
#    print('{0:.6f} '.format(doneEn[j]), end="", file=fileE)

#print(' ', file=fileE)
fileE.close()

