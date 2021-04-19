#!/usr/bin/python

### Calculating glitches for a given orientation ###

import math as m
import sys

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
    if abs((_hx+K0x)**2+(_hy+K0y)**2+(_hz+K0z)**2-K*K) < (tolerance*K)**2:
        return True
    else: 
        return False

# all in A, stupid cubic cell
#?cell_a = 3.57 #?3.559
cell_a = 3.57325 # to match length of K vector
#+omega = 3.9*m.pi/float(180)
#+phi = 4.1*m.pi/float(180)
omega = 4.125*m.pi/float(180)
phi = 3.87*m.pi/float(180)

cell_a = float(sys.argv[1])
omega = float(sys.argv[2])*m.pi/float(180)
phi = float(sys.argv[3])*m.pi/float(180)

print(cell_a, omega, phi)

aR = 1/cell_a
bR = 1/cell_a
cR = 1/cell_a

Emin = 10.
Emax = 20.
numE = 5000
tolerance = 0.08

donehkls = []            # hkls for which energy is found
donehkls.append((0,0,0))
doneEn = []              # corresponding energies
doneEn.append(0)

print(CalcK0(1, omega, phi))

#kcalc1 = ( 0.99520993,  0.08922576, -0.03985014)
#kcalc = (-kcalc1[0], -kcalc1[1], -kcalc1[2])

for eN in range(0, numE+1):
        # current energy and K=1/lambda
        energ = Emax - float(Emax-Emin)*float(numE-eN)/float(numE);
        K = energ/conv_fact
        # boundary for calculation - restricted by 2K sphere - this is 2K cube
        NReci = int(2.*cell_a*K+0.1);
        # Calculate current direction of incident beam for current omega and phi
        K0x, K0y, K0z = CalcK0(K, omega, phi)
#-        print (K0x,K0y,K0z)
        # test of analit solution
#        K0x = K*kcalc[0]
#        K0y = K*kcalc[1]
#        K0z = K*kcalc[2]
#-        print (K0x,K0y,K0z)
#-        exit()
        totalInt = 0.
        # loop over all possibly excited reflections
        for hi in range(-NReci, NReci+1):
          for ki in range(-NReci, NReci+1):
            for li in range(-NReci, NReci+1):
                #skip already found reflections and 000
                if ((hi,ki,li) in donehkls): 
                    continue
                # selection rules for Di unit cell
                if not(CheckSelectionRules(hi,ki,li)):
                    continue
                # calculating vector H for each reciprocal point
                hx, hy, hz = CalcH(aR, bR, cR, hi, ki, li)
                # reject outside 2K sphere
		if m.sqrt(hx*hx+hy*hy+hz*hz) > 2.*K: 
                    continue
                # MAIN function. Checks how far each recipr. point is from current Ewald sphere. Distance should be within K*tolerance
                if IsNBeam(hx,hy,hz,K0x,K0y,K0z,K,tolerance):
#debug                    print(energ, NReci, hi, ki, li, hx, hy, hz, K0x, K0y, K0z)
                # actually here better save current energy in some array - for fitting later

                    theta = abs(0.5*m.pi-AnglVec(hx,hy,hz,K0x,K0y,K0z))
                    dSpace = cell_a/m.sqrt(hi*hi+ki*ki+li*li)
                    eTrue = conv_fact/(2*dSpace*m.sin(theta))
                    donehkls.append((hi,ki,li))
                    doneEn.append(eTrue)
                    totalInt += 1
#        print('{0:.3f}  {1:.2f}'.format(energ, totalInt))
# here should be fitting of exp. data to saved array of energies
# for fitting, for example, the distance Sum((Et-Ed)**2) might be minimized
# print best fi and alf and then fixfi,alf and run like it is now to output the spectrum

for i in range(1, len(doneEn)):
    print('{0:d},{1:d},{2:d}  {3:.5f}'.format(donehkls[i][0],donehkls[i][1],donehkls[i][2],doneEn[i]))
