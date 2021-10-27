// Fitting experimental spectrum with the simulated one

// This version compares experimental energies to simulated. 
// It calculated distances from the experimental E to the closes found simulated

// (c) N.Klimova with a little help from O.Yefanov

#ifdef _OPENMP
#include <omp.h>
#endif
#include "stdio.h"
#include "string.h"
#include <math.h>
#include <ctype.h>
//#include "stdlib.h"

#define SQR(A) ((A)*(A))
const int MaxStrLen = 255;
const float conv_fact = 12.3984;
struct glitches {
  int h;              // hkls for which energy is found
  int k;
  int l;              // corresponding energies
  double E;
};

void TrimNoCh(char* str1)
{ int cstlen = strlen(str1);
  for (int i=cstlen-1; i>0; i--)
	if (!isalnum(str1[i]) && str1[i]!='-' && str1[i] != '/' && str1[i] != '\\' && str1[i] != '.' && str1[i] != '+')
	  str1[i] = 0;
	else break;
  while ((!isalnum(str1[0]) && str1[0] != '/' && str1[0] != '\\' && str1[0] != '.' && str1[0] != '-' && str1[0] != '+') && strlen(str1)>0)
	for (int i=0; i<strlen(str1); i++)
	  str1[i] = str1[i+1];
}

void strLower1(char *inp)
{ for (int i=0; i<strlen(inp); i++)
	inp[i] = tolower(inp[i]);
}

double AbsVec(double AX, double AY, double AZ)
{
  return sqrt(AX*AX+AY*AY+AZ*AZ);
}

double MultScal(double AX, double AY, double AZ, double BX, double BY, double BZ)
{
  return AX*BX+AY*BY+AZ*BZ;
}

int Rotating(double angle, double AinX,   double AinY,   double AinZ,
                           double BX,     double BY,     double BZ,
                           double* AoutX, double* AoutY, double* AoutZ)
{
//  double AbsB = AbsVec(BX,BY,BZ);
  double AbsB = sqrt(BX*BX+BY*BY+BZ*BZ);
  if (AbsB==0) return -1;
  double Skal = (1.-cos(angle))*(BX*AinX+BY*AinY+BZ*AinZ)/AbsB;
  *AoutX = AinX*cos(angle)+(BX*Skal+(BY*AinZ-BZ*AinY)*sin(angle))/AbsB;
  *AoutY = AinY*cos(angle)+(BY*Skal+(BZ*AinX-BX*AinZ)*sin(angle))/AbsB;
  *AoutZ = AinZ*cos(angle)+(BZ*Skal+(BX*AinY-BY*AinX)*sin(angle))/AbsB;
  return 1.;
}

double AnglVec(double AX, double AY, double AZ, double BX, double BY, double BZ)
{
  double mmm = MultScal(AX,AY,AZ,BX,BY,BZ)/(AbsVec(AX,AY,AZ)*AbsVec(BX,BY,BZ));
  if (mmm >= 1) return 0.;
  else return acos(mmm);
}

void CalcH(int _h, int _k, int _l, double aR, double bR, double cR,
           double* _hx, double* _hy, double* _hz)
{
  *_hx = aR*_h;
  *_hy = bR*_k;
  *_hz = cR*_l;
}

int CalcK0(double K, double omega, double phi, double* k0x, double* k0y, double* k0z)
{
  double _k0x = -K;
  double _k0y = 0.;
  double _k0z = 0.;
  Rotating(omega,_k0x,_k0y,_k0z,0.,1.,0.,&_k0x,&_k0y,&_k0z);
  Rotating(phi,_k0x,_k0y,_k0z,0.,0.,1.,&_k0x,&_k0y,&_k0z);
  *k0x = _k0x;
  *k0y = _k0y;
  *k0z = _k0z;
  return 1;
}

bool CheckSelectionRules(int _h, int _k, int _l)
{
  return (((_h%2==0 && _k%2==0 && _l%2==0) && ((_h+_k+_l)%4==0) || (_h%2!=0 && _k%2!=0 && _l%2!=0)));
}

bool IsNBeam(double _hx, double _hy, double _hz, double k0x, double k0y, double k0z,
             double K, double tolerance)
{
  if (fabs(SQR(_hx+k0x)+SQR(_hy+k0y)+SQR(_hz+k0z)-K*K) < SQR(tolerance*K))
    return true;
  else 
    return false;
}


int CalcGlitches(double omega, double phi, int numE, int numFi, double tolerance,
                 double Emin, double Emax, double cell_a, double aR, double bR, double cR,
                 glitches *donehkls)
{
        int doneN = 1; // 0,0,0 - always skip
        for (int iE=0; iE<numE; iE++)
        { donehkls[iE].h = 0;
          donehkls[iE].k = 0;
          donehkls[iE].l = 0;
          donehkls[iE].E = 0.;
        }

        for (int ie=0; ie<=numE; ie++)
        {
          double energ = Emax - (Emax-Emin)*(double)(numE-ie)/(double)numE;
          double K = energ/conv_fact;
          // boundary for calculation - restricted by 2K sphere - this is 2K cube
          int NReci = round(2.*cell_a*K+0.1);
          // Calculate current direction of incident beam for current omega and phi
          double K0x,K0y,K0z;
          CalcK0(K, omega, phi, &K0x, &K0y, &K0z);
          // loop over all possibly excited reflections
          for (int hi=-NReci; hi<=NReci; hi++)
            for (int ki=-NReci; ki<=NReci; ki++)
              for (int li=-NReci; li<=NReci; li++)
              { bool contC = false;
                for (int i=0; i<doneN; i++)
                  if (donehkls[i].h==hi && donehkls[i].k==ki && donehkls[i].l==li) 
                  { contC = true;
                    break;
                  }
                if (contC) continue;
  
                // selection rules for Di unit cell
                if (!CheckSelectionRules(hi,ki,li))
                  continue;
  
                // calculating vector H for each reciprocal point
                double hx, hy, hz;
                CalcH(hi, ki, li, aR, bR, cR, &hx, &hy, &hz);
  
                // reject outside 2K sphere
                if (sqrt(hx*hx+hy*hy+hz*hz) > 2.*K) 
                  continue;
  
                // MAIN function. Checks how far each recipr. point is from current Ewald sphere. Distance should be within K*tolerance
                if (IsNBeam(hx,hy,hz,K0x,K0y,K0z,K,tolerance))
                {   // calculating exact energy for this reflection
                  double theta = fabs(0.5*M_PI - AnglVec(hx,hy,hz,K0x,K0y,K0z));
                  float dSpace = cell_a/sqrt(hi*hi+ki*ki+li*li);
                  float eTrue = conv_fact/(2*dSpace*sin(theta));
    
                  donehkls[doneN].h = hi;
                  donehkls[doneN].k = ki;
                  donehkls[doneN].l = li;
                  donehkls[doneN].E = eTrue;
                  doneN++;
//                  printf("Glitch with hkl %d, %d, %d at E=%0.3f\n",hi,ki,li,eTrue);
                }
              }
        }

        return doneN;
}

int CalcGlitchesFast(double omega, double phi, int numE, int numFi, double tolerance,
                 double Emin, double Emax, double cell_a, double aR, double bR, double cR,
                 glitches *donehkls)
{
        int doneN = 1; // 0,0,0 - always skip
        for (int iE=0; iE<numE; iE++)
        { donehkls[iE].h = 0;
          donehkls[iE].k = 0;
          donehkls[iE].l = 0;
          donehkls[iE].E = 0.;
        }


          double Kmax = Emax/conv_fact;
          // boundary for calculation - restricted by 2K sphere - this is 2K cube
          int NReci = round(2.*cell_a*Kmax+0.1);
          // Calculate current direction of incident beam for current omega and phi
          double K0x,K0y,K0z;   //ort
          CalcK0(1, omega, phi, &K0x, &K0y, &K0z);
          // loop over all possibly excited reflections
          for (int hi=-NReci; hi<=NReci; hi++)
            for (int ki=-NReci; ki<=NReci; ki++)
              for (int li=-NReci; li<=NReci; li++)
              { bool contC = false;
                for (int i=0; i<doneN; i++)
                  if (donehkls[i].h==hi && donehkls[i].k==ki && donehkls[i].l==li) 
                  { contC = true;
                    break;
                  }
                if (contC) continue;
  
                // selection rules for Di unit cell
                if (!CheckSelectionRules(hi,ki,li))
                  continue;
  
                // calculating vector H for each reciprocal point
                double hx, hy, hz;
                CalcH(hi, ki, li, aR, bR, cR, &hx, &hy, &hz);
  
                // reject outside 2K sphere
                float Hm = sqrt(hx*hx+hy*hy+hz*hz);
                if (Hm > 2.*Kmax) 
                  continue;
  
                // MAIN calculation
//old
                float ang = AnglVec(hx,hy,hz,K0x,K0y,K0z);
//old
                float eTrue1 = 0.5*conv_fact*Hm/(cos(-ang));

                float eTrue = -0.5*conv_fact*Hm*Hm/MultScal(hx,hy,hz,K0x,K0y,K0z) ;

                if (eTrue >= Emin && eTrue <= Emax) 
                { 
                  donehkls[doneN].h = hi;
                  donehkls[doneN].k = ki;
                  donehkls[doneN].l = li;
                  donehkls[doneN].E = eTrue;
                  doneN++;
//                  printf("k0: %0.3f %0.3f %0.3f\n", K0x,K0y,K0z);
//                  printf("Glitch with hkl %d, %d, %d at E=%0.3f alt E=%0.3f\n",hi,ki,li,eTrue,eTrue1);
                }
              }

        return doneN;
}


int main(int argc, char* argv[])
{
  printf("Fitting...\n");

  if (argc<2)
  { printf("ERROR! A text file with enegries is needed!\n");
    printf("EXAMPLE: ny_fit energies.txt\n");
    return 0;
  }

  double tolerance = 0.08;
  double Emin = 0.;
  double Emax = 0.;
  double EStep = 0.005;
  double ucStep = 0.0001;
  double ucStart = 3.5725;
  double ucEnd = 3.5740;
  double angStep = 0.005*M_PI/180.;
  double omStart = 4.10*M_PI/180.;
  double omEnd = 4.14*M_PI/180.;
  double phStart = 3.86*M_PI/180.;
  double phEnd = 3.90*M_PI/180.;

  char fName[MaxStrLen];

  // Reading configuration from file
  char confName[MaxStrLen];
  strcpy(confName,argv[1]);
  FILE* fileC = fopen(confName,"rt");  // command line!
  if (fileC==NULL)
  { printf("Config file %s not found!\n",confName);
    return 0;
  }
                               
  char separat = '=';
  char parameter[MaxStrLen];
  char value[MaxStrLen];

  while (!feof(fileC))
  { fgets(parameter,MaxStrLen,fileC);
    if (strchr(parameter,'=')==NULL) continue;
    if (parameter[0]=='#') continue;
    strcpy(value,strchr(parameter,separat)+1);
    *strchr(parameter,separat) = 0;
    TrimNoCh(parameter);
    TrimNoCh(value);
    strLower1(parameter);
    if (strcmp(parameter,"filename")==0)
      sscanf(value, "%s", &fName);
    else if (strcmp(parameter,"tolerance")==0)
      sscanf(value, "%lf", &tolerance);
    else if (strcmp(parameter,"emin")==0)
      sscanf(value, "%lf", &Emin);
    else if (strcmp(parameter,"emax")==0)
      sscanf(value, "%lf", &Emax);
    else if (strcmp(parameter,"estep")==0)
      sscanf(value, "%lf", &EStep);
    else if (strcmp(parameter,"ucmin")==0)
      sscanf(value, "%lf", &ucStart);
    else if (strcmp(parameter,"ucmax")==0)
      sscanf(value, "%lf", &ucEnd);
    else if (strcmp(parameter,"ucstep")==0)
      sscanf(value, "%lf", &ucStep);
    else if (strcmp(parameter,"ommin")==0)
      sscanf(value, "%lf", &omStart);
    else if (strcmp(parameter,"ommax")==0)
      sscanf(value, "%lf", &omEnd);
    else if (strcmp(parameter,"phmin")==0)
      sscanf(value, "%lf", &phStart);
    else if (strcmp(parameter,"phmax")==0)
      sscanf(value, "%lf", &phEnd);
    else if (strcmp(parameter,"angstep")==0)
      sscanf(value, "%lf", &angStep);
  }
  fclose(fileC);

  // reading experimental spectrum
  FILE* fileE = fopen(fName,"rt");
  if (fileE==NULL)
  { printf("List of energies %s not found!\n",fName);
    return 0;
  }

  int numFi = 0;
  char str1[MaxStrLen];
  while (!feof(fileE))
  { fgets(str1,MaxStrLen,fileE);
    numFi++;
  }
  rewind(fileE);
  double* energ_list = new double[numFi];

  numFi=0;
  while (!feof(fileE))
  { fgets(str1,MaxStrLen,fileE);
    if (feof(fileE)) break;
    sscanf(str1,"%lf",&energ_list[numFi]);
    numFi++;
  }
  fclose(fileE);
  if (numFi<1) 
  { printf("Not found any energies in the file! Exiting...\n");
    return 0;
  }
  printf("Found %d energies in the list %s\n",numFi,fName);

  // If the energy range is not set - use from the file
  if (Emax<0.1)
  { Emax = 0;
    Emin = 1e10;
    for (int i=0; i<numFi; i++)
    { if (energ_list[i] < Emin) Emin = energ_list[i];
      if (energ_list[i] > Emax) Emax = energ_list[i];
    }
    Emin = (int)Emin;
    Emax = (int)Emax + 1;
  }


  printf("Fitting parameters:\n");
  printf("Filename with experimental energies: %s\n",fName);
  printf("Energy tolerance=%0.3f, Emin=%0.2f, Emax=%0.2f, Estep=%0.4f\n",tolerance, Emin, Emax, EStep);
  printf("UC min=%0.4f, max=%0.4f, step=%0.4f\n",ucStart, ucEnd, ucStep);
  printf("Omega min=%0.3f, max=%0.3f, Phi min=%0.3f, max=%0.3f, Step=%0.3f\n",omStart, omEnd, phStart, phEnd, angStep);
  omStart *= M_PI/180.;
  omEnd *= M_PI/180.;
  phStart *= M_PI/180.;
  phEnd *= M_PI/180.;
  angStep *= M_PI/180.;

  int numUC = round((ucEnd-ucStart)/ucStep) + 1;
  int numOm = round((omEnd-omStart)/angStep) + 1;
  int numPh = round((phEnd-phStart)/angStep) + 1;
  int numE =  round((Emax-Emin)/EStep) + 1;
//  numE = 1000;

  printf("%d %d %d %d\n", numUC, numOm, numPh, numE);

#ifdef _OPENMP
  printf("OMP is on. Found %d threads\n",omp_get_max_threads());
#endif

  double *cUCs = new double[numUC];
  double *omegas = new double[numUC];
  double *phis = new double[numUC];
  double *gDiss = new double[numUC];

  // loop over possible UCs
#pragma omp parallel for shared(cUCs, omegas, phis, gDiss)
  for (int iuc=0; iuc<numUC; iuc++)  // loop for UC
  {
    double cUC = ucStart;
    if (numUC > 1) cUC += (ucEnd-ucStart)*(double)iuc/(double)(numUC-1);
    double cell_a = cUC; 
    double cell_b = cell_a;
    double cell_c = cell_a;
    //#XViz aR = -1/cell_a
    double aR = 1/cell_a;
    double bR = 1/cell_b;
    double cR = 1/cell_c;

    glitches *donehkls = new glitches[numE];

    double goodDist = 100000000; //?????
    int goodNum = 0; 
    double goodPhi = 0.;
    double goodOmega = 0.;

    // loops over angles
    for (int iph=0; iph<numPh; iph++)
      for (int iom=0; iom<numOm; iom++)
      {
        double omega = omStart;
        if (numOm > 1) omega += (omEnd-omStart)*(double)iom/(double)(numOm-1);
        double phi = phStart;
        if (numPh > 1) phi += (phEnd-phStart)*(double)iph/(double)(numPh-1);
  
        int doneN = CalcGlitchesFast(omega, phi, numE, numFi, tolerance,
                     Emin, Emax, cell_a, aR, bR, cR, donehkls);

        // here calculating the "distance" in energies between exp and fount at these omega and phi
        double total_dist = 0;
        int num_dist = 0;
        for (int ifi=0; ifi<numFi; ifi++)
        { if (energ_list[ifi] <= Emin || energ_list[ifi] >= Emax) continue;
          double dist0 = 10000;
          for (int iD=0; iD<doneN; iD++)
          {
            double dist = SQR(donehkls[iD].E - energ_list[ifi]);
            if (dist < dist0)
              dist0 = dist;
          }
          total_dist += dist0;
          num_dist++;
        }
        total_dist = sqrt(total_dist/num_dist)*1e3;

        if (total_dist < goodDist)
        { goodDist = total_dist;
          goodNum = numFi;
          goodOmega = omega*180./M_PI;
          goodPhi = phi*180./M_PI;
        }
  
      }

    printf("UC: %0.5f, Omega: %0.3f, Phi: %0.3f, Aver.Sq.Dist: %0.5feV, N/D: %0.3f\n", cUC, goodOmega, goodPhi, goodDist, goodNum/goodDist);
//    printf("UC: %0.5f, Omega: %0.3f, Phi: %0.3f, Dist: %0.3e, Num: %d, N/D: %0.0f\n", cUC, goodOmega, goodPhi, sqrt(goodDist), goodNum, goodNum/goodDist);

    cUCs[iuc] = cUC;
    omegas[iuc] = goodOmega;
    phis[iuc] = goodPhi;
//    gDiss[iuc] = goodNum/goodDist;
    gDiss[iuc] = goodDist;

    delete[] donehkls;

  }

  double bestDist = 1e10; 
  double bestPhi = 0.;
  double bestOmega = 0.;
  double bestUC = 0.;
  for (int iuc=0; iuc<numUC; iuc++)  // loop for UC
    if (gDiss[iuc] < bestDist)
    { bestDist = gDiss[iuc];
      bestOmega = omegas[iuc];
      bestPhi = phis[iuc];
      bestUC = cUCs[iuc];
    }
  delete[] cUCs;
  delete[] omegas;
  delete[] phis;
  delete[] gDiss;

  printf("BEST: UC: %0.5f, Omega: %0.3f, Phi: %0.3f, Aver.Sq.Dist: %0.5feV\n", bestUC, bestOmega, bestPhi, bestDist);

  // output the best enegries
  glitches *donehkls = new glitches[numE];
  int doneN = CalcGlitches(bestOmega*M_PI/180., bestPhi*M_PI/180., numE, numFi, tolerance,
                     Emin, Emax, bestUC, 1/bestUC, 1/bestUC, 1/bestUC, donehkls);
  FILE* outfile;
  // Saving calculated indexed glitch spectrum
//  sprintf(str1,"%s_energies",fName);
//  outfile = fopen(str1,"wt");
//  for (int iE=1; iE<doneN; iE++)
//    fprintf(outfile,"%d,%d,%d\t%0.5f\n",donehkls[iE].h,donehkls[iE].k,donehkls[iE].l,donehkls[iE].E);
//  fclose(outfile);

  // saving indexed experimental glitch spectrum
  sprintf(str1,"%s_indexed",fName);
  outfile = fopen(str1,"wt");
  for (int ifi=0; ifi<numFi; ifi++)
  { if (energ_list[ifi] <= Emin || energ_list[ifi] >= Emax) continue;
    double dist0 = 10000;
    int num0 = 0;
    for (int iD=0; iD<doneN; iD++)
    { double dist = SQR(donehkls[iD].E - energ_list[ifi]);
      if (dist < dist0)
      { dist0 = dist;
        num0 = iD;
      }
    }
    fprintf(outfile,"%d,%d,%d\t%0.5f\n",donehkls[num0].h,donehkls[num0].k,donehkls[num0].l,energ_list[ifi]);
  }
  fclose(outfile);


  delete[] donehkls;

  delete[] energ_list;
  return 1;
}
//---------------------------------------------------------------------------
