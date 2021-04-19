#include "stdio.h"
#include "string.h"
#include <math.h>
//#include "stdlib.h"

#define SQR(A) ((A)*(A))
const int MaxStrLen = 255;

int roundCpp(double xd)
{ //double _xd = xd + 103079215104.5;
  //return ((int*)&_xd)[0] >> 16;
  double t = ((xd) + 6755399441055744.0);
  return *((int *)(&t));
}

int roundMy(float x)
{ double x1 = x;
  return roundCpp(x1);
};

float medianCutoff(float arr[], size_t n, float cutoff)  // cutoff from 0 to 1
{ if (n<1) return 0;
  const float maxval = 1e30;
  int upperlim = roundMy(cutoff*((double)n+0.01));
  float _minim=arr[0];
  for (int jj=0; jj<upperlim; jj++)
  { _minim = maxval;
    int cmin = 0;
    for (int j=0; j<n; j++)
      if (arr[j]<_minim)
      { _minim = arr[j];
        cmin = j;
      }
    arr[cmin] = maxval;
  }
  return _minim;
}



bool MedianFilter1D(float* inpAr, int numX, int radX, float badVal)
{
  float* arr1 = new float[(2*radX+1)];
  float* arrOld = new float[numX];
  for (int i=0; i<numX; i++)
    arrOld[i] = inpAr[i];

  for (int xi=0; xi<numX; xi++)
    { int nuel = 0;
      for (int bxi=-radX; bxi<=radX; bxi++)
          if (xi+bxi>=0 && xi+bxi<numX)
            if (arrOld[xi+bxi]>badVal+1)
            { arr1[nuel] = arrOld[xi+bxi];
              nuel++;
            }
      if (nuel>0)
        inpAr[xi] = medianCutoff(arr1, nuel, 0.5);
      else inpAr[xi] = badVal;
    }
  delete[] arr1;
  delete[] arrOld;
  return true;
}

int main(int argc, char* argv[])
{
  printf("Processing glitch spectrum...\n");
  int radX = 1; //make a parameter ?!
  float intCut = 0.98;
  float derCut = 1.5;


  if (argc<2)
  { printf("provide file name!\n");
    return 0;
  }
//  printf("trying to open file %s\n",argv[1]);

  if (argc>2)
    radX = atoi(argv[2]);
  if (argc>3)
    derCut = atof(argv[3]);

  FILE* infil = fopen(argv[1],"rt");
  if (infil==NULL)
  { printf("File %s could not be opened!\n",argv[1]);
    return 0;
  }
  char str1[MaxStrLen];
  int numEn = 0;
  while (!feof(infil))
  { fgets(str1,MaxStrLen,infil);
    numEn++;
  }                                                
  rewind(infil);
  float* energ_list = new float[numEn];
  float* intens_list = new float[numEn];
  float* deriv_list = new float[numEn];

  numEn=0;
  while (!feof(infil))
  { fgets(str1,MaxStrLen,infil);
    if (feof(infil)) break;
    if (sscanf(str1,"%f\t%f",&energ_list[numEn],&intens_list[numEn]) > 1)
      numEn++;
  }
  fclose(infil);
  if (numEn<1) 
  { printf("Not found any energies in the file! Exiting...\n");
    return 0;
  }
  printf("Found %d energies in the list %s\n",numEn,argv[1]);

  //smoothing
  MedianFilter1D(intens_list, numEn, radX, -1);

  // normalizing
  for (int i=0; i<numEn; i++)
    deriv_list[i] = intens_list[i];
  MedianFilter1D(deriv_list, numEn, 10, -1);
  for (int i=0; i<numEn; i++)
    if (fabs(deriv_list[i])>1e-10)
      intens_list[i] /= deriv_list[i];


  //derivative
  deriv_list[0] = 0;
  for (int i=1; i<numEn; i++)
    deriv_list[i] = (intens_list[i]-intens_list[i-1])/(energ_list[i]-energ_list[i-1]);

  //output
  sprintf(str1,"%s_der",argv[1]);
  FILE* outfil = fopen(str1,"wt");
  sprintf(str1,"%s_mod",argv[1]);
  FILE* outfil1 = fopen(str1,"wt");
  for (int i=0; i<numEn; i++)
  { fprintf(outfil,"%f\t%f\n",energ_list[i],deriv_list[i]);
    fprintf(outfil1,"%f\t%f\n",energ_list[i],intens_list[i]);
  }
  fclose(outfil);
  fclose(outfil1);

  // Max intens
  float maxI = 0;
  for (int i=0; i<numEn; i++)
    if (maxI < intens_list[i]) maxI = intens_list[i];

  // Outputing only energies bigger than derCut
  sprintf(str1,"%s_glitches_plot",argv[1]);
  outfil = fopen(str1,"wt");
  sprintf(str1,"%s_glitches",argv[1]);
  outfil1 = fopen(str1,"wt");
  int numGlitches = 0;
  for (int i=3; i<numEn-1; i++)
    if (intens_list[i-1]<intCut)
    if //((deriv_list[i-1] < -derCut && deriv_list[i] > derCut) ||
       ((deriv_list[i-2] < -derCut && deriv_list[i] > derCut) ||
        (deriv_list[i-2] < -derCut && deriv_list[i+1] > derCut) ||
        (deriv_list[i-3] < -derCut && deriv_list[i+1] > derCut))
    { float cE = 0.5*(energ_list[i-1] + energ_list[i-1]);
//      if (fabs(deriv_list[i-1])<1e-5) cE = energ_list[i-1];
//      else if (fabs(deriv_list[i])<1e-5) cE = energ_list[i];
      fprintf(outfil,"%f\t%f\n",cE-0.0001,0);
      fprintf(outfil,"%f\t%f\n",cE,maxI);
      fprintf(outfil1,"%f\n",cE);
      fprintf(outfil,"%f\t%f\n",cE+0.0001,0);
      numGlitches++;
      i += 3;
    }
  fclose(outfil);
  fclose(outfil1);

  printf("Found %d glitches\n",numGlitches);

  delete[] energ_list;
  delete[] intens_list;
  delete[] deriv_list;

  return 1;
}
//---------------------------------------------------------------------------
