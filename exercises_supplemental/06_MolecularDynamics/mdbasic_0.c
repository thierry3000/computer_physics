
/* [[pr_02_1 - all pairs, two dimensions]] */


/*********************************************************************

  (C) 2004  D. C. Rapaport

  This software is copyright material accompanying the book
  "The Art of Molecular Dynamics Simulation", 2nd edition,
  by D. C. Rapaport, published by Cambridge University Press (2004).

**********************************************************************/


#define NDIM  2

#include "include/in_mddefs.h"

#include "include/in_rand.c"
#include "include/in_errexit.c"


// for c++ plotting
#include <cairo/cairo.h>

// definition (macros) are found in the file in_vdefs.h

typedef struct {
  VecR r, rv, ra;
} Mol;


Mol *mol;
VecR region, vSum;
VecI initUcell;
real deltaT, density, rCut, temperature, timeNow, uSum, velMag, vvSum;
Prop kinEnergy, totEnergy;
int moreCycles, nMol, stepAvg, stepCount, stepEquil, stepLimit, stepSnap;
real virSum;
Prop pressure;
FILE * datafile_main;


// input paramter from mdbasic.in (NameR -> real; NameI -> integer)
NameList nameList[] = {
  NameR (deltaT),	// Zeitschritt
  NameR (density),	// Dichte der Teilchen
  NameI (initUcell),	// Anzahl der Teilchen in x und y Richtung
  NameI (stepAvg),	// Anzahl der Schritte, nachdem die Observablen ausgegeben werden
  NameI (stepEquil),
  NameI (stepLimit),	// Anzahl der Schritte insgesamt -> Programmabbruch
  NameI (stepSnap),	// Anzahl der Schritte, nachdem ein Bild erzeugt wird
  NameR (temperature),
};

#include "include/in_namelist.c"

int main (int argc, char **argv)
{
  GetNameList (argc, argv);
  datafile_main = fopen("output.dat","w");
  PrintNameList (stdout);		//gibt Parameter an die Standardausgabe aus 
  PrintNameList (datafile_main);	//gibt Parameter in die Datei "output.dat" aus
  SetParams ();
  SetupJob ();
  OutputPng ();
  moreCycles = 1;
  while (moreCycles) {
    SingleStep ();
    if (stepCount % stepSnap == 0) OutputPng ();
    if (stepCount >= stepLimit) moreCycles = 0;
  }
}


void SingleStep ()
{
  ++ stepCount;
  timeNow = stepCount * deltaT;
  LeapfrogStep (1);
  ApplyBoundaryCond ();
  ComputeForces ();
  LeapfrogStep (2);
  EvalProps ();
  AccumProps (1);
  if (stepCount % stepAvg == 0) {
    AccumProps (2);
    PrintSummary (stdout);
    PrintSummary (datafile_main);
    AccumProps (0);
  }
}

void SetupJob ()
{
  AllocArrays ();
  stepCount = 0;
  InitCoords ();
  InitVels ();
  InitAccels ();
  AccumProps (0);
}

void SetParams ()
{
  rCut = pow (2., 1./6.);  // = 1.122462048
  VSCopy (region, 1. / sqrt (density), initUcell);
  nMol = VProd (initUcell);
  velMag = sqrt (NDIM * (1. - 1. / nMol) * temperature);
}

void AllocArrays ()
{
  AllocMem (mol, nMol, Mol);
}

void ComputeForces ()
{
  VecR dr;
  real fcVal, rr, rrCut, rri, rri3;
  int j1, j2, n;

  rrCut = Sqr (rCut);
  DO_MOL VZero (mol[n].ra);
  uSum = 0.;
  virSum = 0.;
  for (j1 = 0; j1 < nMol - 1; j1 ++) {
    for (j2 = j1 + 1; j2 < nMol; j2 ++) {
      VSub (dr, mol[j1].r, mol[j2].r);
      VWrapAll (dr);
      rr = VLenSq (dr);
      if (rr < rrCut) {
        rri = 1. / rr;
        rri3 = Cube (rri);
        fcVal = 48. * rri3 * (rri3 - 0.5) * rri;
        VVSAdd (mol[j1].ra, fcVal, dr);
        VVSAdd (mol[j2].ra, - fcVal, dr);
        uSum += 4. * rri3 * (rri3 - 1.) + 1.;
        virSum += fcVal * rr;
      }
    }
  }
}


void LeapfrogStep (int part)
{
  int n;

  if (part == 1) {
    DO_MOL {
      VVSAdd (mol[n].rv, 0.5 * deltaT, mol[n].ra);
      VVSAdd (mol[n].r, deltaT, mol[n].rv);
    }
  } else {
    DO_MOL VVSAdd (mol[n].rv, 0.5 * deltaT, mol[n].ra);
  }
}


void ApplyBoundaryCond ()
{
  int n;

  DO_MOL VWrapAll (mol[n].r);
}


void InitCoords ()
{
  VecR c, gap;
  int n, nx, ny;

  VDiv (gap, region, initUcell);
  n = 0;
  for (ny = 0; ny < initUcell.y; ny ++) {
    for (nx = 0; nx < initUcell.x; nx ++) {
      VSet (c, nx + 0.5, ny + 0.5);
      VMul (c, c, gap);
      VVSAdd (c, -0.5, region);
      mol[n].r = c;
      ++ n;
    }
  }
}


void InitVels ()
{
  int n;

  VZero (vSum);
  DO_MOL {
    VRand (&mol[n].rv);
    VScale (mol[n].rv, velMag);
    VVAdd (vSum, mol[n].rv);
  }
  DO_MOL VVSAdd (mol[n].rv, - 1. / nMol, vSum);
}


void InitAccels ()
{
  int n;

  DO_MOL VZero (mol[n].ra);
}


void EvalProps ()
{
  real vv, vvMax;
  int n;

  VZero (vSum);
  vvSum = 0.;
  DO_MOL {
    VVAdd (vSum, mol[n].rv);
    vv = VLenSq (mol[n].rv);
    vvSum += vv;
  }
  kinEnergy.val = 0.5 * vvSum / nMol;
  totEnergy.val = kinEnergy.val + uSum / nMol;
  pressure.val = density * (vvSum + virSum) / (nMol * NDIM);
}


void AccumProps (int icode)
{
  if (icode == 0) {
    PropZero (totEnergy);
    PropZero (kinEnergy);
    PropZero (pressure);
  } else if (icode == 1) {
    PropAccum (totEnergy);
    PropAccum (kinEnergy);
    PropAccum (pressure);
  } else if (icode == 2) {
    PropAvg (totEnergy, stepAvg);
    PropAvg (kinEnergy, stepAvg);
    PropAvg (pressure, stepAvg);
  }
}


void PrintSummary (FILE *fp)
{
  fprintf (fp,
     "%5d %8.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f\n",
     stepCount, timeNow, VCSum (vSum) / nMol, PropEst (totEnergy),
     PropEst (kinEnergy), PropEst (pressure));
  fflush (fp);
}


void OutputPng () 
{
int k;
int image_factor=10;
char number[10],filename[40]="moviedata/Step";

sprintf(number,"%06d",stepCount);
strcat(filename,number);
strcat(filename, ".png");

for ( k = 0; k < nMol ; k++) {
  }
	cairo_surface_t *surface = cairo_image_surface_create (CAIRO_FORMAT_ARGB32, image_factor*region.x, image_factor*region.y);
	cairo_t *cr = cairo_create (surface);

	cairo_set_source_rgb (cr, 1.0, 1.0, 1.0);
	cairo_set_line_width(cr, 1);
	cairo_rectangle(cr, 0, 0, image_factor*region.x, image_factor*region.y);
	cairo_fill(cr);
	cairo_set_source_rgb (cr, 0.0, 0.0, 0.0);

for ( k = 0; k < nMol ; k++) {
	cairo_arc(cr, image_factor*(0.5*region.x + mol[k].r.x), image_factor*(0.5*region.y + mol[k].r.y), image_factor*0.5, 0, 2*M_PI);
	cairo_fill(cr);
  }
        cairo_destroy (cr);
        cairo_surface_write_to_png (surface, filename);
        cairo_surface_destroy (surface);
}
