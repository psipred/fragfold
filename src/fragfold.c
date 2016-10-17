/***************************************************
 *                    FRAGFOLD                     *
 *         By David T. Jones    April 2004         *
 *             UCL Bioinformatics Group            *
 ***************************************************/

/* This software is Copyright (C) 1995, 2005 David T. Jones */
/* This software may not be distributed, copied or used without the permission of the Copyright holders */

#define FragfoldVersion	"4.80"
#define Last_Edit_Date	"5th Sept 2014"

/*
 * This program attempts to fold a sequence into a plausible tertiary
 * structure by joining together protein fragments.
 */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <stdarg.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include <fcntl.h>
#include <time.h>
#ifdef WIN32
#include <windows.h>
#include <io.h>
#include <process.h>
#else
#include <unistd.h>
#include <sys/time.h>
#include <sys/times.h>
#endif

#define NPAIRS 2

/* Define SAVETRAJ to save intermediate models */
#define noSAVETRAJ

/* Define RASMOL to monitor progress with RASMOL */
#define noRASMOL

/* SUPERSEC defined to use supersecondary fragments */
#define SUPERSEC

#define ELITIST

#define MAXSEQLEN 2000

#define SRDCUTOFF 15.0F
#define LRDCUTOFF 15.0F

#define SRTOPOMAX 8

/* Definitions for potentials */

#define TOPOMAX 11
#define INTERVALS 20
#define NACC 5
#define OOIMIN 6
#define OOIMAX 30
#define OOIDIV 1
#define NRR 11
#define RRWIDTH (1.0F)
#define RTCONST (0.582F)

#define OOICUTOFF (10.0F)


/* Constants for calculating CB coords (Prot. Eng. Vol.2 p.121) */
#define	TETH_ANG 0.9128F
#define CACBDIST 1.538F


/* Utility definitions */

#define FALSE 0
#define TRUE 1
#define BIG (1000000)
#define VBIG (1e32F)
#define PI (3.1415927F)

#define SQR(x) ((x)*(x))
#define MIN(x,y) (((x)<(y))?(x):(y))
#define MAX(x,y) (((x)>(y))?(x):(y))
#define CH alloc_verify(), printf("Heap OK at line : %d.\n",__LINE__);
#define veczero(v) memset(v, 0, sizeof(v))
#define vecadd(a,b,c) (a[0]=b[0]+c[0],a[1]=b[1]+c[1],a[2]=b[2]+c[2])
#define vecsub(a,b,c) (a[0]=b[0]-c[0],a[1]=b[1]-c[1],a[2]=b[2]-c[2])
#define vecscale(a,b,c) ((a[0]=b[0]*(c)),(a[1]=b[1]*(c)),(a[2]=b[2]*(c)))
#define vecprod(a,b,c) ((a[0]=b[1]*c[2]-b[2]*c[1]),(a[1]=b[2]*c[0]-b[0]*c[2]),(a[2]=b[0]*c[1]-b[1]*c[0]))
#define dotprod(a,b) (a[0]*b[0]+a[1]*b[1]+a[2]*b[2])
#define veccopy(a,b) ((a[0]=b[0]),(a[1]=b[1]),(a[2]=b[2]))
#define distsq(a,b) (SQR(a[0]-b[0])+SQR(a[1]-b[1])+SQR(a[2]-b[2]))
#define dist(a,b) (sqrtf(distsq(a,b)))

char alnfname[160], confname[160], hbfname[160], datadir[160];

int vtrange=BIG, rrcbflag;

float     srwt = 1.0, lrwt = 1.0, solvwt = 1.0, hbwt = 1.0, compactwt = 1.0, stericwt = 1.0, dswt = 0.0, srrwt = 0.0, lrrwt = 0.0, targwt = 0.0;

/* Simulated annealing parameter (total steps) */
int MAXSTEPS = 1000000;

/* Simulated annealing parameter (initial "temperature") */
float INITEMP = 0.75;

/* Simulated annealing parameter (initial "temperature") */
float TRATIO = 0.75;

/* Parameter for thermodynamic annealing */
float KAVALUE = 100000.0;

/* Max. number of supersecondary fragments per residue */
int MAXFRAGS = 5;

/* Max. number of fixed length fragments per residue */
int MAXFRAGS2  = 25;

/* Min. length of fixed length fragments per residue */
int MINFRAGS2LEN  = 9;

/* Max. length of fixed length fragments per residue */
int MAXFRAGS2LEN  = 16;

/* Refinement search start */
int REFSTART = 0;

/* Refinement search end */
int REFEND = MAXSEQLEN;

/* Refinement search range (max RMSD away from target) */
float REFRANGE = 3.0;

/* Max CPU time */
int MAXTIME = 999999999;

char     buf[MAXSEQLEN], **seq;
char     tplt_ss[MAXSEQLEN];
short    tcbooi[MAXSEQLEN], relacc[MAXSEQLEN], gaps[MAXSEQLEN], svrflg[MAXSEQLEN], segidx[MAXSEQLEN], chnbrk[MAXSEQLEN];
int      firstid, seqlen, nsvr, nseqs;
float    cn_energy, econtrib[MAXSEQLEN], e_min = VBIG;
unsigned short chcount[MAXSEQLEN];

float    param[10], exprad, expsurf, expmean, ecmin, ecmax;

float    ftotwt;

/* Verbosity flag */
int verboseflg = TRUE;

/* Known disulphide list */
int n_ds = 0;
int ds_from[10], ds_to[10];

char *outpdbn = "fold.pdb";
char *mqapcmd = NULL;

const char     *rnames[] =
{
    "ALA", "ARG", "ASN", "ASP", "CYS",
    "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO",
    "SER", "THR", "TRP", "TYR", "VAL",
    "UNK", "UNK"
};

enum aacodes
{
    ALA, ARG, ASN, ASP, CYS,
    GLN, GLU, GLY, HIS, ILE,
    LEU, LYS, MET, PHE, PRO,
    SER, THR, TRP, TYR, VAL,
    GAP, UNK
};

const char     *atmnames[] =
{
    "CA", "CB", "O ", "N ", "C "
};

enum atmcodes
{
    CAATOM, CBATOM, OATOM, NATOM, CATOM
};

enum sstypes
{
  COIL = 1, HELIX = 2, STRAND = 4
};

const char     *sscodes = ".CH.E";

enum PAIRS
{
    CA_CA, CB_CB, CB_N, N_CB, CB_O, O_CB, O_N, N_O
};

enum MODES
{
    FOLDM, REFINEM, MODELM, EVALM
};

const int       pairs[8][2] =
{
    {CAATOM, CAATOM},
    {CBATOM, CBATOM},
    {CBATOM, NATOM},
    {NATOM, CBATOM},
    {CBATOM, OATOM},
    {OATOM, CBATOM},
    {OATOM, NATOM},
    {NATOM, OATOM}
};

const float mindist[TOPOMAX][4][4] =
{
  {
    { 2.71, 2.28, 3.14, 1.91 },
    { 2.56, 2.73, 3.15, 2.05 },
    { 1.86, 1.73, 1.65, 1.86 },
    { 2.55, 1.39, 3.16, 2.29 },
  },
  {
    { 4.38, 3.62, 3.38, 3.32 },
    { 4.12, 3.68, 2.66, 3.26 },
    { 3.08, 2.61, 2.30, 2.37 },
    { 4.18, 3.55, 2.71, 3.47 },
  },
  {
    { 3.89, 2.90, 3.04, 3.84 },
    { 2.61, 2.32, 2.26, 3.03 },
    { 2.79, 2.06, 2.53, 2.63 },
    { 3.58, 2.69, 2.73, 3.51 },
  },
  {
    { 3.80, 2.77, 2.68, 3.61 },
    { 2.98, 2.49, 1.97, 3.53 },
    { 2.67, 1.73, 2.67, 2.49 },
    { 3.45, 2.33, 2.63, 3.46 },
  },
  {
    { 3.50, 2.51, 3.06, 3.22 },
    { 2.54, 2.50, 1.69, 2.79 },
    { 2.80, 2.03, 2.76, 2.53 },
    { 3.16, 2.94, 2.61, 3.37 },
  },
  {
    { 3.75, 2.78, 3.07, 3.66 },
    { 2.45, 2.22, 2.02, 2.61 },
    { 2.98, 1.83, 3.07, 2.66 },
    { 3.33, 2.90, 2.79, 3.72 },
  },
  {
    { 3.78, 2.45, 3.16, 3.51 },
    { 2.81, 2.83, 2.57, 2.89 },
    { 2.89, 1.83, 3.08, 2.72 },
    { 3.57, 3.25, 2.69, 3.80 },
  },
  {
    { 3.96, 3.12, 3.00, 3.71 },
    { 3.09, 2.51, 2.90, 2.87 },
    { 3.01, 1.86, 3.09, 2.66 },
    { 3.87, 3.41, 2.70, 3.75 },
  },
  {
    { 3.71, 3.30, 3.01, 3.80 },
    { 3.49, 2.64, 1.88, 3.36 },
    { 3.09, 2.62, 3.01, 2.69 },
    { 3.37, 2.52, 2.68, 3.56 },
  },
  {
    { 3.91, 2.63, 3.10, 3.81 },
    { 2.72, 2.62, 2.13, 2.84 },
    { 2.98, 1.93, 3.17, 2.65 },
    { 3.53, 2.34, 2.74, 3.73 },
  },
  {
    { 3.96, 2.86, 3.05, 3.69 },
    { 2.62, 2.80, 2.11, 2.90 },
    { 3.09, 1.85, 3.05, 2.72 },
    { 3.96, 3.24, 2.70, 3.60 },
  }
};

const float maxdist[TOPOMAX][4][4] =
{
  {
    { 4.11, 5.27, 6.18, 2.86 },
    { 5.11, 6.42, 7.40, 3.92 },
    { 3.68, 4.92, 5.81, 2.57 },
    { 5.25, 6.40, 7.24, 3.97 },
  },
  {
    { 7.71, 8.60, 9.55, 6.37 },
    { 8.63, 9.80, 10.78, 7.38 },
    { 7.06, 8.13, 9.03, 5.84 },
    { 8.87, 9.75, 10.64, 7.46 },
  },
  {
    { 11.00, 12.22, 12.86, 9.80 },
    { 12.15, 13.23, 13.86, 10.91 },
    { 10.46, 11.48, 11.98, 9.16 },
    { 12.10, 13.28, 14.00, 10.89 },
  },
  {
    { 14.52, 15.42, 16.19, 13.15 },
    { 15.52, 16.53, 17.22, 14.21 },
    { 13.65, 14.75, 15.29, 12.52 },
    { 15.60, 16.53, 17.41, 14.33 },
  },
  {
    { 17.75, 18.79, 19.68, 16.50 },
    { 18.96, 20.08, 20.66, 17.63 },
    { 16.98, 18.00, 18.86, 15.73 },
    { 19.00, 20.06, 20.78, 17.68 },
  },
  {
    { 21.15, 22.38, 23.02, 19.93 },
    { 22.09, 23.00, 23.59, 20.83 },
    { 20.18, 21.25, 21.92, 18.91 },
    { 22.32, 23.55, 24.07, 21.09 },
  },
  {
    { 24.59, 25.49, 26.34, 23.27 },
    { 25.54, 26.60, 27.03, 24.28 },
    { 23.59, 24.67, 25.39, 22.42 },
    { 25.69, 26.81, 27.46, 24.39 },
  },
  {
    { 27.98, 28.73, 29.73, 26.64 },
    { 28.76, 29.69, 30.34, 27.52 },
    { 26.89, 27.92, 28.55, 25.51 },
    { 28.97, 29.96, 30.88, 27.75 },
  },
  {
    { 31.29, 32.26, 33.10, 29.99 },
    { 32.35, 33.51, 33.44, 31.08 },
    { 30.22, 31.15, 32.11, 29.00 },
    { 32.38, 33.34, 34.30, 31.19 },
  },
  {
    { 34.39, 35.43, 35.84, 33.20 },
    { 34.98, 35.86, 36.51, 34.06 },
    { 33.60, 34.57, 34.84, 32.36 },
    { 35.76, 36.78, 36.96, 34.55 },
  },
  {
    { 37.75, 38.83, 39.41, 36.52 },
    { 38.74, 39.60, 40.25, 37.35 },
    { 36.54, 37.68, 38.52, 35.32 },
    { 38.50, 39.55, 39.50, 37.26 },
  }
};

/* DSSP residue solvent accessible area in GGXGG extended pentapeptide */
const float     resacc[22] =
{
    113.0, 253.0, 167.0, 167.0, 140.0, 199.0, 198.0, 88.0, 194.0, 178.0,
    179.0, 215.0, 194.0, 226.0, 151.0, 134.0, 148.0, 268.0, 242.0, 157.0,
    88.0, 88.0
};

/* Amino acid molecular weights */
const float     molwt[22] =
{
    89.09, 174.20, 132.12, 133.10, 121.15, 146.15, 147.13, 75.07, 155.16, 131.17,
    131.17, 146.19, 149.21, 165.19, 115.13, 105.09, 119.12, 204.24, 181.19, 117.15,
    89.09, 89.09
};

/* Amino acid radii */
const float     aaradii[22] =
{
    0.77, 2.38, 1.45, 1.43, 1.22, 1.75, 1.77, 0.58, 1.78, 1.56,
    1.54, 2.08, 1.80, 1.90, 1.25, 1.08, 1.24, 2.21, 2.13, 1.29,
    0.77, 0.77
};

/* Matrix of pairwise minimum CB-CB distances */
const float cbmin[20][20] =
{
    {   3.99, 3.83, 3.62, 3.76, 3.83, 3.83, 3.65, 2.43, 4.46, 3.85, 3.89, 3.84, 3.72, 4.26, 3.53, 3.61, 3.52, 4.71, 3.61, 3.79 },
    {	3.83, 4.84, 4.13, 3.78, 4.97, 4.35, 3.86, 2.61, 4.26, 4.28, 3.90, 4.18, 5.51, 3.83, 4.03, 4.19, 4.07, 4.86, 3.79, 4.90 },
    {	3.62, 4.13, 4.53, 3.77, 4.19, 4.34, 4.18, 2.81, 4.35, 5.29, 3.67, 4.31, 4.43, 4.25, 3.65, 4.04, 4.09, 4.18, 3.85, 4.14 },
    {	3.76, 3.78, 3.77, 4.02, 3.76, 3.90, 4.54, 2.71, 4.10, 4.36, 3.88, 4.13, 4.36, 3.84, 3.70, 3.63, 3.85, 3.47, 4.13, 4.07 },
    {	3.83, 4.97, 4.19, 3.76, 3.51, 4.18, 5.02, 2.34, 4.73, 4.20, 4.53, 5.38, 3.94, 4.30, 4.18, 3.93, 3.97, 5.54, 4.00, 4.49 },
    {	3.83, 4.35, 4.34, 3.90, 4.18, 3.78, 3.68, 2.95, 4.17, 4.52, 4.31, 3.92, 4.20, 4.20, 4.20, 4.09, 3.87, 4.89, 4.28, 4.41 },
    {	3.65, 3.86, 4.18, 4.54, 5.02, 3.68, 4.00, 3.43, 3.89, 5.04, 4.06, 4.15, 5.49, 4.08, 3.48, 3.87, 4.50, 4.12, 3.86, 4.19 },
    {	2.43, 2.61, 2.81, 2.71, 2.34, 2.95, 3.43, 1.92, 3.29, 3.61, 3.24, 3.00, 4.17, 3.26, 2.85, 2.53, 2.61, 3.32, 2.87, 3.32 },
    {	4.46, 4.26, 4.35, 4.10, 4.73, 4.17, 3.89, 3.29, 5.07, 5.32, 3.70, 3.96, 5.23, 4.37, 3.86, 4.07, 4.34, 5.27, 5.03, 4.30 },
    {	3.85, 4.28, 5.29, 4.36, 4.20, 4.52, 5.04, 3.61, 5.32, 4.04, 3.89, 4.17, 4.11, 4.65, 3.97, 4.06, 4.71, 4.45, 4.74, 4.26 },
    {	3.89, 3.90, 3.67, 3.88, 4.53, 4.31, 4.06, 3.24, 3.70, 3.89, 4.35, 3.83, 3.85, 4.21, 4.15, 3.59, 4.29, 3.81, 3.92, 4.26 },
    {	3.84, 4.18, 4.31, 4.13, 5.38, 3.92, 4.15, 3.00, 3.96, 4.17, 3.83, 4.24, 5.20, 4.48, 4.11, 3.62, 4.03, 4.58, 3.75, 4.19 },
    {	3.72, 5.51, 4.43, 4.36, 3.94, 4.20, 5.49, 4.17, 5.23, 4.11, 3.85, 5.20, 6.40, 4.10, 4.20, 4.53, 3.92, 5.36, 3.99, 4.07 },
    {	4.26, 3.83, 4.25, 3.84, 4.30, 4.20, 4.08, 3.26, 4.37, 4.65, 4.21, 4.48, 4.10, 4.13, 3.96, 3.87, 4.17, 4.00, 4.68, 4.07 },
    {	3.53, 4.03, 3.65, 3.70, 4.18, 4.20, 3.48, 2.85, 3.86, 3.97, 4.15, 4.11, 4.20, 3.96, 3.74, 3.85, 3.97, 4.40, 3.94, 4.26 },
    {	3.61, 4.19, 4.04, 3.63, 3.93, 4.09, 3.87, 2.53, 4.07, 4.06, 3.59, 3.62, 4.53, 3.87, 3.85, 4.21, 3.76, 5.51, 3.61, 3.47 },
    {	3.52, 4.07, 4.09, 3.85, 3.97, 3.87, 4.50, 2.61, 4.34, 4.71, 4.29, 4.03, 3.92, 4.17, 3.97, 3.76, 4.19, 4.72, 4.32, 4.31 },
    {	4.71, 4.86, 4.18, 3.47, 5.54, 4.89, 4.12, 3.32, 5.27, 4.45, 3.81, 4.58, 5.36, 4.00, 4.40, 5.51, 4.72, 6.43, 3.74, 4.93 },
    {	3.61, 3.79, 3.85, 4.13, 4.00, 4.28, 3.86, 2.87, 5.03, 4.74, 3.92, 3.75, 3.99, 4.68, 3.94, 3.61, 4.32, 3.74, 3.89, 4.13 },
    {	3.79, 4.90, 4.14, 4.07, 4.49, 4.41, 4.19, 3.32, 4.30, 4.26, 4.26, 4.19, 4.07, 4.07, 4.26, 3.47, 4.31, 4.93, 4.13, 4.34 }
};

/* Matrix of pairwise minimum side-chain centroid-centroid distances */
const float cgmin[20][20] =
{
    { 3.62, 3.93, 3.77, 3.72, 4.10, 4.03, 3.88, 3.70, 3.98, 4.02, 4.06, 4.18, 4.15, 3.77, 4.06, 3.68, 3.85, 3.84, 3.68, 3.85 },
    { 3.93, 4.40, 4.29, 3.88, 5.07, 4.38, 4.11, 4.15, 4.35, 4.48, 4.35, 4.81, 4.82, 4.31, 4.42, 3.93, 4.15, 4.30, 4.07, 4.32 },
    { 3.77, 4.29, 4.13, 4.03, 4.72, 4.34, 4.29, 3.91, 4.47, 4.49, 4.38, 4.32, 4.70, 4.45, 4.38, 3.92, 4.07, 4.28, 4.29, 4.30 },
    { 3.72, 3.88, 4.03, 4.06, 4.84, 4.36, 4.42, 3.79, 4.36, 4.48, 4.38, 3.95, 4.84, 4.68, 4.28, 3.58, 3.82, 5.32, 4.89, 4.42 },
    { 4.10, 5.07, 4.72, 4.84, 2.68, 4.97, 4.79, 4.13, 5.05, 4.55, 4.58, 5.42, 4.82, 4.40, 4.58, 4.50, 4.46, 5.51, 4.73, 4.34 },
    { 4.03, 4.38, 4.34, 4.36, 4.97, 4.38, 4.53, 4.15, 4.62, 4.65, 4.53, 4.45, 4.83, 4.35, 4.42, 3.98, 4.19, 4.58, 4.30, 4.47 },
    { 3.88, 4.11, 4.29, 4.42, 4.79, 4.53, 4.42, 4.05, 4.56, 4.49, 4.44, 4.05, 4.97, 4.49, 4.30, 3.78, 3.85, 5.00, 4.53, 4.26 },
    { 3.70, 4.15, 3.91, 3.79, 4.13, 4.15, 4.05, 3.82, 4.21, 4.21, 4.23, 4.24, 4.25, 3.73, 4.18, 3.59, 3.82, 3.95, 3.74, 4.06 },
    { 3.98, 4.35, 4.47, 4.36, 5.05, 4.62, 4.56, 4.21, 4.68, 4.55, 4.42, 4.71, 4.97, 4.81, 4.42, 4.17, 4.36, 4.86, 4.10, 4.44 },
    { 4.02, 4.48, 4.49, 4.48, 4.55, 4.65, 4.49, 4.21, 4.55, 4.43, 4.51, 4.81, 4.59, 4.17, 4.61, 4.30, 4.39, 4.47, 4.18, 4.34 },
    { 4.06, 4.35, 4.38, 4.38, 4.58, 4.53, 4.44, 4.23, 4.42, 4.51, 4.54, 4.65, 4.55, 4.28, 4.60, 4.10, 4.35, 4.23, 4.20, 4.36 },
    { 4.18, 4.81, 4.32, 3.95, 5.42, 4.45, 4.05, 4.24, 4.71, 4.81, 4.65, 4.88, 4.99, 4.32, 4.84, 4.18, 4.36, 4.35, 4.19, 4.67 },
    { 4.15, 4.82, 4.70, 4.84, 4.82, 4.83, 4.97, 4.25, 4.97, 4.59, 4.55, 4.99, 4.85, 4.41, 4.84, 4.40, 4.60, 4.98, 4.43, 4.50 },
    { 3.77, 4.31, 4.45, 4.68, 4.40, 4.35, 4.49, 3.73, 4.81, 4.17, 4.28, 4.32, 4.41, 4.46, 4.14, 4.51, 4.36, 5.03, 4.60, 4.07 },
    { 4.06, 4.42, 4.38, 4.28, 4.58, 4.42, 4.30, 4.18, 4.42, 4.61, 4.60, 4.84, 4.84, 4.14, 4.55, 4.11, 4.32, 4.26, 4.00, 4.45 },
    { 3.68, 3.93, 3.92, 3.58, 4.50, 3.98, 3.78, 3.59, 4.17, 4.30, 4.10, 4.18, 4.40, 4.51, 4.11, 3.42, 3.69, 4.61, 4.22, 4.01 },
    { 3.85, 4.15, 4.07, 3.82, 4.46, 4.19, 3.85, 3.82, 4.36, 4.39, 4.35, 4.36, 4.60, 4.36, 4.32, 3.69, 3.87, 5.01, 4.38, 4.26 },
    { 3.84, 4.30, 4.28, 5.32, 5.51, 4.58, 5.00, 3.95, 4.86, 4.47, 4.23, 4.35, 4.98, 5.03, 4.26, 4.61, 5.01, 5.20, 5.12, 4.24 },
    { 3.68, 4.07, 4.29, 4.89, 4.73, 4.30, 4.53, 3.74, 4.10, 4.18, 4.20, 4.19, 4.43, 4.60, 4.00, 4.22, 4.38, 5.12, 4.68, 4.10 },
    { 3.85, 4.32, 4.30, 4.42, 4.34, 4.47, 4.26, 4.06, 4.44, 4.34, 4.36, 4.67, 4.50, 4.07, 4.45, 4.01, 4.26, 4.24, 4.10, 4.21 }
};

const int isacceptor[20][14] =
{
    { 0,0,0,1,0,0,0,0,0,0,0,0,0,0 },
    { 0,0,0,1,0,0,0,0,0,0,0,0,0,0 },
    { 0,0,0,1,0,0,1,0,0,0,0,0,0,0 },
    { 0,0,0,1,0,0,1,1,0,0,0,0,0,0 },
    { 0,0,0,1,0,1,0,0,0,0,0,0,0,0 },
    { 0,0,0,1,0,0,0,1,0,0,0,0,0,0 },
    { 0,0,0,1,0,0,0,1,1,0,0,0,0,0 },
    { 0,0,0,1,0,0,0,0,0,0,0,0,0,0 },
    { 0,0,0,1,0,0,1,0,0,0,0,0,0,0 },
    { 0,0,0,1,0,0,0,0,0,0,0,0,0,0 },
    { 0,0,0,1,0,0,0,0,0,0,0,0,0,0 },
    { 0,0,0,1,0,0,0,0,0,0,0,0,0,0 },
    { 0,0,0,1,0,0,0,0,0,0,0,0,0,0 },
    { 0,0,0,1,0,0,0,0,0,0,0,0,0,0 },
    { 0,0,0,1,0,0,0,0,0,0,0,0,0,0 },
    { 0,0,0,1,0,1,0,0,0,0,0,0,0,0 },
    { 0,0,0,1,0,1,0,0,0,0,0,0,0,0 },
    { 0,0,0,1,0,0,0,0,0,0,0,0,0,0 },
    { 0,0,0,1,0,0,0,0,0,0,0,1,0,0 },
    { 0,0,0,1,0,0,0,0,0,0,0,0,0,0 }
};

const int isdonor[20][14] =
{
    { 1,0,0,0,0,0,0,0,0,0,0,0,0,0 },
    { 1,0,0,0,0,0,0,1,0,1,1,0,0,0 },
    { 1,0,0,0,0,0,0,1,0,0,0,0,0,0 },
    { 1,0,0,0,0,0,0,0,0,0,0,0,0,0 },
    { 1,0,0,0,0,1,0,0,0,0,0,0,0,0 },
    { 1,0,0,0,0,0,0,0,1,0,0,0,0,0 },
    { 1,0,0,0,0,0,0,0,0,0,0,0,0,0 },
    { 1,0,0,0,0,0,0,0,0,0,0,0,0,0 },
    { 1,0,0,0,0,0,1,0,0,1,0,0,0,0 },
    { 1,0,0,0,0,0,0,0,0,0,0,0,0,0 },
    { 1,0,0,0,0,0,0,0,0,0,0,0,0,0 },
    { 1,0,0,0,0,0,0,0,1,0,0,0,0,0 },
    { 1,0,0,0,0,0,0,0,0,0,0,0,0,0 },
    { 1,0,0,0,0,0,0,0,0,0,0,0,0,0 },
    { 0,0,0,0,0,0,0,0,0,0,0,0,0,0 },
    { 1,0,0,0,0,1,0,0,0,0,0,0,0,0 },
    { 1,0,0,0,0,1,0,0,0,0,0,0,0,0 },
    { 1,0,0,0,0,0,0,0,1,0,0,0,0,0 },
    { 1,0,0,0,0,0,0,0,0,0,0,1,0,0 },
    { 1,0,0,0,0,0,0,0,0,0,0,0,0,0 }
};

/* A->AA covalent partners */
const int hbaa_table[20][14] =
{
    { 0,0,0,2,0,0,0,0,0,0,0,0,0,0 },
    { 0,0,0,2,0,0,0,0,0,0,0,0,0,0 },
    { 0,0,0,2,0,0,5,0,0,0,0,0,0,0 },
    { 0,0,0,2,0,0,5,5,0,0,0,0,0,0 },
    { 0,0,0,2,0,4,0,0,0,0,0,0,0,0 },
    { 0,0,0,2,0,0,0,6,0,0,0,0,0,0 },
    { 0,0,0,2,0,0,0,6,6,0,0,0,0,0 },
    { 0,0,0,2,0,0,0,0,0,0,0,0,0,0 },
    { 0,0,0,2,0,0,5,0,0,0,0,0,0,0 },
    { 0,0,0,2,0,0,0,0,0,0,0,0,0,0 },
    { 0,0,0,2,0,0,0,0,0,0,0,0,0,0 },
    { 0,0,0,2,0,0,0,0,0,0,0,0,0,0 },
    { 0,0,0,2,0,0,0,0,0,0,0,0,0,0 },
    { 0,0,0,2,0,0,0,0,0,0,0,0,0,0 },
    { 0,0,0,2,0,0,0,0,0,0,0,0,0,0 },
    { 0,0,0,2,0,4,0,0,0,0,0,0,0,0 },
    { 0,0,0,2,0,4,0,0,0,0,0,0,0,0 },
    { 0,0,0,2,0,0,0,0,0,0,0,0,0,0 },
    { 0,0,0,2,0,0,0,0,0,0,0,10,0,0 },
    { 0,0,0,2,0,0,0,0,0,0,0,0,0,0 }
};

/* D->DD covalent partners */
const int hbdd_table[20][14] =
{
    { 1,0,0,0,0,0,0,0,0,0,0,0,0,0 },
    { 1,0,0,0,0,0,0,6,0,8,8,0,0,0 },
    { 1,0,0,0,0,0,0,5,0,0,0,0,0,0 },
    { 1,0,0,0,0,0,0,0,0,0,0,0,0,0 },
    { 1,0,0,0,0,4,0,0,0,0,0,0,0,0 },
    { 1,0,0,0,0,0,0,0,6,0,0,0,0,0 },
    { 1,0,0,0,0,0,0,0,0,0,0,0,0,0 },
    { 1,0,0,0,0,0,0,0,0,0,0,0,0,0 },
    { 1,0,0,0,0,0,5,0,0,7,0,0,0,0 },
    { 1,0,0,0,0,0,0,0,0,0,0,0,0,0 },
    { 1,0,0,0,0,0,0,0,0,0,0,0,0,0 },
    { 1,0,0,0,0,0,0,0,7,0,0,0,0,0 },
    { 1,0,0,0,0,0,0,0,0,0,0,0,0,0 },
    { 1,0,0,0,0,0,0,0,0,0,0,0,0,0 },
    { 1,0,0,0,0,0,0,0,0,0,0,0,0,0 },
    { 1,0,0,0,0,4,0,0,0,0,0,0,0,0 },
    { 1,0,0,0,0,4,0,0,0,0,0,0,0,0 },
    { 1,0,0,0,0,0,0,0,6,0,0,0,0,0 },
    { 1,0,0,0,0,0,0,0,0,0,0,10,0,0 },
    { 1,0,0,0,0,0,0,0,0,0,0,0,0,0 }
};


/*  RMSD Correlation Substitution Matrix */
const short           rmsdcmat[23][23] =
{
    { 10,  7,  8,  0,  3,  9,  0,  9,  9,  7,  0,  9,  0,  0,  0,  9,  9,  0,  0,  9,  4,  4,  6 },
    {  7, 10,  9,  0,  0,  9,  0,  9,  9,  9,  0,  8,  9,  6,  0,  7,  9,  0,  7,  9,  4,  4,  0 },
    {  8,  9, 10,  9,  9,  9,  7,  9,  9,  9,  0,  9,  8,  9,  0,  0,  1,  0,  9,  0,  9,  8,  9 },
    {  0,  0,  9, 10,  9,  0,  9,  0,  0,  0,  0,  0,  9,  0,  0,  9,  0,  0,  0,  8,  9,  4,  9 },
    {  3,  0,  9,  9, 10,  0,  9,  9,  7,  0,  9,  8,  0,  9,  9,  0,  9,  7,  9,  9,  9,  4,  0 },
    {  9,  9,  9,  0,  0, 10,  8,  8,  0,  0,  0,  9,  9,  9,  9,  9,  9,  3,  0,  0,  4,  9,  7 },
    {  0,  0,  7,  9,  9,  8, 10,  0,  0,  9,  0,  0,  0,  0,  9,  9,  9,  9,  0,  9,  8,  9,  9 },
    {  9,  9,  9,  0,  9,  8,  0, 10,  1,  0,  0,  0,  0,  0,  9,  0,  0,  0,  0,  0,  4,  4,  0 },
    {  9,  9,  9,  0,  7,  0,  0,  1, 10,  9,  9,  9,  0,  9,  0,  8,  9,  0,  9,  8,  4,  0,  6 },
    {  7,  9,  9,  0,  0,  0,  9,  0,  9, 10,  8,  0,  8,  9,  0,  0,  0,  0,  0,  9,  4,  4,  0 },
    {  0,  0,  0,  0,  9,  0,  0,  0,  9,  8, 10,  0,  3,  9,  0,  0,  1,  0,  8,  8,  0,  0,  7 },
    {  9,  8,  9,  0,  8,  9,  0,  0,  9,  0,  0, 10,  9,  9,  9,  0,  9,  0,  0,  9,  4,  4,  1 },
    {  0,  9,  8,  9,  0,  9,  0,  0,  0,  8,  3,  9, 10,  9,  0,  0,  0,  9,  0,  0,  8,  4,  0 },
    {  0,  6,  9,  0,  9,  9,  0,  0,  9,  9,  9,  9,  9, 10,  7,  0,  0,  9,  9,  7,  4,  4,  9 },
    {  0,  0,  0,  0,  9,  9,  9,  9,  0,  0,  0,  9,  0,  7, 10,  9,  9,  0,  0,  9,  0,  9,  0 },
    {  9,  7,  0,  9,  0,  9,  9,  0,  8,  0,  0,  0,  0,  0,  9, 10,  9,  1,  0,  9,  4,  9,  0 },
    {  9,  9,  1,  0,  9,  9,  9,  0,  9,  0,  1,  9,  0,  0,  9,  9, 10,  0,  0,  9,  0,  9,  0 },
    {  0,  0,  0,  0,  7,  3,  9,  0,  0,  0,  0,  0,  9,  9,  0,  1,  0, 10,  9,  8,  0,  6,  9 },
    {  0,  7,  9,  0,  9,  0,  0,  0,  9,  0,  8,  0,  0,  9,  0,  0,  0,  9, 10,  0,  4,  0,  0 },
    {  9,  9,  0,  8,  9,  0,  9,  0,  8,  9,  8,  9,  0,  7,  9,  9,  9,  8,  0, 10,  4,  4,  9 },
    {  4,  4,  9,  9,  9,  4,  8,  4,  4,  4,  0,  4,  8,  4,  0,  4,  0,  0,  4,  4, 10,  6,  0 },
    {  4,  4,  8,  4,  4,  9,  9,  4,  0,  4,  0,  4,  4,  4,  9,  9,  9,  6,  0,  4,  6, 10,  0 },
    {  6,  0,  9,  9,  0,  7,  9,  0,  6,  0,  7,  1,  0,  9,  0,  0,  0,  9,  0,  9,  0,  0, 10 }
};

/* Dirichlet Mixture Model Substitution Matrix (Crooks-Brenner, 2005) */
const short           dmmmat[23][23] =
{
    {  4, -2, -2, -2,  0, -1, -1, -1, -2, -2, -2, -2, -1, -2, -1,  0, -1, -3, -2, -1, -2, -1, -1 },
    { -2,  6, -1, -1, -4,  1,  0, -3,  0, -3, -3,  2, -2, -3, -2, -1, -1, -2, -2, -3, -1,  0, -1 },
    { -2, -1,  7,  1, -3,  0,  0, -1,  0, -5, -4,  0, -3, -4, -2,  0, -1, -3, -2, -4,  4,  0, -1 },
    { -2, -1,  1,  7, -4,  0,  1, -1, -1, -6, -5,  0, -4, -5, -1,  0, -1, -4, -3, -5,  4,  0, -1 },
    {  0, -4, -3, -4, 12, -3, -4, -3, -3, -1, -2, -4, -1, -2, -3, -2, -1, -2, -2,  0, -3, -3, -1 },
    { -1,  1,  0,  0, -3,  6,  1, -2,  0, -3, -3,  1, -2, -3, -1,  0, -1, -2, -2, -3,  0,  3,  0 },
    { -1,  0,  0,  1, -4,  1,  5, -2, -1, -4, -4,  1, -3, -4, -1, -1, -1, -3, -3, -4,  0,  3, -1 },
    { -1, -3, -1, -1, -3, -2, -2,  7, -2, -6, -5, -2, -4, -5, -2, -1, -2, -4, -4, -5, -1, -2, -2 },
    { -2,  0,  0, -1, -3,  0, -1, -2,  9, -3, -3, -1, -2, -1, -2, -1, -1,  0,  0, -3,  0,  0,  0 },
    { -2, -3, -5, -6, -1, -3, -4, -6, -3,  5,  2, -4,  1,  0, -4, -4, -2, -1, -1,  3, -5, -3, -1 },
    { -2, -3, -4, -5, -2, -3, -4, -5, -3,  2,  5, -3,  2,  1, -3, -3, -2, -1, -1,  1, -4, -3, -1 },
    { -2,  2,  0,  0, -4,  1,  1, -2, -1, -4, -3,  5, -2, -4, -1, -1, -1, -3, -3, -3,  0,  1, -1 },
    { -1, -2, -3, -4, -1, -2, -3, -4, -2,  1,  2, -2,  7,  1, -3, -2, -1,  0,  0,  1, -3, -2,  0 },
    { -2, -3, -4, -5, -2, -3, -4, -5, -1,  0,  1, -4,  1,  7, -3, -3, -2,  3,  3,  0, -4, -3, -1 },
    { -1, -2, -2, -1, -3, -1, -1, -2, -2, -4, -3, -1, -3, -3,  8, -1, -2, -3, -3, -3, -1, -1, -1 },
    {  0, -1,  0,  0, -2,  0, -1, -1, -1, -4, -3, -1, -2, -3, -1,  4,  1, -3, -2, -3,  0,  0, -1 },
    { -1, -1, -1, -1, -1, -1, -1, -2, -1, -2, -2, -1, -1, -2, -2,  1,  5, -2, -2, -1, -1, -1,  0 },
    { -3, -2, -3, -4, -2, -2, -3, -4,  0, -1, -1, -3,  0,  3, -3, -3, -2, 12,  3, -2, -3, -2, -1 },
    { -2, -2, -2, -3, -2, -2, -3, -4,  0, -1, -1, -3,  0,  3, -3, -2, -2,  3,  8, -2, -2, -2, -1 },
    { -1, -3, -4, -5,  0, -3, -4, -5, -3,  3,  1, -3,  1,  0, -3, -3, -1, -2, -2,  5, -4, -3, -1 },
    { -2, -1,  4,  4, -3,  0,  0, -1,  0, -5, -4,  0, -3, -4, -1,  0, -1, -3, -2, -4,  5,  0,  0 },
    { -1,  0,  0,  0, -3,  3,  3, -2,  0, -3, -3,  1, -2, -3, -1,  0, -1, -2, -2, -3,  0,  4,  0 },
    { -1, -1, -1, -1, -1,  0, -1, -2,  0, -1, -1, -1,  0, -1, -1, -1,  0, -1, -1, -1,  0,  0, -1 }
};

/* Maximum CB-CB distance tables (> 3.5 s.d.) */
const float condmax[20][20] = {
    {  7.4, 8.8, 8.1, 8.8, 8.3, 8.3, 9.8, 7.4, 8.0, 8.5, 9.4, 9.1,10.0,10.0, 7.7, 7.2, 9.0, 8.8,10.0, 9.9 },
    {  8.8,10.0,10.0,10.0,10.0, 9.8,10.0, 8.8, 9.4, 8.0, 9.5,10.0,10.0, 9.9, 8.7, 9.5,10.0, 9.8,10.0, 9.1 },
    {  8.1,10.0, 8.0, 8.1,10.0, 8.5, 9.5, 8.1, 8.1, 9.2, 7.7, 8.7, 8.6,10.0, 7.4, 9.1, 7.5,10.0,10.0, 7.6 },
    {  8.8,10.0, 8.1, 8.8,10.0, 9.1, 9.4, 8.8, 9.1, 8.6, 9.8,10.0, 9.5, 8.8, 7.4, 7.6, 9.1,10.0, 9.8, 9.9 },
    {  8.3,10.0,10.0,10.0, 6.4, 7.9,10.0, 8.3,10.0, 7.9, 8.8,10.0,10.0, 8.7,10.0,10.0,10.0,10.0,10.0, 7.0 },
    {  8.3, 9.8, 8.5, 9.1, 7.9, 7.9, 9.7, 8.3, 9.9, 9.6,10.0, 9.7, 9.1,10.0, 9.7, 9.9, 9.2,10.0,10.0, 9.7 },
    {  9.8,10.0, 9.5, 9.4,10.0, 9.7, 9.5, 9.8, 9.3, 9.7, 9.3,10.0, 9.6, 9.7, 7.3, 8.9, 9.9, 9.9, 9.9, 8.8 },
    {  7.4, 8.8, 8.1, 8.8, 8.3, 8.3, 9.8, 7.4, 8.0, 8.5, 9.4, 9.1,10.0,10.0, 7.7, 7.2, 9.0, 8.8,10.0, 9.9 },
    {  8.0, 9.4, 8.1, 9.1,10.0, 9.9, 9.3, 8.0, 9.1, 8.9, 9.3, 8.8, 9.0,10.0,10.0, 9.8, 9.4,10.0, 9.4, 8.0 },
    {  8.5, 8.0, 9.2, 8.6, 7.9, 9.6, 9.7, 8.5, 8.9, 9.0, 8.7, 8.4, 9.3, 9.6, 8.0, 9.8, 9.3, 9.7,10.0, 8.7 },
    {  9.4, 9.5, 7.7, 9.8, 8.8,10.0, 9.3, 9.4, 9.3, 8.7, 9.2, 9.6, 9.5,10.0, 9.7, 8.6, 9.1,10.0, 9.6, 8.6 },
    {  9.1,10.0, 8.7,10.0,10.0, 9.7,10.0, 9.1, 8.8, 8.4, 9.6,10.0, 9.0, 9.1, 8.3,10.0, 8.9, 7.7,10.0, 9.5 },
    { 10.0,10.0, 8.6, 9.5,10.0, 9.1, 9.6,10.0, 9.0, 9.3, 9.5, 9.0,10.0,10.0,10.0, 8.2, 9.5, 8.6, 8.7, 9.0 },
    { 10.0, 9.9,10.0, 8.8, 8.7,10.0, 9.7,10.0,10.0, 9.6,10.0, 9.1,10.0,10.0, 9.2, 8.7, 9.9,10.0,10.0, 9.4 },
    {  7.7, 8.7, 7.4, 7.4,10.0, 9.7, 7.3, 7.7,10.0, 8.0, 9.7, 8.3,10.0, 9.2,10.0, 7.2, 8.4,10.0,10.0, 9.3 },
    {  7.2, 9.5, 9.1, 7.6,10.0, 9.9, 8.9, 7.2, 9.8, 9.8, 8.6,10.0, 8.2, 8.7, 7.2, 7.6, 8.3, 9.0, 9.1, 8.1 },
    {  9.0,10.0, 7.5, 9.1,10.0, 9.2, 9.9, 9.0, 9.4, 9.3, 9.1, 8.9, 9.5, 9.9, 8.4, 8.3, 8.1, 9.5,10.0, 9.6 },
    {  8.8, 9.8,10.0,10.0,10.0,10.0, 9.9, 8.8,10.0, 9.7,10.0, 7.7, 8.6,10.0,10.0, 9.0, 9.5,10.0,10.0, 9.8 },
    { 10.0,10.0,10.0, 9.8,10.0,10.0, 9.9,10.0, 9.4,10.0, 9.6,10.0, 8.7,10.0,10.0, 9.1,10.0,10.0, 8.5, 9.1 },
    {  9.9, 9.1, 7.6, 9.9, 7.0, 9.7, 8.8, 9.9, 8.0, 8.7, 8.6, 9.5, 9.0, 9.4, 9.3, 8.1, 9.6, 9.8, 9.1, 8.3 }
};

typedef float  Transform[4][4];
typedef float  Vector[3];

typedef struct
{
    char            *atmnam;
    short           type;
    Vector          pos;
}
SIDECATM;

typedef struct
{
    Vector           n, h, ca, c, o, cb, sc_cg;
    float            phi, psi, omega;
    short            nscats, rotnum, builtflg;
    SIDECATM        *sc;
}
RESATM;

const int nscatoms[20] = {
    1, 7, 4, 4, 2, 5, 5, 0, 6, 4, 4, 5, 4, 7, 3, 2, 3, 10, 8, 3
};

/* Side chain rotamer table */

int nrots[20];

SIDECATM rotamers[20][200][10];
float rotafreq[20][200];

struct chnentry
{
    short           length, *relacc;
    char           *seq, *sstruc;
    RESATM         *chain;
} chnlist[20000];

struct chnidx
{
    short           chain, pos, length;
    float           weight, rmsd;
};

struct flist
{
    struct chnidx  *frags;
    int             nfrags, ntot;
    short           nchn, from, bestlen;
    float           totwt, totbest, rmsd, totsum, totsumsq;
    char           *type;
}
fraglist[MAXSEQLEN + 1], fraglist2[MAXSEQLEN + 1];

struct SSTRUC
{
    short           start, length, type;
}
sstlist[100];

int      nchn, chnres;

float **cb_targ, **cb_last;

FILE    *rfp;

typedef struct
{
    RESATM         *genome;
    float           perfval, selval;
    short           evalflg;
}
Schema;

typedef struct
{
    RESATM         *chn;
    float           t, energy;
}
Replica;

float    diffdist[TOPOMAX][4][4];

const char     *rescodes = "ARNDCQEGHILKMFPSTWYVX";

/* Offsets for residue-specific atom types */
const int resatofs[20] = {
    0, 5, 16, 24, 32, 38, 47, 56, 60, 70, 78, 86, 95, 103, 114, 121, 127, 134, 148, 160
};

/* Table of minimum distances between 167 atom types */
float atomdistsq[167][167][3];

/* DFIRE potential table */
float dfire[167][167][20];

/* Local conformation tables */
float    sr_de[4][4][21][21][INTERVALS][TOPOMAX];

/* LR Residue-Residue interaction matrices */
float    lr_de[4][4][21][21][NRR];

/* Residue accessibility matrices */
float    acc_de[21][NACC];

/* Residue OOI number matrices */
float    ooi_de[21][(OOIMAX - OOIMIN) / OOIDIV + 1];

/* Structure generation arrays */

float    (**dsqmat)[NPAIRS];

float    ***cbcb_potmat, **solv_potmat, **conmat, **dcomat, **hbmat;

RESATM  *curchn, *bestchn, *oldchn, *targchn;

Schema  *curpool, *newpool;
int      poolsize = 5000, seqlen, *samparr, besti;
float    mutrate = 0.5, crosrate = 0.1, maxngen = 200, worst, best, avc_perf;

/* Flags for program modes */
int mod_mode, opt_mode, wt_mode = 1, global_ster_mode = 3, fragsel;


/* Dump a rude message to standard error and exit */
void fail(char *fmt, ...)
{
    va_list ap;
    
    va_start(ap, fmt) ;
    fprintf(stderr, "*** ");
    vfprintf(stderr, fmt, ap);
    fputc('\n', stderr);
    
    exit(-1);
}


/* Return accessibility range */
int
                accrange(int acc)
{
    if (acc < 12)
	return 0;
    else if (acc < 36)
	return 1;
    else if (acc < 44)
	return 2;
    else if (acc < 87)
	return 3;
    return 4;
}


/* Convert AA letter to numeric code (0-22) */
int
                aanum(int ch)
{
    static const int      aacvs[] =
    {
	999, 0, 20, 4, 3, 6, 13, 7, 8, 9, 22, 11, 10, 12, 2,
	22, 14, 5, 1, 15, 16, 22, 19, 17, 22, 18, 21
    };

    return isalpha(ch) ? aacvs[ch & 31] : 22;
}


/* Allocate matrix */
void           *allocmat(int rows, int columns, int size)
{
    int             i;
    void          **p;

    p = malloc(rows * sizeof(void *));

    if (p == NULL)
	fail("allocmat: malloc [] failed!");
    for (i = 0; i < rows; i++)
	if ((p[i] = calloc(columns, size)) == NULL)
	    fail("allocmat: malloc [][] failed!");

    return p;
}


/* Free matrix */
void
                freemat(void *p, int rows)
{
    int             i;

    for (i = 0; i < rows; i++)
	free(((void **) p)[i]);
    free(p);
}


/* Implementation of the WELL1024 PRNG by Panneton, L'Ecuyer & Matsumoto (Period 2^1024-1) */

unsigned int wstate[32], widx;

unsigned int WELL1024(void)
{
    unsigned int a, b, c;

    a = wstate[(widx+31) & 31];
    b = wstate[widx] ^ wstate[(widx+3) & 31] ^ (wstate[(widx+3) & 31]>>8);
    c = wstate[(widx+24) & 31] ^ (wstate[(widx+24) & 31]<<19) ^ wstate[(widx+10) & 31] ^ (wstate[(widx+10) & 31]<<14);

    wstate[widx] = b ^ c;
    widx = (widx + 31) & 31;

    return wstate[widx] = (a^(a<<11)) ^ (b^(b<<7)) ^ (c^(c<<13));
}


/* Generate random hash of integer */
unsigned int lcghash(unsigned int key)
{
    int i;
    
    for (i=0; i<4; i++)
    {
	key = 314527869 * key + 1;
	key ^= key >> 22;
    }

    return key;
}


/*
 * Randomise RNG (initial state must be unique across a large cluster of machines)
 */
void
                randomise(void)
{
#ifdef WIN32
    DWORD i, hash=2166136261, seed1, seed2, seed3, seed4;
    TCHAR  infoBuf[256];
    DWORD  bufCharCount = 256;
    SYSTEMTIME time_struct;
 
    if (!GetComputerName(infoBuf, &bufCharCount))
	fail("Cannot query computer name!");
    
    for (i=0; i<bufCharCount; i++)
	hash = (hash * 16777619) ^ infoBuf[i]; /* FNV hash */
    
    GetLocalTime(&time_struct);
    
    seed1 = getpid();
    seed2 = time_struct.wMilliseconds;
    seed3 = hash;
    seed4 = time(NULL);
#else
    unsigned int i, seed1, seed2, seed3, seed4;
    struct timeval tv;

    if (gettimeofday(&tv, NULL))
	fail("randomise: cannot generate random number seeds!");
    
    seed1 = getpid();
    seed2 = tv.tv_usec;
    seed3 = gethostid();
    seed4 = tv.tv_sec;
#endif

    for (i=0; i<32; i++)
	wstate[i] = lcghash(seed1+i) + lcghash(seed2+i) + lcghash(seed3+i) + lcghash(seed4+i);

    widx = 0;

    printf("Random seeds: %u %u %u %u\n", seed1, seed2, seed3, seed4);
}

/* Generate random number 0<=x<1 */
#define ran0()  (WELL1024()*(1.0/4294967296.0))

/* randint(a,b) : return random integer a <= n <= b */
#define randint(low,high) ((low) + (int)(((high)-(low)+1) * ran0()))


/* 
   Apply a transformation matrix to a point:
   Transform       transform; transformation to apply to the point
   Vector          p;         the point to transformation
   Vector          tp;        the returned point after transformation
 */
void
                transform_point(Transform transform, Vector p, Vector tp)
{
    Vector           temp;

    temp[0] = p[0] + transform[0][3];
    temp[1] = p[1] + transform[1][3];
    temp[2] = p[2] + transform[2][3];
    tp[0] = dotprod(transform[0], temp);
    tp[1] = dotprod(transform[1], temp);
    tp[2] = dotprod(transform[2], temp);
}


/* Calculate transformation matrix for coordinate frame defined by 3 points */
void
                calcxf(Vector p1, Vector p2, Vector p3, Transform xf)
{
    int             i;
    Vector          rx, ry, rz, temp;
    float           m;

    vecsub(rz, p2, p1);
    m = 1.0F / sqrtf(dotprod(rz, rz));
    vecscale(rz, rz, m);
    vecsub(temp, p3, p1);
    vecprod(rx, temp, rz);
    m = 1.0F / sqrtf(dotprod(rx, rx));
    vecscale(rx, rx, m);
    vecprod(ry, rz, rx);
    for (i = 0; i < 3; i++)
    {
	xf[0][i] = rx[i];
	xf[1][i] = ry[i];
	xf[2][i] = rz[i];
	xf[3][i] = 0.0F;
	xf[i][3] = -p1[i];
    }
    xf[3][3] = 1.0F;
}


/* Position atom d relative to abc with bond length, bond angle phi, and torsion angle theta */
void setdihedral(Vector a, Vector b, Vector c, Vector d, float blencd, float theta, float phi)
{
    float m;
    Vector AB, bc, n, nbc, d2;

    /* Unit vector bc */
    vecsub(bc, c, b);
    m = 1.0F / sqrtf(dotprod(bc, bc));
    vecscale(bc, bc, m);

    d2[0] = blencd * cosf(theta);
    d2[1] = blencd * cosf(phi) * sinf(theta);
    d2[2] = blencd * sinf(phi) * sinf(theta);
    
    /* Vector AB */
    vecsub(AB, b, a);

    /* Unit normal */
    vecprod(n, AB, bc);
    m = 1.0F / sqrtf(dotprod(n, n));
    vecscale(n, n, m);

    vecprod(nbc, n, bc);
 
    /* NeRF method by Parsons et al. */
    d[0] = bc[0] * d2[0] + nbc[0] * d2[1] + n[0] * d2[2];
    d[1] = bc[1] * d2[0] + nbc[1] * d2[1] + n[1] * d2[2];
    d[2] = bc[2] * d2[0] + nbc[2] * d2[1] + n[2] * d2[2];
    
    vecadd(d, d, c);
}


/* Splice a fragment from src to dest+didx */
void            splice(RESATM *dest, RESATM *src, int didx, int length)
{
    int             i;
    Transform       xf1, xf2;

/*    printf("splice: %d %d %d\n", didx, length, svrflg[didx]); */
    
    /* Calculate transformation for source fragment based on the first residue (C-alpha at the origin) */
    calcxf(src->n, src->ca, src->c, xf2);

    /* Unless splice is at the beginning we need to transform the residues before insertion point */
    if (didx)
    {
	calcxf(dest[didx].n, dest[didx].ca, dest[didx].c, xf1);
	
	for (i = 0; i < didx; i++)
	{
	    transform_point(xf1, dest[i].n, dest[i].n);
	    transform_point(xf1, dest[i].ca, dest[i].ca);
	    transform_point(xf1, dest[i].c, dest[i].c);
	    transform_point(xf1, dest[i].o, dest[i].o);
	    transform_point(xf1, dest[i].cb, dest[i].cb);
	}
    }
    
    if (didx + length < seqlen)
	calcxf(dest[didx + length - 1].n, dest[didx + length - 1].ca, dest[didx + length - 1].c, xf1);

    for (i = 0; i < length; i++)
    {
	if (didx + i >= seqlen)
	{
	    printf("*** Cannot carry out splice (%d %d) !!\n", didx, length);
	    exit(1);
	}

	transform_point(xf2, src[i].n, dest[didx + i].n);
	transform_point(xf2, src[i].ca, dest[didx + i].ca);
	transform_point(xf2, src[i].c, dest[didx + i].c);
	transform_point(xf2, src[i].o, dest[didx + i].o);
	transform_point(xf2, src[i].cb, dest[didx + i].cb);
    }

    /* If the fragment isn't placed at the end then we must transform the following chain segment */ 
    if (didx + length < seqlen)
    {
	calcxf(dest[didx + length - 1].n, dest[didx + length - 1].ca, dest[didx + length - 1].c, xf2);

	for (i = didx + length; i < seqlen; i++)
	{
	    transform_point(xf1, dest[i].n, dest[i].n);
	    transform_point(xf1, dest[i].ca, dest[i].ca);
	    transform_point(xf1, dest[i].c, dest[i].c);
	    transform_point(xf1, dest[i].o, dest[i].o);
	    transform_point(xf1, dest[i].cb, dest[i].cb);
	}

	for (i = 0; i < didx + length; i++)
	{
	    transform_point(xf2, dest[i].n, dest[i].n);
	    transform_point(xf2, dest[i].ca, dest[i].ca);
	    transform_point(xf2, dest[i].c, dest[i].c);
	    transform_point(xf2, dest[i].o, dest[i].o);
	    transform_point(xf2, dest[i].cb, dest[i].cb);
	}
    }
}


/* Splice a fragment from src to dest+didx keeping rest of dest chain fixed */
void            fixed_splice(RESATM * dest, RESATM * src, int didx, int length)
{
    int             i;
    Transform       xf1, xf2;

/*    printf("fixed_splice: %d %d %d\n", didx, length, svrflg[didx]); */

    /* Calculate transform for source fragment based on the first residue (C-alpha at the origin) */
    calcxf(src->n, src->ca, src->c, xf2);

    /* Calculate transform for dest fragment based on the first residue (C-alpha at the origin) */
    calcxf(dest[didx].n, dest[didx].ca, dest[didx].c, xf1);
	
    for (i = 0; i < seqlen; i++)
    {
	transform_point(xf1, dest[i].n, dest[i].n);
	transform_point(xf1, dest[i].ca, dest[i].ca);
	transform_point(xf1, dest[i].c, dest[i].c);
	transform_point(xf1, dest[i].o, dest[i].o);
	transform_point(xf1, dest[i].cb, dest[i].cb);
    }
    
    for (i = 0; i < length; i++)
    {
/*	printf("occ at %d : %d\n", didx+i, svrflg[didx+i]); */
	if (i > 0)
	    assert (svrflg[didx+i] != 0);
	
	if (didx + i >= seqlen)
	{
	    printf("*** Cannot carry out fixed_splice (didx=%d len=%d) !!\n", didx, length);
	    exit(1);
	}

	transform_point(xf2, src[i].n, dest[didx + i].n);
	transform_point(xf2, src[i].ca, dest[didx + i].ca);
	transform_point(xf2, src[i].c, dest[didx + i].c);
	transform_point(xf2, src[i].o, dest[didx + i].o);
	transform_point(xf2, src[i].cb, dest[didx + i].cb);
    }
}


/* Write PDB file */
void
writemainpdb(char *fname, RESATM * chain, int start, int end)
{
    FILE           *ofp;
    int             i, atomn;

    ofp = fopen(fname, "w");
    if (ofp != NULL)
    {
	fprintf(ofp, "HEADER MAIN CHAIN\n");

	for (atomn = i = start; i <= end; i++)
	{
	    fprintf(ofp, "ATOM   %4d %s %s  %4d    %8.3f%8.3f%8.3f  1.00%6.2f\n",
  ++atomn, " N  ", rnames[seq[0][i]], i + 1, chain[i].n[0], chain[i].n[1], chain[i].n[2], 1.0F + svrflg[i]);
	    fprintf(ofp, "ATOM   %4d %s %s  %4d    %8.3f%8.3f%8.3f  1.00%6.2f\n",
  ++atomn, " CA ", rnames[seq[0][i]], i + 1, chain[i].ca[0], chain[i].ca[1], chain[i].ca[2], 1.0F + svrflg[i]);
	    fprintf(ofp, "ATOM   %4d %s %s  %4d    %8.3f%8.3f%8.3f  1.00%6.2f\n",
  ++atomn, " C  ", rnames[seq[0][i]], i + 1, chain[i].c[0], chain[i].c[1], chain[i].c[2], 1.0F + svrflg[i]);
	    fprintf(ofp, "ATOM   %4d %s %s  %4d    %8.3f%8.3f%8.3f  1.00%6.2f\n",
  ++atomn, " O  ", rnames[seq[0][i]], i + 1, chain[i].o[0], chain[i].o[1], chain[i].o[2], 1.0F + svrflg[i]);
	    if (seq[0][i] != GLY)
		fprintf(ofp, "ATOM   %4d %s %s  %4d    %8.3f%8.3f%8.3f  1.00%6.2f\n",
		++atomn, " CB ", rnames[seq[0][i]], i + 1, chain[i].cb[0], chain[i].cb[1], chain[i].cb[2], 1.0F + svrflg[i]);
	}
	fprintf(ofp, "END\n");

	fclose(ofp);
    }
}


/* Combine fragment taking into account modelling mode and template constraints */
void frag_combine(RESATM * dest, RESATM * src, int didx, int length)
{
    int i;

    for (i=1; i<length; i++)
	if (segidx[didx+i] != segidx[didx+i-1])
	    break;

    if (i != length)
	splice(dest, src, didx, length);
    else if (mod_mode == MODELM)
    {
	if (svrflg[didx] > 1)
	    splice(dest, src, didx, length);
	else if (svrflg[didx] == 1) 
	    fixed_splice(dest, src, didx, length);
#if 0
	else if (!svrflg[didx] && svrflg[didx+1] == 1)
	    fixed_splice(dest, src, didx, length);
#endif
	else
	{
	    puts("SVR MAP:");
	    for (i=0; i<seqlen; i++)
		putchar('0' + svrflg[i]);
	    putchar('\n');
	    printf("%d %d %d %d\n", didx, length, svrflg[didx], svrflg[didx+1]);
	    fail("frag_combine: Bad splice position!");
	}
    }
    else
	splice(dest, src, didx, length);
}


/* Check that copied fragment agrees with sec. struc. */
int chksst(int chain, int from, int to, int len)
{
    int i;

    if (mod_mode == MODELM)
    {
	for (i=0; i<len; i++)
	    if (svrflg[to+i] && !(tplt_ss[to + i] & chnlist[chain].sstruc[from + i]))
		break;
    }
    else
    {
	for (i=0; i<len; i++)
	    if (!(tplt_ss[to + i] & chnlist[chain].sstruc[from + i]))
		break;
    }
    
    return i == len;
}


/* Return potential term for at1(aa-type-a) -> at2(aa-type-b) */
float pairpot(const int a, const int b, const int at1, const int at2, const int t, const float d)
{
    int    k;

    if (t > SRTOPOMAX)
    {
	k = (d-4.0F) / RRWIDTH;
	if (k < 0)
	    k = 0;
	else if (k >= NRR)
	    return 0.0F;

	return lr_de[at1][at2][a][b][k];
    }

    if (t >= 1)
    {
	k = INTERVALS * (d - mindist[t-1][at1][at2]) / diffdist[t-1][at1][at2];
	if (k < 0)
	    k = 0;
	if (k < INTERVALS)
	    return sr_de[at1][at2][a][b][k][t-1];
    }

    return 0.0F;
}


/* Return potential term for CB(i) -> CB(j) (from precomputed matrix with linear interpolation) */
float cbcb_pairpot(const int i, const int j, const float d)
{
    int k, t;
    float df;

    t = j - i;

    if (t > SRTOPOMAX)
    {
	df = (d-4.0F) / RRWIDTH;
	k = (int) df;
	df -= k;

	if (k < 0 || (k == 0 && df <= 0.5F))
	    return cbcb_potmat[i][j][0];
	if (k > NRR-1 || (k == NRR-1 && df >= 0.5F))
	    return cbcb_potmat[i][j][NRR-1];

	if (df < 0.5F)
	    return cbcb_potmat[i][j][k-1] * (0.5F - df) + (0.5F + df) * cbcb_potmat[i][j][k];
	else
	    return cbcb_potmat[i][j][k+1] * (df - 0.5F) + (1.5F - df) * cbcb_potmat[i][j][k];
    }

    if (t >= 1)
    {
	df = INTERVALS * (d - mindist[t-1][CBATOM][CBATOM]) / diffdist[t-1][CBATOM][CBATOM];
	k = (int) df;
	df -= k;

	if (k < 0 || (k == 0 && df <= 0.5F))
	    return cbcb_potmat[i][j][0];
	if (k > INTERVALS-1 || (k == INTERVALS-1 && df >= 0.5F))
	    return cbcb_potmat[i][j][INTERVALS-1];

	if (df < 0.5F)
	    return cbcb_potmat[i][j][k-1] * (0.5F - df) + (0.5F + df) * cbcb_potmat[i][j][k];
	else
	    return cbcb_potmat[i][j][k+1] * (df - 0.5F) + (1.5F - df) * cbcb_potmat[i][j][k];
    }

    return 0.0F;
}


/* Build fast lookup matrix for consensus CB->CB and solvation potentials */
void build_emats()
{
    int   i, j, k, n, t;
    float weight[nseqs], sum;

    for (i = 0; i < nseqs; i++)
	weight[i] = 0.0;

    for (i = 0; i < nseqs; i++)
    {
	for (sum=k=0; k<seqlen; k++)
	    sum += rmsdcmat[seq[i][k]][seq[0][k]];
	weight[i] += sum;
	if (verboseflg)
	    printf("Seq %d sim = %f\n", i, sum);
    }
    
    for (sum = i = 0; i < nseqs; i++)
	sum += weight[i];

    for (i = 0; i < nseqs; i++)
	weight[i] = weight[i] / sum;

    if (verboseflg)
	for (i = 0; i < nseqs; i++)
	    printf("Seq %d weight = %f\n", i+1, weight[i]);

    cbcb_potmat = allocmat(seqlen, seqlen, sizeof(float *));

    solv_potmat = allocmat(seqlen, (OOIMAX - OOIMIN) / OOIDIV + 1, sizeof(float));

    for (i = 0; i < seqlen; i++)
    {
	for (n = 0; n < nseqs; n++)
	{
	    if (seq[n][i] < GAP)
		for (k=0; k < (OOIMAX - OOIMIN) / OOIDIV + 1; k++)
		    solv_potmat[i][k] += ooi_de[seq[n][i]][k] * weight[n];
	    else if (seq[n][i] == GAP)
		for (k=0; k < (OOIMAX - OOIMIN) / OOIDIV + 1; k++)
		    solv_potmat[i][k] += ooi_de[ASP][k] * weight[n];
	}

	for (j = i + 1; j < seqlen; j++)
	{
	    t = j - i;

	    if (t > SRTOPOMAX)
	    {
		if ((cbcb_potmat[i][j] = calloc(NRR, sizeof(float))) == NULL)
		    fail("build_cbcbtab: calloc failed!");

		for (n = 0; n < nseqs; n++)
		    if (seq[n][i] < GAP && seq[n][j] < GAP)
		    {
			for (k=0; k<NRR; k++)
			    cbcb_potmat[i][j][k] += lr_de[CBATOM][CBATOM][seq[n][i]][seq[n][j]][k] * weight[n];
		    }
	    }
	    else
	    {
		if ((cbcb_potmat[i][j] = calloc(INTERVALS, sizeof(float))) == NULL)
		    fail("build_cbcbtab: calloc failed!");

		for (n = 0; n < nseqs; n++)
		    if (seq[n][i] < GAP && seq[n][j] < GAP)
			for (k=0; k<INTERVALS; k++)
			    cbcb_potmat[i][j][k] += sr_de[CBATOM][CBATOM][seq[n][i]][seq[n][j]][k][t-1] * weight[n];
	    }
	}
    }
}


/* Calculate dfire potential between two atoms */
float dfirepot(int atmtyp1, int atmtyp2, float distsq)
{
    int k;

    if (distsq >= SQR(14.5F))
	return 0.0F;

    if (distsq < SQR(2.0F))
	k = 0;
    else if (distsq < SQR(8.0F))
	k = 1 + (sqrtf(distsq)-2.0F) / 0.5F;
    else
	k = 13 + (sqrtf(distsq)-8.0F);

    return dfire[atmtyp1][atmtyp2][k];
}


/* Read in rotamers file */
void readrots(void)
{
    int aa, i, atomn = 0;
    float x, y, z, perc;
    char atmnam[10], buf[160], rotalibname[160];
    FILE *rfp;
    extern char    *getenv();

    if (getenv("FDATA_DIR"))
    {
	strcpy(rotalibname, getenv("FDATA_DIR"));
	if (rotalibname[strlen(rotalibname)-1] != '/')
	    strcat(rotalibname, "/");
    }
    else
	rotalibname[0] = '\0';

    strcat(rotalibname, "rotalib.dat");
    
    if (!(rfp = fopen(rotalibname, "r")))
	fail("Cannot read rotalib.dat!");

    while (!feof(rfp))
    {
	if (!fgets(buf, 160, rfp))
	    break;

	if (buf[0] == '*')
	{
	    atomn = 0;
	    nrots[aa]++;
	    continue;
	}

	if (buf[1] == ' ')
	{
	    aa = aanum(buf[0]);
	    if (aa >= 20)
		fail("readrots: Unknown aa type!");
	    if (sscanf(buf+2, "%f", &perc) != 1)
		fail("readrots: Bad file format!");
	    rotafreq[aa][nrots[aa]] = perc;
	    continue;
	}
	
	if (sscanf(buf+3, "%f%f%f", &x, &y, &z) != 3)
	    break;

	if (aa != GLY && atomn < nscatoms[aa])
	{
	    rotamers[aa][nrots[aa]][atomn].pos[0] = x;
	    rotamers[aa][nrots[aa]][atomn].pos[1] = y;
	    rotamers[aa][nrots[aa]][atomn].pos[2] = z;
	    sscanf(buf, "%s", atmnam);
	    if ((rotamers[aa][nrots[aa]][atomn].atmnam = strdup(atmnam)) == NULL)
		fail("readrots: Out of memory!");
	    switch(atmnam[0])
	    {
	    case 'C':
		rotamers[aa][nrots[aa]][atomn].type = 0;
		break;
	    case 'N':
		rotamers[aa][nrots[aa]][atomn].type = 1;
		break;
	    case 'O':
		rotamers[aa][nrots[aa]][atomn].type = 2;
		break;
	    case 'S':
		rotamers[aa][nrots[aa]][atomn].type = 3;
		break;
	    default:
		puts(buf);
		fail("readrots: Unknown atom type in rotamer file!");
	    }
	}
	atomn++;
    }

    if (verboseflg)
	for (i=0; i<20; i++)
	    printf("%c %d rotamers\n", rescodes[i], nrots[i]);
}


/* Construct side chains */
void
buildschn(RESATM *chn, int first, int last)
{
    int             aa, i, j;
    Vector           rx, ry, rz, temp, vsum;
    float          m;

    for (i=first; i<=last; i++)
    {
	aa = seq[0][i];

	if (aa >= 20)
	    fail("buildschn: cannot build side chain for unknown residue type!");
	
	vecsub(rz, chn[i].cb, chn[i].ca);
	m = 1.0F / sqrtf(dotprod(rz, rz));
	vecscale(rz, rz, m);
	vecsub(temp, chn[i].n, chn[i].ca);
	vecprod(rx, temp, rz);
	m = 1.0F / sqrtf(dotprod(rx, rx));
	vecscale(rx, rx, m);
	vecprod(ry, rz, rx);

	if (chn[i].sc == NULL && nscatoms[aa])
	{
	    if ((chn[i].sc = malloc(sizeof(SIDECATM) * nscatoms[aa])) == NULL)
		fail("buildschn: cannot allocate memory for side chains!");
	    chn[i].nscats = nscatoms[aa];
	}

	/* Calculate side chain centroid position */
	if (nscatoms[aa])
	{
	    veczero(chn[i].sc_cg);
	    for (j=0; j<nscatoms[aa]; j++)
	    {
		chn[i].sc[j] = rotamers[aa][chn[i].rotnum][j];
		vecscale(vsum, rx, chn[i].sc[j].pos[0]);
		vecscale(temp, ry, chn[i].sc[j].pos[1]);
		vecadd(vsum, vsum, temp);
		vecscale(temp, rz, chn[i].sc[j].pos[2]);
		vecadd(vsum, vsum, temp);
		vecadd(chn[i].sc[j].pos, vsum, chn[i].ca);
		vecadd(chn[i].sc_cg, chn[i].sc_cg, chn[i].sc[j].pos);
	    }
	
	    vecscale(chn[i].sc_cg, chn[i].sc_cg, 1.0F/nscatoms[aa]);
	}
	else
	    veccopy(chn[i].sc_cg, chn[i].ca);
    }
}


/* Main chain hydrogen bond evaluation */
float           calchb(float (**dsqmat)[NPAIRS], RESATM * chn)
{
    int             i, j, ndonor[MAXSEQLEN], nacceptor[MAXSEQLEN], nhbs=0, bestj;
    float           hb_energy = 0.0F, costheta, e, ebest;
    float           dnosq;
    Vector          c_nvec, ca_nvec, n_hvec;

    /* Build backbone amide hydrogens (Pauling & Corey rule) */
    for (i = 0; i < seqlen; i++)
    {
	ndonor[i] = nacceptor[i] = 0;

	/* Generate main chain amide H position bisecting C-N-CA (mean of C->N and CA->N directions) */
	if (i)
	{
	    vecsub(c_nvec, chn[i].n, chn[i-1].c);
	    vecscale(c_nvec, c_nvec, 1.0F / dist(chn[i-1].c, chn[i].n));
	    vecsub(ca_nvec, chn[i].n, chn[i].ca);
	    vecscale(ca_nvec, ca_nvec, 1.0F / dist(chn[i].n, chn[i].ca));
	    
	    vecadd(n_hvec, c_nvec, ca_nvec);

	    vecscale(n_hvec, n_hvec, 1.0F / sqrtf(dotprod(n_hvec, n_hvec)));

	    vecadd(chn[i].h, chn[i].n, n_hvec);
	}
    }
    
    for (i = 0; i < seqlen; i++)
	if (seq[0][i] != PRO)
	{
	    /* Find best acceptor for this donor */
	    bestj = -1;
	    ebest = VBIG;
	    for (j = 0; j < seqlen; j++)
		if (abs(j-i) > 2 && dsqmat[i][j][CA_CA] < 81.0F && nacceptor[j] < 2)
		{
		    dnosq = distsq(chn[i].n, chn[j].o);

		    if (dnosq > 25.0F)
			continue;

		    if (i)
		    {
			/* Dot product N->H . O->H (N.B. N-H vector is already unit length) */
			costheta = ((chn[i].h[0] - chn[i].n[0]) * (chn[i].h[0] - chn[j].o[0]) + (chn[i].h[1] - chn[i].n[1]) * (chn[i].h[1] - chn[j].o[1]) + (chn[i].h[2] - chn[i].n[2]) * (chn[i].h[2] - chn[j].o[2])) / dist(chn[i].h, chn[j].o);
			
			/* DREIDING hydrogen bond potential with cos^2 weighting */
			e = 5.5F * (5.0F * powf(8.41F / dnosq, 6.0F) - 6.0F * powf(8.41F / dnosq, 5.0F)) * SQR(costheta);
		    }
		    else
		    {
			/* Non-directional potential for first N-H */
			e = 5.5F * (5.0F * powf(8.41F / dnosq, 6.0F) - 6.0F * powf(8.41F / dnosq, 5.0F));
		    }

		    if (e < ebest)
		    {
			ebest = e;
			bestj = j;
		    }
		}
		    
	    if (ebest < -0.5F)
	    {
		ndonor[i]++;
		nacceptor[bestj]++;
		nhbs++;
		
		if (mod_mode == REFINEM || (abs(bestj-i) > 9 && (tplt_ss[i] & STRAND) && (tplt_ss[bestj] & STRAND)))
		{
		    if (abs(bestj-i) <= vtrange)
		    {
			if (hbfname[0])
			    hb_energy += ebest * hbmat[i][bestj];
			else
			    hb_energy += ebest;
		    }
		    
		    econtrib[i] += ebest * hbwt;
		    econtrib[bestj] += ebest * hbwt;
		}
	    }
	}

    /* printf("NHBs = %d\n", nhbs); */
    
    return hb_energy;
}


/* Rapid RMSD calculation - Ref. Douglas L. Theobald (2005) Acta Crystallographica A 61(4):478-480. */
static void
CalcQuarticCoeffs(RESATM *struc1, RESATM *struc2, const int len, double *coeff)
{
    double         Sxx, Sxy, Sxz, Syx, Syy, Syz, Szx, Szy, Szz;
    double         Szz2, Syy2, Sxx2, Sxy2, Syz2, Sxz2, Syx2, Szy2, Szx2,
                   SyzSzymSyySzz2, Sxx2Syy2Szz2Syz2Szy2, Sxy2Sxz2Syx2Szx2,
                   SxzpSzx, SyzpSzy, SxypSyx, SyzmSzy,
                   SxzmSzx, SxymSyx, SxxpSyy, SxxmSyy;
    double         x1, x2, y1, y2, z1, z2;
    int            i;

    Sxx = Sxy = Sxz = Syx = Syy = Syz = Szx = Szy = Szz = 0.0;

    for (i = 0; i < len; i++)
	if (!svrflg[i])
	{
	    x1 = struc1[i].ca[0];
	    y1 = struc1[i].ca[1];
	    z1 = struc1[i].ca[2];
	    x2 = struc2[i].ca[0];
	    y2 = struc2[i].ca[1];
	    z2 = struc2[i].ca[2];
	    	    
	    Sxx += (x1 * x2);
	    Sxy += (x1 * y2);
	    Sxz += (x1 * z2);
	    
	    Syx += (y1 * x2);
	    Syy += (y1 * y2);
	    Syz += (y1 * z2);
	    
	    Szx += (z1 * x2);
	    Szy += (z1 * y2);
	    Szz += (z1 * z2);  
	}
    
    Sxx2 = Sxx * Sxx;
    Syy2 = Syy * Syy;
    Szz2 = Szz * Szz;

    Sxy2 = Sxy * Sxy;
    Syz2 = Syz * Syz;
    Sxz2 = Sxz * Sxz;

    Syx2 = Syx * Syx;
    Szy2 = Szy * Szy;
    Szx2 = Szx * Szx;

    SyzSzymSyySzz2 = 2.0*(Syz*Szy - Syy*Szz);
    Sxx2Syy2Szz2Syz2Szy2 = Syy2 + Szz2 - Sxx2 + Syz2 + Szy2;

    coeff[2] = -2.0 * (Sxx2 + Syy2 + Szz2 + Sxy2 + Syx2 + Sxz2 + Szx2 + Syz2 + Szy2);
    coeff[1] = 8.0 * (Sxx*Syz*Szy + Syy*Szx*Sxz + Szz*Sxy*Syx - Sxx*Syy*Szz - Syz*Szx*Sxy - Szy*Syx*Sxz);

    SxzpSzx = Sxz+Szx;
    SyzpSzy = Syz+Szy;
    SxypSyx = Sxy+Syx;
    SyzmSzy = Syz-Szy;
    SxzmSzx = Sxz-Szx;
    SxymSyx = Sxy-Syx;
    SxxpSyy = Sxx+Syy;
    SxxmSyy = Sxx-Syy;
    Sxy2Sxz2Syx2Szx2 = Sxy2 + Sxz2 - Syx2 - Szx2;

    coeff[0] = Sxy2Sxz2Syx2Szx2 * Sxy2Sxz2Syx2Szx2
             + (Sxx2Syy2Szz2Syz2Szy2 + SyzSzymSyySzz2) * (Sxx2Syy2Szz2Syz2Szy2 - SyzSzymSyySzz2)
             + (-(SxzpSzx)*(SyzmSzy)+(SxymSyx)*(SxxmSyy-Szz)) * (-(SxzmSzx)*(SyzpSzy)+(SxymSyx)*(SxxmSyy+Szz))
             + (-(SxzpSzx)*(SyzpSzy)-(SxypSyx)*(SxxpSyy-Szz)) * (-(SxzmSzx)*(SyzmSzy)-(SxypSyx)*(SxxpSyy+Szz))
             + (+(SxypSyx)*(SyzpSzy)+(SxzpSzx)*(SxxmSyy+Szz)) * (-(SxymSyx)*(SyzmSzy)+(SxzpSzx)*(SxxpSyy+Szz))
             + (+(SxypSyx)*(SyzmSzy)+(SxzmSzx)*(SxxmSyy-Szz)) * (-(SxymSyx)*(SyzpSzy)+(SxzmSzx)*(SxxpSyy-Szz));
}


/* Evaluates the Newton-Raphson correction for the Horn quartic */
static double
eval_horn_NR_corrxn(const double *c, const double x)
{
    double x2 = x*x;
    double b = (x2 + c[2])*x;
    double a = b + c[1];

    return((a*x + c[0])/(2.0*x2*x + b + a));
}


/* Evaluates the Horn quartic for coefficients c and given x. */
double
eval_horn_quart(const double *c, const double x)
{
    return(((x*x + c[2])*x + c[1])*x + c[0]);
}


/* Evaluates the derivative of the Horn quartic for
   coefficients c and given x. */
double
eval_horn_quart_deriv(const double *c, const double x)
{
    return(2.0*(2.0*x*x + c[2])*x + c[1]);
}


/* Newton-Raphson root finding */
static double
QCProot(double *coeff, double guess, const double delta)
{
    int             i;
    double          oldg;

    for (i = 0; i < 50; ++i)
    {
        oldg = guess;
        guess -= eval_horn_NR_corrxn(coeff, guess);
    
        if (fabs(guess - oldg) < fabs(delta*guess))
            return(guess);
    }

    fprintf(stderr,
            "\n\n ERROR21: Newton-Raphson root-finding in \'QCProot()\' did not converge \n");

    exit(EXIT_FAILURE);
}


/* Calculate the inner product of some coordinates.
   This is the same as the squared radius of gyration without normalization
   for the number of atoms. */
static double
CoordsInnerProd(RESATM *struc, const int len)
{
    int             i;
    double          sum;

    sum = 0.0;
    for (i = 0; i < len; ++i)
	if (!svrflg[i])
	    sum += dotprod(struc[i].ca, struc[i].ca);

    return(sum);
}


static void
CenterCoords(RESATM *struc, const int len)
{
    int             i, n;
    double          xav, yav, zav;

    xav = yav = zav = 0.0;
    for (n = i = 0; i < len; i++)
	if (!svrflg[i])
	{
	    xav += struc[i].ca[0];
	    yav += struc[i].ca[1];
	    zav += struc[i].ca[2];
	    n++;
	}

    xav /= n;
    yav /= n;
    zav /= n;

    for (i = 0; i < len; i++)
    {
	struc[i].n[0] -= xav;
	struc[i].n[1] -= yav;
	struc[i].n[2] -= zav;
	struc[i].ca[0] -= xav;
	struc[i].ca[1] -= yav;
	struc[i].ca[2] -= zav;
	struc[i].c[0] -= xav;
	struc[i].c[1] -= yav;
	struc[i].c[2] -= zav;
	struc[i].o[0] -= xav;
	struc[i].o[1] -= yav;
	struc[i].o[2] -= zav;
	struc[i].cb[0] -= xav;
	struc[i].cb[1] -= yav;
	struc[i].cb[2] -= zav;
    }
}

/* returns the sum of the squared deviations (sumdev^2) between the two structures
   rmsd = sqrt(sumdev^2 / atom_num)
*/
double
QuatCharPoly(RESATM *struc1, RESATM *struc2, const int len, double *coeff)
{
    double         innerprod;
    double         lambdamax;

    innerprod = CoordsInnerProd(struc1, len) + CoordsInnerProd(struc2, len);
    CalcQuarticCoeffs(struc1, struc2, len, coeff);
    lambdamax = QCProot(coeff, 0.5 * innerprod, 1e-6);

    return innerprod - (2.0 * lambdamax);
}


/* Return RMSD between two chains */
double calcrmsd(RESATM *struc1, RESATM *struc2, int len)
{
    double sumdev2, coeff[3];

    sumdev2 = QuatCharPoly(struc1, struc2, len, coeff);

    if (sumdev2 < -1e4)
	fail("calcrmsd: QuatCharPoly rounding error!");

    sumdev2 = MAX(0.0, sumdev2);

    return sqrt(sumdev2 / len);
}


double calcsqdev(RESATM *struc1, RESATM *struc2, int len)
{
    double sumdev2, coeff[3];

    sumdev2 = QuatCharPoly(struc1, struc2, len, coeff);

    if (sumdev2 < -1e4)
	fail("calcrmsd: QuatCharPoly rounding error!");

    return MAX(0.0, sumdev2);
}


/* Return total CPU time from first call in seconds */
#ifdef WIN32
float           cputime()
{
    return clock() / (float) CLK_TCK;
}
#else
float           cputime()
{
    float           tot;
    static float    ftime = -1.0, tickspersec = -1.0;
    struct tms      t;

    if (tickspersec < 0.0F)
	tickspersec = sysconf(_SC_CLK_TCK);

    (void) times(&t);

    tot = (float) (t.tms_utime + t.tms_stime) / tickspersec;

    if (ftime < 0.0F)
	ftime = tot;

    return tot - ftime;
}
#endif

/* Return total walltime from first call in seconds */
int walltime()
{
    static time_t t0;
    
    if (!t0)
	t0 = time(NULL);
    
    return (int)(time(NULL) - t0);
}

double sr_sum, sr_sumsq, lr_sum, lr_sumsq, comp_sum, comp_sumsq;
double steric_sum, steric_sumsq, solv_sum, solv_sumsq, ds_sum, ds_sumsq, srr_sum, lrr_sum, srr_sumsq, lrr_sumsq;
double hb_sum, hb_sumsq, targ_sum, targ_sumsq;

float           last_sr, last_lr, last_compact, last_steric, last_hbond, last_solv, last_ds, last_srr, last_lrr, last_rmsd;

#define evsteric(a, b, dsqthresh) ((dsq = SQR(a[0] - b[0]) + SQR(a[1] - b[1]) + SQR(a[2] - b[2])) < dsqthresh ? dsqthresh - dsq : 0.0F)

/* Compute steric energy contribution */
float calc_steric(float (**dsqmat)[NPAIRS], RESATM * chn, int residx, int ster_mode)
{
    int             i, j, k, l, t, t2, atmtyp1, atmtyp2;
    float           steric;
    float           dsq, radsum;
    float          *at1, *at2;

    steric = 0.0F;

    for (i = 0; i < seqlen; i++)
    {
	for (j = i+1; j < MIN(seqlen, i+vtrange+1); j++)
	{
	    if (residx >= 0 && i != residx && j != residx)
		continue;

	    if (mod_mode == MODELM && !svrflg[i] && !svrflg[j])
		continue;
	    
	    t = j - i;
	    
	    if (!ster_mode && t > 1)
	    {
		radsum = (seq[0][i] != GLY && seq[0][j] != GLY) ? 3.8F : 3.5F;
		if (dsqmat[i][j][CA_CA] < SQR(radsum))
		    steric += dsqmat[i][j][CA_CA] - SQR(radsum);
		steric += evsteric(chn[i].cb, chn[j].cb, SQR(cbmin[seq[0][i]][seq[0][j]]));
	    }
	    if (ster_mode == 1 && t > 1)
	    {
		radsum = (seq[0][i] != GLY && seq[0][j] != GLY) ? 3.8F : 3.5F;
		if (dsqmat[i][j][CA_CA] < SQR(radsum))
		    steric += dsqmat[i][j][CA_CA] - SQR(radsum);
		steric += evsteric(chn[i].sc_cg, chn[j].sc_cg, SQR(cgmin[seq[0][i]][seq[0][j]]));
	    }
	    else if (ster_mode == 2)
	    {
		t2 = MIN(2, t);

		/* SC-SC */
		if (dsqmat[i][j][CB_CB] < 196.0F)
		    for (k = 0; k < nscatoms[seq[0][i]]; k++)
		    {
			atmtyp1 = resatofs[seq[0][i]]+k+4;
			
			for (l = 0; l < nscatoms[seq[0][j]]; l++)
			{
			    atmtyp2 = resatofs[seq[0][j]]+l+4;
			    steric += evsteric(chn[i].sc[k].pos, chn[j].sc[l].pos, atomdistsq[atmtyp1][atmtyp2][t2]);
			}
		    }
		
		if (dsqmat[i][j][CA_CA] < 169.0F)
		{
		    /* SC-MC */
		    for (k = 0; k < nscatoms[seq[0][i]]; k++)
		    {
			atmtyp1 = resatofs[seq[0][i]]+k+4;
			
			atmtyp2 = resatofs[seq[0][j]];

			steric += evsteric(chn[i].sc[k].pos, chn[j].n, atomdistsq[atmtyp1][atmtyp2][t2]);
			
			atmtyp2++;

			steric += evsteric(chn[i].sc[k].pos, chn[j].ca, atomdistsq[atmtyp1][atmtyp2][t2]);
			
			atmtyp2++;

			steric += evsteric(chn[i].sc[k].pos, chn[j].c, atomdistsq[atmtyp1][atmtyp2][t2]);
			
			atmtyp2++;

			steric += evsteric(chn[i].sc[k].pos, chn[j].o, atomdistsq[atmtyp1][atmtyp2][t2]);
		    }
		    
		    /* MC-SC */
		    for (l = 0; l < nscatoms[seq[0][j]]; l++)
		    {
			atmtyp2 = resatofs[seq[0][j]]+l+4;
			
			atmtyp1 = resatofs[seq[0][i]];
			
			steric += evsteric(chn[i].n, chn[j].sc[l].pos, atomdistsq[atmtyp1][atmtyp2][t2]);

			atmtyp1++;

			steric += evsteric(chn[i].ca, chn[j].sc[l].pos, atomdistsq[atmtyp1][atmtyp2][t2]);

			atmtyp1++;

			steric += evsteric(chn[i].c, chn[j].sc[l].pos, atomdistsq[atmtyp1][atmtyp2][t2]);

			atmtyp1++;

			steric += evsteric(chn[i].o, chn[j].sc[l].pos, atomdistsq[atmtyp1][atmtyp2][t2]);
		    }
		    
		    /* MC-MC */
		    if (dsqmat[i][j][CA_CA] < 81.0F)
			for (k = 0; k < 4; k++)
			{
			    atmtyp1 = resatofs[seq[0][i]] + k;
			    
			    switch (k)
			    {
			    case 0:
				at1 = chn[i].n;
				break;
				
			    case 1:
				at1 = chn[i].ca;
				break;
				
			    case 2:
				at1 = chn[i].c;
				break;
				
			    case 3:
				at1 = chn[i].o;
				break;
			    }
			    
			    for (l = 0; l < 4; l++)
			    {
				atmtyp2 = resatofs[seq[0][j]] + l;
				
				switch (l)
				{
				case 0:
				    at2 = chn[j].n;
				    break;
				    
				case 1:
				    at2 = chn[j].ca;
				    break;
				    
				case 2:
				    at2 = chn[j].c;
				    break;
				    
				case 3:
				    at2 = chn[j].o;
				    break;
				}
				
				steric += evsteric(at1, at2, atomdistsq[atmtyp1][atmtyp2][t2]);
			    }
			}
		}
	    }
	    else if (ster_mode == 3)
	    {
		/* DFIRE potential */
		
		/* SC-SC */
		if (dsqmat[i][j][CB_CB] < 625.0F)
		    for (k = 0; k < nscatoms[seq[0][i]]; k++)
			for (l = 0; l < nscatoms[seq[0][j]]; l++)
			{
			    atmtyp1 = resatofs[seq[0][i]]+k+4;
			    atmtyp2 = resatofs[seq[0][j]]+l+4;
			    dsq = distsq(chn[i].sc[k].pos, chn[j].sc[l].pos);
			    steric += dfirepot(atmtyp1, atmtyp2, dsq);
#if 0
			    if (dfirepot(atmtyp1, atmtyp2, dsq) > 0.01F)
				printf("%d %d %f %f\n", i, j, sqrtf(dsq), dfirepot(atmtyp1, atmtyp2, d));
#endif
			}
		
		if (dsqmat[i][j][CA_CA] < 576.0F && t > 1)
		{
		    /* SC-MC */
		    for (k = 0; k < nscatoms[seq[0][i]]; k++)
		    {
			atmtyp1 = resatofs[seq[0][i]]+k+4;
			
			atmtyp2 = resatofs[seq[0][j]];
			dsq = distsq(chn[i].sc[k].pos, chn[j].n);
			steric += dfirepot(atmtyp1, atmtyp2, dsq);
			
			atmtyp2++;
			dsq = distsq(chn[i].sc[k].pos, chn[j].ca);
			steric += dfirepot(atmtyp1, atmtyp2, dsq);
			
			atmtyp2++;
			dsq = distsq(chn[i].sc[k].pos, chn[j].c);
			steric += dfirepot(atmtyp1, atmtyp2, dsq);
			
			atmtyp2++;
			dsq = distsq(chn[i].sc[k].pos, chn[j].o);
			steric += dfirepot(atmtyp1, atmtyp2, dsq);
		    }
		    
		    /* MC-SC */
		    for (l = 0; l < nscatoms[seq[0][j]]; l++)
		    {
			atmtyp2 = resatofs[seq[0][j]]+l+4;
			
			atmtyp1 = resatofs[seq[0][i]];
			dsq = distsq(chn[j].sc[l].pos, chn[i].n);
			steric += dfirepot(atmtyp1, atmtyp2, dsq);
			
			atmtyp1++;
			dsq = distsq(chn[j].sc[l].pos, chn[i].ca);
			steric += dfirepot(atmtyp1, atmtyp2, dsq);
			
			atmtyp1++;
			dsq = distsq(chn[j].sc[l].pos, chn[i].c);
			steric += dfirepot(atmtyp1, atmtyp2, dsq);
			
			atmtyp1++;
			dsq = distsq(chn[j].sc[l].pos, chn[i].o);
			steric += dfirepot(atmtyp1, atmtyp2, dsq);
		    }
		    
		    /* MC-MC */
		    if (dsqmat[i][j][CA_CA] < 400.0F)
			for (k = 0; k < 4; k++)
			{
			    atmtyp1 = resatofs[seq[0][i]] + k;
			    
			    switch (k)
			    {
			    case 0:
				at1 = chn[i].n;
				break;
				
			    case 1:
				at1 = chn[i].ca;
				break;
				
			    case 2:
				at1 = chn[i].c;
				break;
				
			    case 3:
				at1 = chn[i].o;
				break;
			    }
			    
			    for (l = 0; l < 4; l++)
			    {
				atmtyp2 = resatofs[seq[0][j]] + l;
				
				switch (l)
				{
				case 0:
				    at2 = chn[j].n;
				    break;
				    
				case 1:
				    at2 = chn[j].ca;
				    break;
				    
				case 2:
				    at2 = chn[j].c;
				    break;
				    
				case 3:
				    at2 = chn[j].o;
				    break;
				}
				
				dsq = distsq(at1, at2);
				steric += dfirepot(atmtyp1, atmtyp2, dsq);
			    }
			}
		}
	    }
	}
    }

    return steric;
}


float e_model(float (**dsqmat)[NPAIRS], RESATM * chn, int debug)
{
    int             i, j, k, t, n_sr, n_lr;
    float           p, sr, lr, compact, steric, hbond, solv, ds, srrcon, lrrcon, lclosure;
    float           rmsd, dviol;
    short           dspart[MAXSEQLEN];

    sr = lr = compact = steric = solv = ds = srrcon = lrrcon = rmsd = hbond = lclosure = 0.0F;
    n_sr = n_lr = 0;

    if (dswt > 0.0F)
    {
	if (n_ds)
	    for (i = 0; i < n_ds; i++)
	    {
		ds += SQR(sqrtf(dsqmat[ds_from[i]][ds_to[i]][CA_CA]) - 5.7F);
		ds += SQR(sqrtf(dsqmat[ds_from[i]][ds_to[i]][CB_CB]) - 4.0F);
	    }
	else
	{
	    for (i = 0; i < seqlen; i++)
		dspart[i] = 0;
	    for (i = 0; i < seqlen; i++)
		if (seq[0][i] == CYS && !dspart[i])
		{
		    for (j = i+4; j < seqlen; j++)
			if (seq[0][j] == CYS && !dspart[j])
			{
			    p = SQR(sqrtf(dsqmat[i][j][CA_CA]) - 5.7F);
			    p += SQR(sqrtf(dsqmat[i][j][CB_CB]) - 4.0F);
			    if (p < 1.0F)
			    {
				ds -= 1.0F;
				dspart[i] = dspart[j] = 1;
				break;
			    }
			}
		}
	}
    }
    
    for (i = 0; i < seqlen; i++)
	econtrib[i] = 0.0F;

    if (hbwt > 0.0F)
	hbond = calchb(dsqmat, chn);

/*    printf("Ehb = %f\n", hbond); */

    for (i = 0; i < seqlen; i++)
    {
	for (j = i+1; j < seqlen; j++)
	{
	    if (mod_mode == MODELM && !svrflg[i] && !svrflg[j])
		continue;
	    
	    t = j - i;
	    
	    if (t > vtrange)
		continue;

	    if (t > 4 && dsqmat[i][j][CA_CA] <= 36.0F)
	        compact -= 1.0F;

	    if (mod_mode == FOLDM && t == 1 && (dsqmat[i][j][CA_CA] < 7.29F || dsqmat[i][j][CA_CA] > 16.0F))
	    {
		printf("Invalid virtual bond distance (%d-%d : %5.2f)!\n", i, j, sqrtf(dsqmat[i][j][CA_CA]));
		/* exit(1); */
	    }

	    if (t == 1)
		lclosure += SQR(sqrtf(dsqmat[i][j][CA_CA]) - 3.8F);

	    if (rrcbflag > 1)
	    {
		if (dcomat[i][j] > 0.0F && conmat[i][j] > 8.0F && conmat[i][j] < 14.0F)
		{
		    if (abs(j-i) <= 23)
			srrcon += SQR(sqrtf(dsqmat[i][j][CA_CA]) - conmat[i][j]) / dcomat[i][j];
		    else
			lrrcon += SQR(sqrtf(dsqmat[i][j][CA_CA]) - conmat[i][j]) / dcomat[i][j];
		}
	    }
	    else
	    {
		if (conmat && conmat[i][j] != 0.0F && j-i > 5)
		{
		    /* short-range contacts */
		    if (srrwt > 0.0F && abs(j-i) < 23)
		    {
			if ((tplt_ss[i] & STRAND) && (tplt_ss[j] & STRAND))
			{
			    dviol = sqrtf(dsqmat[i][j][CA_CA]) - 6.5F;
			}
			else
			    dviol = sqrtf(dsqmat[i][j][CB_CB]) - condmax[seq[0][i]][seq[0][j]];
			
			if (dviol <= 0.0F)
			    srrcon += conmat[i][j];
			else
			    srrcon += conmat[i][j] * expf(-dviol*dviol);
//			srrcon += (conmat[i][j] * expf(-dviol*dviol)) - (conmat[i][j] * dviol / (dviol + 8.0));
		    }
		    else if (srrwt < 0.0F && abs(j-i) < 23)
		    { 
			if (rrcbflag)
			    dviol = sqrtf(dsqmat[i][j][CB_CB]) - dcomat[i][j];
			else
			    dviol = sqrtf(dsqmat[i][j][CA_CA]) - dcomat[i][j];
			
			if (dviol < 0.0F)
			    srrcon += conmat[i][j];
		    }
                    /* long-range contacts */		
		    if (lrrwt > 0.0F && j-i >= 23)
		    { 
			if ((tplt_ss[i] & STRAND) && (tplt_ss[j] & STRAND))
			{
			    dviol = sqrtf(dsqmat[i][j][CA_CA]) - 6.5F;
			}
			else
			    dviol = sqrtf(dsqmat[i][j][CB_CB]) - condmax[seq[0][i]][seq[0][j]];
			
			if (dviol <= 0.0F)
			    lrrcon += conmat[i][j];
			else
			    lrrcon += conmat[i][j] * expf(-dviol*dviol);
//			lrrcon += (conmat[i][j] * expf(-dviol*dviol)) - (conmat[i][j] * dviol / (dviol + 8.0));
		    }
		    else if (lrrwt < 0.0F && j-i >= 23)
		    { 
			if (rrcbflag)
			    dviol = sqrtf(dsqmat[i][j][CB_CB]) - dcomat[i][j];
			else
			    dviol = sqrtf(dsqmat[i][j][CA_CA]) - dcomat[i][j];
			
			if (dviol < 0.0F)
			    lrrcon += conmat[i][j];
		    }
		}
	    }
	    
	    if (dsqmat[i][j][CB_CB] < SQR(SRDCUTOFF))
	    {
		p = cbcb_pairpot(i, j, sqrtf(dsqmat[i][j][CB_CB]));

		if (t <= SRTOPOMAX)
		{
		    sr -= p;
		    n_sr++;
		}
		else if (dsqmat[i][j][CB_CB] < SQR(LRDCUTOFF))
		{
		    lr -= p;
		    n_lr++;
		}
	    }
	}
	
	p = solv_potmat[i][MAX(0, MIN(OOIMAX, tcbooi[i]) - OOIMIN) / OOIDIV];
	solv -= p;

	/* Relate coordinate error to burial of residue */
	econtrib[i] = -tcbooi[i];
    }

    steric = calc_steric(dsqmat, chn, -1, global_ster_mode);

//    lclosure = sqrtf(lclosure/(seqlen-1));

    if (cb_targ != NULL)
    {
#if 0
	for (rmsd=k=i=0; i<seqlen; i++)
	    if (!svrflg[i])
		for (j=i+1; j<seqlen; j++)
		    if (!svrflg[j] && cb_targ[i][j] > 0.0)
		    {
			rmsd += SQR(cb_targ[i][j] - sqrtf(dsqmat[i][j][CB_CB]));
			k++;
		    }
	rmsd = sqrt(rmsd / k);
#else
	rmsd = calcrmsd(targchn, chn, seqlen);
#endif
	last_rmsd = rmsd;
    }
    
    for (ecmin = VBIG, k = 0; k < seqlen; k++)
	if (econtrib[k] < ecmin)
	    ecmin = econtrib[k];

    for (ecmax = -VBIG, k = 0; k < seqlen; k++)
	if (econtrib[k] > ecmax)
	    ecmax = econtrib[k];

    if (debug > 0 && verboseflg)
    {
	printf("\nEsr = %g\nElr = %g\n", sr, lr);
	printf("Esolv = %g\nEcompact = %g\nEhb = %g\nEsteric = %g\nDisulf = %g\nSRcontacts = %g\nLRcontacts = %g\nLoop = %g\n", solv, compact, hbond, steric, ds, srrcon, lrrcon, lclosure);
	if (cb_targ != NULL)
	{
	    printf("RMSD = %f\n", rmsd);
#if 0
	    if (fabsf(rmsd - lastrmsd) > 0.5F)
	    {
		if (rmsd > 0.1F)
		    printf("# %g %g %g %g %g %g %g %g %g %g %g %g %g\n", rmsd, sr, srwt, lr, lrwt, solv, solvwt, compact, compactwt, hbond, hbwt, steric, stericwt);
		lastrmsd = rmsd;
	    }
#endif
	}

	if (targwt <= 0.0F)
	    printf("TOT = %f\n", srwt * sr + lrwt * lr + solvwt * solv + compactwt * compact + hbwt * hbond + stericwt * steric + dswt * ds + fabs(srrwt) * srrcon + fabs(lrrwt) * lrrcon + fabsf(targwt) * (lclosure + rmsd));
	else
	    printf("TOT = %f\n", srwt * sr + lrwt * lr + solvwt * solv + compactwt * compact + hbwt * hbond + stericwt * steric + dswt * ds + fabs(srrwt) * srrcon + fabs(lrrwt) * lrrcon + fabsf(targwt) * (lclosure + ((rmsd  > 1.0F && rmsd < REFRANGE) ? 0.0F : fabsf(rmsd - 0.5F * (REFRANGE + 1.0F)))));

	fflush(stdout);
    }

    sr_sum += sr;
    sr_sumsq += SQR(sr);
    lr_sum += lr;
    lr_sumsq += SQR(lr);
    comp_sum += compact;
    comp_sumsq += SQR((double)compact);
    steric_sum += steric;
    steric_sumsq += SQR(steric);
    solv_sum += solv;
    solv_sumsq += SQR(solv);
    hb_sum += hbond;
    hb_sumsq += SQR(hbond);
    ds_sum += ds;
    ds_sumsq += SQR(ds);
    srr_sum += srrcon;
    srr_sumsq += SQR(srrcon);
	lrr_sum += lrrcon;
    lrr_sumsq += SQR(lrrcon);
    targ_sum += rmsd;
    targ_sumsq += SQR(rmsd);

    last_sr = sr;
    last_lr = lr;
    last_compact = compact;
    last_steric = steric;
    last_hbond = hbond;
    last_solv = solv;
    last_srr = srrcon;
	last_lrr = lrrcon;
    last_ds = ds;

    if (targwt <= 0.0)
	return srwt * sr + lrwt * lr + solvwt * solv + compactwt * compact + hbwt * hbond + stericwt * steric + dswt * ds + fabs(srrwt) * srrcon + fabs(lrrwt) * lrrcon + fabsf(targwt) * (lclosure + rmsd);
    else
	return srwt * sr + lrwt * lr + solvwt * solv + compactwt * compact + hbwt * hbond + stericwt * steric + dswt * ds + fabs(srrwt) * srrcon + fabs(lrrwt) * lrrcon + fabsf(targwt) * (lclosure + ((rmsd > 1.0F && rmsd < REFRANGE) ? 0.0F : fabsf(rmsd - 0.5F * (REFRANGE + 1.0F))));
}


/* Byte-swap for little-endian systems */
void byteswap(void *f, register int n)
{
    register unsigned int *vp;

    vp = (unsigned int *)f;

    while (n--)
    {
        *vp = ((*vp >> 24) & 0xFF) | ((*vp >> 8) & 0x0000FF00) | ((*vp << 8) & 0x00FF0000) | ((*vp << 24) & 0xFF000000);
	vp++;
    }
}

struct atomhashentry
{
    unsigned sernum;
    int residx;
} atomhash[16777216];

#define qhash(x) (2654435761U * (unsigned)(x))

/* Generate distance matrix and Ooi numbers - return TRUE immediately if chain overlaps itself */
int calcdist(const RESATM * chain, const int check_clash)
{
    int             i, j, t, tt;
    unsigned hash;
    float dsq;
    static const int tpc[6] = {4, 5, 6, 3, 7, 2}; /* List of most likely clashing seq. separations in order */
    static unsigned sernum;
    
    if (check_clash)
    {
	sernum++;
	
	/* 1st stage clash check by geometric hashing */
	for (i = 0; i < seqlen; i++)
	{
	    hash = (qhash((int)chain[i].ca[0] / 4) + qhash((int)chain[i].ca[1] / 4) + qhash((int)chain[i].ca[2] / 4)) & 16777215U;
	    if (atomhash[hash].sernum == sernum)
	    {
		if (atomhash[hash].residx < i-1)
		{
		    dsq = distsq(chain[i].ca, chain[atomhash[hash].residx].ca);
		    if (dsq < 14.44F)
			if ((seq[0][i] != GLY && seq[0][atomhash[hash].residx] != GLY) || dsq < 12.25F)
			    return TRUE;
		}
	    }
	    else
	    {
		atomhash[hash].sernum = sernum;
		atomhash[hash].residx = i;
	    }
	}

	/* 2nd stage clash check in sequence separation order of clash likelihood */
	for (tt = 2; tt < seqlen; tt++)
	{
	    if (tt >= 8)
		t = tt;
	    else
		t = tpc[tt-2];

	    for (i = 0; i < seqlen-t; i++)
	    {
		j = i+t;

		dsqmat[i][j][CA_CA] = distsq(chain[i].ca, chain[j].ca);
		if (dsqmat[i][j][CA_CA] < 14.44F)
		    if ((seq[0][i] != GLY && seq[0][j] != GLY) || dsqmat[i][j][CA_CA] < 12.25F)
			return TRUE;
	    }
	}

	for (i = 0; i < seqlen-1; i++)
	    dsqmat[i][i+1][CA_CA] = dsqmat[i+1][i][CA_CA] = distsq(chain[i].ca, chain[i+1].ca);
    }
    else
	for (i = 0; i < seqlen; i++)
	    for (j = i + 1; j < seqlen; j++)
		dsqmat[i][j][CA_CA] = dsqmat[j][i][CA_CA] = distsq(chain[i].ca, chain[j].ca);

    for (i = 0; i < seqlen; i++)
	tcbooi[i] = 0;

    for (i = 0; i < seqlen; i++)
	for (j = i + 1; j < seqlen; j++)
	{	    
	    dsqmat[i][j][CB_CB] = dsqmat[j][i][CB_CB] = distsq(chain[i].cb, chain[j].cb);

	    /* Compute CB-Ooi numbers */
	    if (dsqmat[i][j][CB_CB] <= SQR(OOICUTOFF))
	    {
		tcbooi[i]++;
		tcbooi[j]++;
	    }
	}
    
#if 0
    for (i = 0; i < seqlen; i++)
	for (j = 0; j < seqlen; j++)
	{
	    dsqmat[i][j][CB_N] = dsqmat[j][i][N_CB] = (SQR(chain[i].cb[0] - chain[j].n[0]) + SQR(chain[i].cb[1] - chain[j].n[1]) + SQR(chain[i].cb[2] - chain[j].n[2]));
	    dsqmat[i][j][CB_O] = dsqmat[j][i][O_CB] = (SQR(chain[i].cb[0] - chain[j].o[0]) + SQR(chain[i].cb[1] - chain[j].o[1]) + SQR(chain[i].cb[2] - chain[j].o[2]));
	}
#endif

    return FALSE;
}


/* Count severe CA-CA clashes */
int clashcount(const RESATM * chain, int length)
{
    int             i, j, nclash=0;
    float dsq;

    for (i = 0; i < length; i++)
	for (j = i + 2; j < length; j++)
	{
	    dsq = distsq(chain[i].ca, chain[j].ca);
 	    if (dsq < 12.25F)
		nclash++;
	}

    return nclash;
}


/* Copy chain segment with checking for side chain atom memory allocation */
void chaincpy(RESATM *dest, RESATM *src, int n)
{
    SIDECATM        *oldsc;

    while (n--)
    {
	if (dest->sc != NULL && src->sc == NULL)
	    continue;
	oldsc = dest->sc;
	memcpy(dest, src, sizeof(RESATM));
	dest->sc = oldsc;

	if (src->sc != NULL)
	{
	    if (dest->sc == NULL)
		if ((dest->sc = malloc(sizeof(SIDECATM) * dest->nscats)) == NULL)
		    fail("chaincpy: cannot allocate memory for side chains!");
	
	    memcpy(dest->sc, src->sc, dest->nscats * sizeof(SIDECATM));
	}
	
	dest++;
	src++;
    }
}


void
                readchn(void)
{
    FILE           *dfp, *tfp;
    char            brkid[10], brkidlst[100000], tdbname[512], buf[512], tdbloc[160], torname[160], seq[MAXSEQLEN], sstruc[MAXSEQLEN];
    int             i, n, helix = 0, strand = 0, tot = 0, relacc[MAXSEQLEN];
    RESATM          chain[MAXSEQLEN];
    extern char    *getenv();

    if (getenv("FDATA_DIR"))
    {
	strcpy(torname, getenv("FDATA_DIR"));
	if (torname[strlen(torname)-1] != '/')
	    strcat(torname, "/");
    }
    else
	torname[0] = '\0';

    strcat(torname, "tor.lst");
    if (!(tfp = fopen(torname, "r")))
	fail("readchn: cannot open tor.lst!");

    brkidlst[0] = '\0';
    while (!feof(tfp))
    {
	if (fscanf(tfp, "%s", brkid) != 1)
	    break;
	if (brkid[0] == '#')
	    continue;
	strcat(brkidlst, brkid);
    }

    fclose(tfp);

    /* Read coords from TDB file */
    if (getenv("FTDB_DIR"))
    {
	strcpy(tdbname, getenv("FTDB_DIR"));
	if (tdbname[strlen(tdbname)-1] != '/')
	    strcat(tdbname, "/");
    }
    else
	tdbname[0] = '\0';
    
    strcat(tdbname, "tdb.dat");
    dfp = fopen(tdbname, "r");
    if (!dfp)
	fail("readchn: cannot open tdb.dat file!");

    if (!fgets(buf, 160, dfp))
	fail("readchn: bad tdb.dat file!");

    while (!feof(dfp))
    {
	if (buf[0] != '#')
	    fail("readchn: bad tdb.dat file!");

	if (sscanf(buf, "%*s%s", brkid) != 1)
	    break;

	n = 0;
	
	while (!feof(dfp))
	{
	    if (!fgets(buf, 512, dfp))
		break;

	    if (buf[0] == '#')
		break;

	    if (brkid[0] != '-' && !strstr(brkidlst, brkid))
		continue;
	    
	    if (n >= MAXSEQLEN)
		fail("readchn: MAXSEQLEN exceeded!");
	    
	    if (buf[7] == 'H' || buf[7] == 'G')
	    {
		helix++;
		sstruc[n] = HELIX;
	    }
	    else if (toupper(buf[7]) == 'E' || toupper(buf[7]) == 'A' || toupper(buf[7]) == 'P')
	    {
		strand++;
		sstruc[n] = STRAND;
	    }
	    else
		sstruc[n] = COIL;
	    tot++;
	    sscanf(buf + 9, "%d", &relacc[n]);
	    seq[n] = aanum(buf[5]);
	    if (sscanf(buf + 13, "%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f", &chain[n].phi, &chain[n].psi, &chain[n].omega, chain[n].n, chain[n].n + 1, chain[n].n + 2, chain[n].ca, chain[n].ca + 1, chain[n].ca + 2, chain[n].c, chain[n].c + 1, chain[n].c + 2, chain[n].o, chain[n].o + 1, chain[n].o + 2, chain[n].cb, chain[n].cb + 1, chain[n].cb + 2) != 18)
		fail("readchn: bad tdb file!");
	    chain[n].sc = NULL;
	    chain[n].rotnum = 0;
	    
	    n++;
	}

#if 0
/* Check for mistakes in chain geometry - use only when creating new tdb.dat file */
	if (n < 20)
	    continue;

	for (i = 0; i < n; i++)
	    if (seq[i] >= 20)
		break;
	
	if (i != n)
	{
	    printf("Unknown residue in chain!\n");
	    continue;
	}
	
	for (i = 0; i < n - 1; i++)
	{
	    bondl = dist(chain[i].c, chain[i + 1].n);
	    if (bondl > 1.45F)
		break;

#if 0
	    printf("! ");
	    for (j=1; j<8; j++)
	    {
		if (i+j < n)
		{
		    printf(" %7.3f", dist(chain[i].ca, chain[i + j].ca));
		}
		else
		    printf(" %7.3f", 0.0);
	    }
	    putchar('\n');
#endif
	}

	if (i != n - 1)
	{
	    printf("Bondlen = %f!\n", bondl);
	    continue;
	}

	/* Check backbone geometry */
	for (i = 0; i < n; i++)
	{
	    bondl = dist(chain[i].n, chain[i].c);
	    if (fabsf(bondl - 2.46F) / 0.0475F > 3.0F)
	    {
		printf("Bad N..C dist at %d = %f!\n", i+1, bondl);
		break;
	    }
	    bondl = dist(chain[i].n, chain[i].ca);
	    if (fabsf(bondl - 1.46F) / 0.0162F > 3.0F)
	    {
		printf("Bad N-CA dist at %d = %f!\n", i+1, bondl);
		break;
	    }
	    bondl = dist(chain[i].ca, chain[i].c);
	    if (fabsf(bondl - 1.52F) / 0.0159F > 3.0F)
	    {
		printf("Bad CA-C dist at %d = %f!\n", i+1, bondl);
		break;
	    }
	    bondl = dist(chain[i].c, chain[i].o);
	    if (fabsf(bondl - 1.23F) / 0.0103F > 4.0F)
	    {
		printf("Bad C=O dist at %d = %f!\n", i+1, bondl);
		break;
	    }
	}

	if (i != n)
	    continue;
	
	for (i = 0; i < n - 1; i++)
	{
	    bondl = dist(chain[i].ca, chain[i + 1].ca);

	    /* Check virtual bond distance for cis and trans peptides */
	    if (bondl < 2.8F || bondl > 3.9F || (bondl >= 3.2F && bondl < 3.7F))
	    {
		printf("Bad virtual bond (%f A)!\n", bondl);
		break;
	    }
	}

	if (i != n-1)
	    continue;

	if (clashcount(chain, n) > 0)
	{
	    printf("Bad CA-CA bumps in %s chain!\n", brkid);
	    continue;
	}

	printf("OKCHAIN: %s\n", brkid);
#endif

	chnlist[nchn].length = n;
	if (!(chnlist[nchn].chain = malloc(n * sizeof(RESATM))) ||
	    !(chnlist[nchn].relacc = malloc(n * sizeof(short))) ||
	    !(chnlist[nchn].seq = malloc(n)) ||
	    !(chnlist[nchn].sstruc = malloc(n)))
	        fail("readchn: malloc failed!");

	for (i = 0; i < n; i++)
	{
	    chaincpy(chnlist[nchn].chain + i, chain + i, 1);
	    chnlist[nchn].relacc[i] = relacc[i];
	    chnlist[nchn].seq[i] = seq[i];
	    chnlist[nchn].sstruc[i] = sstruc[i];
	}

	nchn++;
	chnres += n;
    }
    
    fclose(dfp);

    if (verboseflg)
    {
	printf("%d chains read (%d residues)\n", nchn, chnres);
	printf("%5.1f %% helix\n", 100.0 * helix / tot);
	printf("%5.1f %% strand\n", 100.0 * strand / tot);
    }
}


/* Generate fully extended chain */
void
                fullyextended(RESATM *chn)
{
    int             i;
    float           sx, sy;
    Vector          cca, nca, xx, yy;

    /* Generate fully extended chain to start without clashes */
    chn[0].n[0] = chn[0].n[1] = chn[0].n[2] = 0.0;
    chn[0].ca[0] = 1.45;
    chn[0].ca[1] = chn[0].ca[2] = 0.0;
    chn[0].c[0] = chn[0].ca[0] + cos(1.2) * 1.52;
    chn[0].c[1] = sin(1.2) * 1.52;
    chn[0].c[2] = 0.0;
    
    setdihedral(chn[0].n, chn[0].ca, chn[0].c, chn[0].o, 1.23, 58.9 * PI / 180.0, 0.0);
    
    for (i=1; i<seqlen; i++)
    {
	setdihedral(chn[i-1].n, chn[i-1].ca, chn[i-1].c, chn[i].n, 1.33, 64.4 * PI / 180.0, PI);
	setdihedral(chn[i-1].ca, chn[i-1].c, chn[i].n, chn[i].ca, 1.45, 58.1 * PI / 180.0, PI);
	setdihedral(chn[i-1].ca, chn[i-1].c, chn[i].n, chn[i].h, 1.0, 58.1 * PI / 180.0, 0.0);
	setdihedral(chn[i-1].c, chn[i].n, chn[i].ca, chn[i].c, 1.52, 68.8 * PI / 180.0, PI);
	setdihedral(chn[i].n, chn[i].ca, chn[i].c, chn[i].o, 1.23, 58.9 * PI / 180.0, 0.0);
    }
    
    /* Generate CB positions */
    for (i=0; i<seqlen; i++)
    {
	vecsub(nca, chn[i].ca, chn[i].n);
	vecsub(cca, chn[i].ca, chn[i].c);
	vecadd(xx, nca, cca);
	vecprod(yy, nca, cca);
	sx = CACBDIST * cosf(TETH_ANG) / sqrtf(dotprod(xx, xx));
	sy = CACBDIST * sinf(TETH_ANG) / sqrtf(dotprod(yy, yy));
	vecscale(xx, xx, sx);
	vecscale(yy, yy, sy);
	vecadd(chn[i].cb, chn[i].ca, xx);
	vecadd(chn[i].cb, chn[i].cb, yy);
    }
}

/* Estimate coordinate error from residue burial */
#define errorest(x) (0.806452F * MAX(x, -OOIMAX) + 25.0F)

/* Write PDB file */
void
                writepdb(RESATM * chain, float e)
{
    FILE           *ofp;
    int             i,j,atomn;
    float rmsdiff;
    char outpdbn2[160];
    static  int modeln;

    ofp = fopen(outpdbn, "w");
    if (ofp != NULL)
    {
	if (mqapcmd == NULL)
	    fprintf(ofp, "HEADER Potential terms: %f %f %f %f %f %f %f %f %f %f\n", e, last_sr, last_lr, last_solv, last_compact, last_hbond, last_steric, last_ds, last_srr, last_lrr);
	else
	    fprintf(ofp, "HEADER %f\n", e);

	for (atomn = i = 0; i < seqlen; i++)
	{
	    fprintf(ofp, "ATOM   %4d %s %s  %4d    %8.3f%8.3f%8.3f  1.00%6.2f\n",
		    ++atomn, " N  ", rnames[seq[0][i]], i + 1, chain[i].n[0], chain[i].n[1], chain[i].n[2], errorest(econtrib[i]));
#if 0
	    if (hbwt > 0.0 && i > 0 && seq[0][i] != PRO)
		fprintf(ofp, "ATOM   %4d %s %s  %4d    %8.3f%8.3f%8.3f  1.00%6.2f\n",
			++atomn, " H  ", rnames[seq[0][i]], i + 1, chain[i].h[0], chain[i].h[1], chain[i].h[2], errorest(econtrib[i]));
#endif
	    fprintf(ofp, "ATOM   %4d %s %s  %4d    %8.3f%8.3f%8.3f  1.00%6.2f\n",
		    ++atomn, " CA ", rnames[seq[0][i]], i + 1, chain[i].ca[0], chain[i].ca[1], chain[i].ca[2], errorest(econtrib[i]));
	    fprintf(ofp, "ATOM   %4d %s %s  %4d    %8.3f%8.3f%8.3f  1.00%6.2f\n",
		    ++atomn, " C  ", rnames[seq[0][i]], i + 1, chain[i].c[0], chain[i].c[1], chain[i].c[2], errorest(econtrib[i]));
	    fprintf(ofp, "ATOM   %4d %s %s  %4d    %8.3f%8.3f%8.3f  1.00%6.2f\n",
		    ++atomn, " O  ", rnames[seq[0][i]], i + 1, chain[i].o[0], chain[i].o[1], chain[i].o[2], errorest(econtrib[i]));
	    if (seq[0][i] != GLY)
		fprintf(ofp, "ATOM   %4d %s %s  %4d    %8.3f%8.3f%8.3f  1.00%6.2f\n",
			++atomn, " CB ", rnames[seq[0][i]], i + 1, chain[i].cb[0], chain[i].cb[1], chain[i].cb[2], errorest(econtrib[i]));
	    if (chain[i].sc && seq[0][i] != GLY && seq[0][i] != ALA)
	    {
		for (j=1; j<nscatoms[seq[0][i]]; j++)
		    fprintf(ofp, "ATOM   %4d  %-3.3s %s  %4d    %8.3f%8.3f%8.3f  1.00%6.2f\n",
			    ++atomn, chain[i].sc[j].atmnam, rnames[seq[0][i]], i + 1, chain[i].sc[j].pos[0], chain[i].sc[j].pos[1], chain[i].sc[j].pos[2], errorest(econtrib[i]));
	    }
	}
	fprintf(ofp, "END\n");

	fclose(ofp);
    }

/* Save RMSD/energy data */

#if 0
    for (rmsdiff = i = 0; i < seqlen; i++)
	for (j = i + 1; j < seqlen; j++)
	    rmsdiff += SQR(sqrtf(dsqmat[i][j][CB_CB]) - cb_last[i][j]);

    if (sqrt(2.0*rmsdiff/seqlen/(seqlen-1)) < 0.5)
	return;

    for (i = 0; i < seqlen; i++)
	for (j = i + 1; j < seqlen; j++)
	    cb_last[i][j] = sqrtf(dsqmat[i][j][CB_CB]);

    ofp = fopen("energy.log", "a");
    if (!ofp)
	fail("Cannot write to energy log!");
    fprintf(ofp, "%f %f %f %f %f %f %f %f %f %f %f %f ", last_sr, srwt, last_lr, lrwt, last_compact, compactwt, last_steric, stericwt, last_hbond, hbwt, last_solv, solvwt);
    fclose(ofp);
    system("justrmsd target.pdb fold.pdb | tail -1 >> energy.log");
#endif

#ifdef SAVETRAJ
#if 0
    if (mod_mode == FOLDM)
    {
	for (rmsdiff = i = 0; i < seqlen; i++)
	    for (j = i + 1; j < seqlen; j++)
		rmsdiff += SQR(sqrtf(dsqmat[i][j][CB_CB]) - cb_targ[i][j]);
	
	if (sqrt(2.0*rmsdiff/seqlen/(seqlen-1)) > 4.0)
	    return;
    }
#endif
    for (rmsdiff = i = 0; i < seqlen; i++)
	for (j = i + 1; j < seqlen; j++)
	    rmsdiff += SQR(sqrtf(dsqmat[i][j][CB_CB]) - cb_last[i][j]);

    if (modeln && sqrt(2.0*rmsdiff/seqlen/(seqlen-1)) < 0.5)
	return;

    for (i = 0; i < seqlen; i++)
	for (j = i + 1; j < seqlen; j++)
	    cb_last[i][j] = sqrtf(dsqmat[i][j][CB_CB]);

    sprintf(outpdbn2, "TR%04d_%s", ++modeln, outpdbn);
    printf("Saving %s ...\n", outpdbn2);
    ofp = fopen(outpdbn2, "w");
    if (ofp != NULL)
    {
	if (mqapcmd == NULL)
	    fprintf(ofp, "HEADER Potential terms: %f %f %f %f %f %f %f %f\n", e, last_sr, last_lr, last_compact, last_steric, last_hbond, last_solv, last_ds);
	else
	    fprintf(ofp, "HEADER %f\n", e);
	for (atomn = i = 0; i < seqlen; i++)
	{
	    fprintf(ofp, "ATOM   %4d %s %s  %4d    %8.3f%8.3f%8.3f  1.00%6.2f\n",
		    ++atomn, " N  ", rnames[seq[0][i]], i + 1, chain[i].n[0], chain[i].n[1], chain[i].n[2], errorest(econtrib[i]));
#if 0
	    if (hbwt > 0.0)
		fprintf(ofp, "ATOM   %4d %s %s  %4d    %8.3f%8.3f%8.3f  1.00%6.2f\n",
			++atomn, " H  ", rnames[seq[0][i]], i + 1, chain[i].h[0], chain[i].h[1], chain[i].h[2], errorest(econtrib[i]));
#endif
	    fprintf(ofp, "ATOM   %4d %s %s  %4d    %8.3f%8.3f%8.3f  1.00%6.2f\n",
		    ++atomn, " CA ", rnames[seq[0][i]], i + 1, chain[i].ca[0], chain[i].ca[1], chain[i].ca[2], errorest(econtrib[i]));
	    fprintf(ofp, "ATOM   %4d %s %s  %4d    %8.3f%8.3f%8.3f  1.00%6.2f\n",
		    ++atomn, " C  ", rnames[seq[0][i]], i + 1, chain[i].c[0], chain[i].c[1], chain[i].c[2], errorest(econtrib[i]));
	    fprintf(ofp, "ATOM   %4d %s %s  %4d    %8.3f%8.3f%8.3f  1.00%6.2f\n",
		    ++atomn, " O  ", rnames[seq[0][i]], i + 1, chain[i].o[0], chain[i].o[1], chain[i].o[2], errorest(econtrib[i]));
	    if (seq[0][i] != GLY)
		fprintf(ofp, "ATOM   %4d %s %s  %4d    %8.3f%8.3f%8.3f  1.00%6.2f\n",
			++atomn, " CB ", rnames[seq[0][i]], i + 1, chain[i].cb[0], chain[i].cb[1], chain[i].cb[2], errorest(econtrib[i]));
	    if (seq[0][i] != GLY && seq[0][i] != ALA)
	    {
		for (j=1; j<nscatoms[seq[0][i]]; j++)
		    fprintf(ofp, "ATOM   %4d  %-3.3s %s  %4d    %8.3f%8.3f%8.3f  1.00%6.2f\n",
			    ++atomn, chain[i].sc[j].atmnam, rnames[seq[0][i]], i + 1, chain[i].sc[j].pos[0], chain[i].sc[j].pos[1], chain[i].sc[j].pos[2], errorest(econtrib[i]));
	    }
	}
	fprintf(ofp, "END\n");

	fclose(ofp);
    }
#endif

    return;
}

/* Write PDB file and callout to external MQAP */
float
                runmqap(RESATM * chain)
{
    FILE           *ofp, *mqfp;
    int             i,j,atomn, fd;
    float result;
    char temppdbn[160], cmdstr[160];
    static  int modeln;

    strcpy(temppdbn, "/tmp/fragftmpXXXXXX");
    
    if ((fd = mkstemp(temppdbn)) < 0)
	fail("Cannot create temp file!");
    
    ofp = fdopen(fd, "w");
    if (ofp == NULL)
	fail("Cannot open temp file for writing!");

    for (atomn = i = 0; i < seqlen; i++)
    {
	fprintf(ofp, "ATOM   %4d %s %s  %4d    %8.3f%8.3f%8.3f  1.00%6.2f\n",
		++atomn, " N  ", rnames[seq[0][i]], i + 1, chain[i].n[0], chain[i].n[1], chain[i].n[2], 0.0);
	fprintf(ofp, "ATOM   %4d %s %s  %4d    %8.3f%8.3f%8.3f  1.00%6.2f\n",
		++atomn, " CA ", rnames[seq[0][i]], i + 1, chain[i].ca[0], chain[i].ca[1], chain[i].ca[2], 0.0);
	fprintf(ofp, "ATOM   %4d %s %s  %4d    %8.3f%8.3f%8.3f  1.00%6.2f\n",
		++atomn, " C  ", rnames[seq[0][i]], i + 1, chain[i].c[0], chain[i].c[1], chain[i].c[2], 0.0);
	fprintf(ofp, "ATOM   %4d %s %s  %4d    %8.3f%8.3f%8.3f  1.00%6.2f\n",
		++atomn, " O  ", rnames[seq[0][i]], i + 1, chain[i].o[0], chain[i].o[1], chain[i].o[2], 0.0);
	if (seq[0][i] != GLY)
	    fprintf(ofp, "ATOM   %4d %s %s  %4d    %8.3f%8.3f%8.3f  1.00%6.2f\n",
		    ++atomn, " CB ", rnames[seq[0][i]], i + 1, chain[i].cb[0], chain[i].cb[1], chain[i].cb[2], 0.0);
	if (chain[i].sc && seq[0][i] != GLY && seq[0][i] != ALA)
	{
	    for (j=1; j<nscatoms[seq[0][i]]; j++)
		fprintf(ofp, "ATOM   %4d  %-3.3s %s  %4d    %8.3f%8.3f%8.3f  1.00%6.2f\n",
			++atomn, chain[i].sc[j].atmnam, rnames[seq[0][i]], i + 1, chain[i].sc[j].pos[0], chain[i].sc[j].pos[1], chain[i].sc[j].pos[2], 0.0);
	}
    }
    fprintf(ofp, "END\n");
    
    fclose(ofp);

    sprintf(cmdstr, "%s %s", mqapcmd, temppdbn);

    mqfp = popen(cmdstr, "r");

    if (!mqfp)
	fail("Cannot execute MQAP command!");

    if (fscanf(mqfp, "%f", &result) != 1)
	fail("Cannot interpret output from MQAP!");

    fclose(mqfp);

    unlink(temppdbn);

    return result;
}


/* Write PDB file */
void
                writepdbensemble(RESATM * chain, float e)
{
    FILE           *ofp;
    int             i,j,atomn;

    ofp = fopen("ensemble.pdb", "a");
    if (ofp != NULL)
    {
	fprintf(ofp, "HEADER   *** Etot = %g ***\n", e);
	for (atomn = i = 0; i < seqlen; i++)
	{
	    fprintf(ofp, "ATOM   %4d %s %s  %4d    %8.3f%8.3f%8.3f  1.00%6.2f\n",
		    ++atomn, " N  ", rnames[seq[0][i]], i + 1, chain[i].n[0], chain[i].n[1], chain[i].n[2], errorest(econtrib[i]));
#if 0
	    if (hbwt > 0.0 && i > 0 && seq[0][i] != PRO)
		fprintf(ofp, "ATOM   %4d %s %s  %4d    %8.3f%8.3f%8.3f  1.00%6.2f\n",
			++atomn, " H  ", rnames[seq[0][i]], i + 1, chain[i].h[0], chain[i].h[1], chain[i].h[2], errorest(econtrib[i]));
#endif
	    fprintf(ofp, "ATOM   %4d %s %s  %4d    %8.3f%8.3f%8.3f  1.00%6.2f\n",
		    ++atomn, " CA ", rnames[seq[0][i]], i + 1, chain[i].ca[0], chain[i].ca[1], chain[i].ca[2], errorest(econtrib[i]));
	    fprintf(ofp, "ATOM   %4d %s %s  %4d    %8.3f%8.3f%8.3f  1.00%6.2f\n",
		    ++atomn, " C  ", rnames[seq[0][i]], i + 1, chain[i].c[0], chain[i].c[1], chain[i].c[2], errorest(econtrib[i]));
	    fprintf(ofp, "ATOM   %4d %s %s  %4d    %8.3f%8.3f%8.3f  1.00%6.2f\n",
		    ++atomn, " O  ", rnames[seq[0][i]], i + 1, chain[i].o[0], chain[i].o[1], chain[i].o[2], errorest(econtrib[i]));
	    if (seq[0][i] != GLY)
		fprintf(ofp, "ATOM   %4d %s %s  %4d    %8.3f%8.3f%8.3f  1.00%6.2f\n",
			++atomn, " CB ", rnames[seq[0][i]], i + 1, chain[i].cb[0], chain[i].cb[1], chain[i].cb[2], errorest(econtrib[i]));
	    if (chain[i].sc && seq[0][i] != GLY && seq[0][i] != ALA)
	    {
		for (j=1; j<nscatoms[seq[0][i]]; j++)
		    fprintf(ofp, "ATOM   %4d  %-3.3s %s  %4d    %8.3f%8.3f%8.3f  1.00%6.2f\n",
			    ++atomn, chain[i].sc[j].atmnam, rnames[seq[0][i]], i + 1, chain[i].sc[j].pos[0], chain[i].sc[j].pos[1], chain[i].sc[j].pos[2], errorest(econtrib[i]));
	    }
	}
	fprintf(ofp, "END\n");

	fclose(ofp);
    }
}


/* Evaluate given set of coordinates */
float
                eval(RESATM * t, int reinit)
{
    int             i;
    float           e, mq, rmsd;
    static unsigned int ncalls;

    if (!reinit)
	ncalls++;

    if (cb_targ)
	CenterCoords(t, seqlen);

    buildschn(t, 0, seqlen-1);

    if (mqapcmd != NULL)
    {
	mq = runmqap(t);
	rmsd = calcrmsd(targchn, t, seqlen);
	e = mq + fabsf(targwt) * ((rmsd < REFRANGE) ? 0.0F : fabsf(rmsd - REFRANGE));
    }
    else
	e = e_model(dsqmat, t, FALSE);
    if (e < e_min - 0.01F && !reinit)
    {
	chaincpy(bestchn, t, seqlen);

	e_min = e;

	if (mqapcmd == NULL)
	    e = e_model(dsqmat, t, TRUE);
	else
	    printf("RMSD = %f MQAP = %f  E = %f\n", rmsd, mq, e);

	if (verboseflg)
	{
	    for (i=0; i<seqlen; i++)
		if (bestchn[i].rotnum < 10)
		    printf("%d", bestchn[i].rotnum);
		else
		    putchar('>');
	    puts("\n");

	    printf("%f conformations per CPU second.\n", (float) ncalls / cputime());
	}

	if (fabsf(e - e_min) > 0.0001F * fabsf(e))
	{
	    printf("!!! %f %f\n", e, e_min);
	    printf("%f\n", e_model(dsqmat, t, FALSE));
	    printf("%f\n", e_model(dsqmat, t, FALSE));
	    printf("%f\n", e_model(dsqmat, t, FALSE));
	    printf("%f\n", e_model(dsqmat, t, FALSE));
	    printf("%f\n", e_model(dsqmat, t, TRUE));
	    fail("eval: inconsistent results!");
	}

	if (verboseflg)
	    writepdb(t, e);

	if (rfp)
	{
	    fprintf(rfp, "zap\n");
	    fflush(rfp);
	    fprintf(rfp, "load %s\n", outpdbn);
	    fflush(rfp);
	    fprintf(rfp, "restrict *.CA\n");
	    fflush(rfp);
	    fprintf(rfp, "backbone 20\n");
	    fflush(rfp);
	}
    }
    return e;
}

void sc_cheat(float (**dsqmat)[NPAIRS], RESATM * chn, int sco_ster_mode)
{
    int             i;

    for (i=0; i<seqlen; i++)
	chn[i].rotnum = targchn[i].rotnum;
    buildschn(chn, 0, seqlen-1);
}

/* Iteratively optimize side-chain packing */
void sc_opt(float (**dsqmat)[NPAIRS], RESATM * chn, int sco_ster_mode)
{
    int             i, nr, oldnr, optflg, pass=1;
    float e, e_min;

    buildschn(chn, 0, seqlen-1);

    printf("sc_opt : before = %f\n", e_min = calc_steric(dsqmat, chn, -1, sco_ster_mode));

    do
    {
	for (optflg=i=0; i<seqlen; i++)
	    if (nrots[seq[0][i]] > 1)
	    {
		oldnr = chn[i].rotnum;
		for (nr=0; nr<nrots[seq[0][i]]; nr++)
		{
		    chn[i].rotnum = nr;
		    buildschn(chn, i, i);
		    e = calc_steric(dsqmat, chn, -1, sco_ster_mode);
		    if (e < e_min - 0.001F)
		    {
			e_min = e;
			optflg = TRUE;
			oldnr = nr;
		    }
		}

		chn[i].rotnum = oldnr;
		buildschn(chn, i, i);
	    }
    } while (optflg && pass <= 10);

    printf("sc_opt : after = %f\n", e_min);
}

/* Iteratively optimize side-chain packing */
void testnewsc_opt(float (**dsqmat)[NPAIRS], RESATM * chn, int sco_ster_mode)
{
    int             i, nr, oldnr, optflg;
    float e, e_min;

    buildschn(chn, 0, seqlen-1);

    printf("sc_opt : before = %f\n", calc_steric(dsqmat, chn, -1, sco_ster_mode));

    do
    {
	for (optflg=i=0; i<seqlen; i++)
	    if (nrots[seq[0][i]] > 1)
	    {
		e_min = calc_steric(dsqmat, chn, i, sco_ster_mode) + 0.001F;
		oldnr = chn[i].rotnum;
		e_min += -4.5F * logf(rotafreq[seq[0][i]][oldnr]/(100.0F/nrots[seq[0][i]]));
		printf("e_min = %f\n", e_min);
		for (nr=0; nr<nrots[seq[0][i]]; nr++)
		{
		    chn[i].rotnum = nr;
		    buildschn(chn, i, i);
		    e = calc_steric(dsqmat, chn, i, sco_ster_mode);
		    e += -4.5F * logf(rotafreq[seq[0][i]][nr]/(100.0F/nrots[seq[0][i]]));
/*		    printf("e_min = %f e = %f\n", e_min, e); */
		    if (e < e_min - 0.001F)
		    {
			e_min = e;
			optflg = TRUE;
			oldnr = nr;
			printf("e = %f!!!\n", e);
		    }
		}

		chn[i].rotnum = oldnr;
		buildschn(chn, i, i);
	    }
    } while (optflg);

    printf("sc_opt : after = %f\n", calc_steric(dsqmat, chn, -1, sco_ster_mode));
}


/* Ask Boltzmann if we should accept this energy change! */
int             boltzmann(float de, float t)
{
    if (t > 0.0F && de >= 0.0F)
	return ran0() < exp(-de / t);
    else
	return de < 0.0F;
}


/* Anneal side-chain packing */
void sc_anneal(float (**dsqmat)[NPAIRS], RESATM * chn, int sco_ster_mode)
{
    int             i, nr, oldnr, nfail;
    float e, old_e, t;

    buildschn(chn, 0, seqlen-1);

    printf("sc_anneal : before = %f\n", calc_steric(dsqmat, chn, -1, sco_ster_mode));

    nfail = 0;
    t = 10.0;
    
    do
    {
	/* Set a rotamer to a random state (biased to most common rotamers) */
	i = randint(0, seqlen-1);
	
	if (nrots[seq[0][i]] > 1)
	{
	    old_e = calc_steric(dsqmat, chn, i, sco_ster_mode);
	    oldnr = chn[i].rotnum;
	    old_e += -4.5F * logf(rotafreq[seq[0][i]][oldnr]/(100.0F/nrots[seq[0][i]]));

	    nr = randint(0, nrots[seq[0][i]]-1);
	    chn[i].rotnum = nr;
	    buildschn(chn, i, i);
	    e = calc_steric(dsqmat, chn, i, sco_ster_mode);
	    e += -4.5F * logf(rotafreq[seq[0][i]][nr]/(100.0F/nrots[seq[0][i]]));

	    if (!boltzmann(e - old_e, t))
	    {
		chn[i].rotnum = oldnr;
		buildschn(chn, i, i);
		nfail++;
	    }
	}

	t *= 0.95;
    } while (nfail < 100000 && t > 0);

    printf("sc_anneal : after = %f\n", calc_steric(dsqmat, chn, -1, sco_ster_mode));
}


/* Locate secondary structure elements */
int
                loc_sst(char *tsstruc, int len)
{
    int             i, istart, l, nsst = 0;

    for (i = 0; i < len; i++)
	if (tsstruc[i] != COIL)
	{
	    istart = i;
	    while (i < len && tsstruc[i] == tsstruc[istart])
		i++;
	    l = i - istart;
	    if (l < 3)
		continue;
	    sstlist[nsst].start = istart;
	    sstlist[nsst].length = l;
	    switch (tsstruc[istart])
	    {
	    case HELIX:
		sstlist[nsst].type = HELIX;
		break;
	    case STRAND:
		sstlist[nsst].type = STRAND;
		break;
	    default:
		printf("Code = %d\n", tsstruc[istart]);
		fail("loc_sst: unknown secondary structure code!");
	    }
	    nsst++;
	}

    return nsst;
}



/* Build supersecondary fragment lookup table */
void
                mkfragtb(void)
{
    int             i, j, k, l, m, n, ns, npair, nrms, len, sim, nclose, nfrag,
	nsst, from, to, ncon, ooi[MAXSEQLEN];
    float         **cb_d, pair, seqsim, tot, av_n = 0.0F, minwt, rmsd, rmsdtot,
	minrms, maxrms, av, sd;
    char           *type;
    
    if (verboseflg)
	puts("Building supersecondary motif lists...");
    
    for (i = 0; i < seqlen; i++)
    {
	fraglist[i].totbest = -1000.0F;
	fraglist[i].type = "";
	fraglist[i].frags = malloc(MAXFRAGS * sizeof(struct chnidx));

	if (fraglist[i].frags == NULL)
	    fail("mkfragtb: cannot allocate fraglist!");
    }

    for (i = 0; i < nchn; i++)
    {
	cb_d = allocmat(chnlist[i].length, chnlist[i].length, sizeof(float));
	
	for (j = 0; j < chnlist[i].length; j++)
	    ooi[j] = 0;

	for (j = 0; j < chnlist[i].length; j++)
	    for (k = j + 1; k < chnlist[i].length; k++)
	    {
		cb_d[j][k] = cb_d[k][j] = dist(chnlist[i].chain[j].cb, chnlist[i].chain[k].cb);
		if (j-i > 1 && cb_d[j][k] <= OOICUTOFF)
		{
		    ooi[j]++;
		    ooi[k]++;
		}
	    }
	
	nsst = loc_sst(chnlist[i].sstruc, chnlist[i].length);
/*	printf("nsst = %d\n", nsst); */
	for (n = 0; n < nsst - 1; n++)
	{
	    from = to = -1;

	    if (sstlist[n].type == HELIX && sstlist[n + 1].type == HELIX)
	    {
		ncon = 0;
		for (j = 0; j < sstlist[n].length; j++)
		    for (k = 0; k < sstlist[n + 1].length; k++)
			if (cb_d[sstlist[n].start + j][sstlist[n + 1].start + k] < 8.0F)
			    ncon++;
		if (ncon > 4)
		{
		    if (verboseflg)
			printf("Alpha hairpin at %d\n", sstlist[n].start);
		    type = "ALPHA HAIRPIN";
		    from = sstlist[n].start;
		    to = sstlist[n + 1].start + sstlist[n + 1].length - 1;
		}
		else
		{
		    if (verboseflg)
			printf("Alpha corner at %d\n", sstlist[n].start);
		    type = "ALPHA CORNER";
		    from = sstlist[n].start;
		    to = sstlist[n + 1].start + sstlist[n + 1].length - 1;
		}
	    }
	    else if (sstlist[n].type == STRAND && sstlist[n + 1].type == STRAND)
	    {
		ncon = 0;
		for (j = 0; j < sstlist[n].length; j++)
		    for (k = 0; k < sstlist[n + 1].length; k++)
			if (cb_d[sstlist[n].start + j][sstlist[n + 1].start + k] < 7.0F)
			    ncon++;
		if (ncon > 3)
		{
		    if (verboseflg)
			printf("Beta hairpin at %d\n", sstlist[n].start);
		    type = "BETA HAIRPIN";
		    from = sstlist[n].start;
		    to = sstlist[n + 1].start + sstlist[n + 1].length - 1;
		}
#if 1
		else
		{
		    if (verboseflg)
			printf("Beta corner at %d\n", sstlist[n].start);
		    type = "BETA CORNER";
		    from = sstlist[n].start;
		    to = sstlist[n + 1].start + sstlist[n + 1].length - 1;
		}
#endif
	    }
	    else if (n < nsst - 2 && sstlist[n].type == STRAND && sstlist[n + 1].type == HELIX && sstlist[n + 2].type == STRAND)
	    {
		ncon = 0;
		for (j = 0; j < sstlist[n].length; j++)
		    for (k = 0; k < sstlist[n + 2].length; k++)
			if (cb_d[sstlist[n].start + j][sstlist[n + 2].start + k] < 9.0F)
			    ncon++;
		if (ncon > 2)
		{
		    if (verboseflg)
			printf("BAB unit at %d\n", sstlist[n].start);
		    type = "BAB UNIT";
		    from = sstlist[n].start;
		    to = sstlist[n + 2].start + sstlist[n + 2].length - 1;
		}
		else
		{
		    if (verboseflg)
			printf("Split BAB unit at %d\n", sstlist[n].start);
		    type = "Split BAB UNIT";
		    from = sstlist[n].start;
		    to = sstlist[n + 2].start + sstlist[n + 2].length - 1;
		}
	    }
	    else if (sstlist[n].type == STRAND && sstlist[n + 1].type == HELIX)
	    {
		ncon = 0;
		for (j = 0; j < sstlist[n].length; j++)
		    for (k = 0; k < sstlist[n + 1].length; k++)
			if (cb_d[sstlist[n].start + j][sstlist[n + 1].start + k] < 7.0F)
			    ncon++;
		if (ncon > 2)
		{
		    if (verboseflg)
			printf("BA unit at %d\n", sstlist[n].start);
		    type = "BA UNIT";
		    from = sstlist[n].start;
		    to = sstlist[n + 1].start + sstlist[n + 1].length - 1;
		}
#if 0
		else
		{
		    if (verboseflg)
			printf("BA corner at %d\n", sstlist[n].start);
		    type = "BA CORNER";
		    from = sstlist[n].start;
		    to = sstlist[n + 1].start + sstlist[n + 1].length - 1;
		}
#endif
	    }
	    else if (sstlist[n].type == HELIX && sstlist[n + 1].type == STRAND)
	    {
		ncon = 0;
		for (j = 0; j < sstlist[n].length; j++)
		    for (k = 0; k < sstlist[n + 1].length; k++)
			if (cb_d[sstlist[n].start + j][sstlist[n + 1].start + k] < 7.0F)
			    ncon++;
		if (ncon > 2)
		{
		    if (verboseflg)
			printf("AB unit at %d\n", sstlist[n].start);
		    type = "AB UNIT";
		    from = sstlist[n].start;
		    to = sstlist[n + 1].start + sstlist[n + 1].length - 1;
		}
#if 0
		else
		{
		    if (verboseflg)
			printf("AB corner at %d\n", sstlist[n].start);
		    type = "AB CORNER";
		    from = sstlist[n].start;
		    to = sstlist[n + 1].start + sstlist[n + 1].length - 1;
		}
#endif
	    }

	    if (from >= 0)
	    {
		len = to - from + 1;
		if (len > 1)
		    for (k = 0; k <= seqlen - len; k++)
		    {
			/* Only pick fragments within SVRs */
			if (mod_mode == MODELM)
			{
			    for (l = 1; l < len; l++)
				if (!svrflg[k+l])
				    break;

			    if (l != len)
				continue;
			}
			
			for (sim = l = 0; l < len; l++)
			    if (tplt_ss[k + l] & chnlist[i].sstruc[from + l])
				sim++;

			if (sim < len-1)
			    continue;
			
			for (rmsd = pair = seqsim = 0.0F, npair = nrms = l = 0; l < len; l++)
			{
			    seqsim += dmmmat[chnlist[i].seq[from + l]][seq[0][k + l]];
			    
			    for (m = l + 1; m < len; m++)
			    {
				if (cb_targ != NULL)
				{
				    if (cb_targ[k + l][k + m] > 0.0)
				    {
					rmsd += SQR(cb_d[from + l][from + m] - cb_targ[k + l][k + m]);
					nrms++;
				    }
				}
				
				for (ns = 0; ns < n_ds; ns++)
				    if (ds_from[ns] == k+l && ds_to[ns] == k+m && cb_d[from+l][from+m] > 6.0F)
					pair -= 1000.0F;

				if (cb_d[from + l][from + m] < 15.0F)
				{
				    pair += cbcb_pairpot(k + l, k + m, cb_d[from + l][from + m]);
				    npair++;
				}
			    }
			}
				
/*		    tot = ran0(); */
/*		    tot = pair; */
			
		    tot = 0.1F * seqsim / len + pair / MAX(1, npair);

		    /* Make use of target structure for selecting fragments if modelling */
		    if (fragsel && nrms > 1)
			tot = -rmsd / nrms;

#if 0
		    for (l = 0; l < len; l++)
			printf("%c %f - %c %f\n", rescodes[seq[0][k+l]], bestchn[k+l].cb[0], rescodes[chnlist[i].seq[from+l]], chnlist[i].chain[from+l].cb[0]);

		    printf("%d %d %d %d %d %f\n", k, len, seqlen, from, to, sqrt(-tot));
#endif

#if 0
		    printf("%g %g %g\n", seqsim, pair, tot);
#endif
		    
		    if (tot > fraglist[k].totbest)
		    {
			fraglist[k].totbest = tot;
			fraglist[k].bestlen = len;
			fraglist[k].nchn = i;
			fraglist[k].from = from;
			fraglist[k].type = type;
			fraglist[k].rmsd = rmsd > 0.0F ? sqrt(rmsd / nrms) : 0.0F;
		    }

		    if (fraglist[k].nfrags < MAXFRAGS)
		    {
			fraglist[k].frags[fraglist[k].nfrags].chain = i;
			fraglist[k].frags[fraglist[k].nfrags].pos = from;
			fraglist[k].frags[fraglist[k].nfrags].weight = tot;
			fraglist[k].frags[fraglist[k].nfrags].rmsd = rmsd > 0.0F ? sqrt(rmsd / nrms) : 0.0F;
			fraglist[k].frags[fraglist[k].nfrags].length = len;
			fraglist[k].totwt += tot;
			fraglist[k].totsum += tot;
			fraglist[k].totsumsq += SQR(tot);
			fraglist[k].nfrags++;
			fraglist[k].ntot++;
		    }
		    else
		    {
			minwt = fraglist[k].frags[0].weight;
			j = 0;
			for (l = 1; l < MAXFRAGS; l++)
			    if (fraglist[k].frags[l].weight < minwt)
			    {
				minwt = fraglist[k].frags[l].weight;
				j = l;
			    }
			if (minwt < tot)
			{
			    fraglist[k].frags[j].chain = i;
			    fraglist[k].frags[j].pos = from;
			    fraglist[k].frags[j].length = len;
			    fraglist[k].totwt -= minwt;
			    fraglist[k].frags[j].rmsd = rmsd > 0.0F ? sqrt(rmsd / nrms) : 0.0F;
			    fraglist[k].frags[j].weight = tot;
			    fraglist[k].totwt += tot;
			    fraglist[k].totsum += tot;
			    fraglist[k].totsumsq += SQR(tot);
			    fraglist[k].ntot++;
			}
		    }
		}
	    }
	}
	freemat(cb_d, chnlist[i].length);
    }

    for (rmsdtot = nrms = av_n = nfrag = nclose = i = 0; i < seqlen; i++)
	if (fraglist[i].nfrags)
	{
	    for (len = j = 0; j < fraglist[i].nfrags; j++)
		len = MAX(fraglist[i].frags[j].length, len);
	    av = fraglist[i].totsum / fraglist[i].ntot;
	    sd = sqrt(fraglist[i].totsumsq / fraglist[i].ntot - SQR(av));
	    if (verboseflg)
		printf("%4d %2d %4d %4d %4d %8.3f %8.3f %5.2f %-14s ", i, fraglist[i].bestlen, fraglist[i].nchn, fraglist[i].from, fraglist[i].nfrags, fraglist[i].totbest, (fraglist[i].totbest - av) / sd, fraglist[i].rmsd, fraglist[i].type);
	    minrms = 1000.0F;
	    maxrms = 0.0F;
	    for (j = 0; j < fraglist[i].nfrags; j++)
	    {
		minrms = MIN(minrms, fraglist[i].frags[j].rmsd);
		maxrms = MAX(maxrms, fraglist[i].frags[j].rmsd);
		if (fraglist[i].frags[j].rmsd < 2.0F)
		    nclose++;
		nfrag++;
	    }
	    if (verboseflg)
	    {
		printf(" %4.1f %4.1f", minrms, maxrms);
		putchar('\n');
	    }
	    av_n += fraglist[i].nfrags;
	    rmsdtot += minrms;
	    nrms++;
	}

    if (verboseflg)
    {
	printf("Average Min RMSD = %f\n", rmsdtot / nrms);
	printf("p(RMSD < 2.0) = %f\n", (float) nclose / nfrag);
	printf("Average frags per position = %f\n", av_n / seqlen);
    }
}

/* Build fixed length fragment lookup table */
void
                mkfragtb2(void)
{
    int             i, j, k, l, m, ns, npair, nrms, len, fraglen, from, to, ooi[MAXSEQLEN], nfrag, nclose;
    float         **cb_d, pair, seqsim, tot, av_n = 0.0F, minwt, rmsd, rmsdtot, minrms, maxrms, av, sd;
    char           *type;
    Transform       xf1, xf2;
    Vector point1, point2;

    if (verboseflg)
	puts("Building fixed length fragment lists...");

    for (i = 0; i < seqlen; i++)
    {
	fraglist2[i].totbest = -1000.0F;
	fraglist2[i].type = "";
	fraglist2[i].frags = malloc(MAXFRAGS2 * sizeof(struct chnidx));

	if (fraglist2[i].frags == NULL)
	    fail("mkfragtb: cannot allocate fraglist2!");
    }

    for (i = 0; i < nchn; i++)
    {
	cb_d = allocmat(chnlist[i].length, chnlist[i].length, sizeof(float));
	
	for (j = 0; j < chnlist[i].length; j++)
	    ooi[j] = 0;

	for (j = 0; j < chnlist[i].length; j++)
	    for (k = j + 1; k < chnlist[i].length; k++)
	    {
		cb_d[j][k] = cb_d[k][j] = dist(chnlist[i].chain[j].cb, chnlist[i].chain[k].cb);
		if (j-i > 1 && cb_d[j][k] <= OOICUTOFF)
		{
		    ooi[j]++;
		    ooi[k]++;
		}
	    }

	for (fraglen=MINFRAGS2LEN; fraglen<=MAXFRAGS2LEN; fraglen++)
	    for (from = 0; from <= chnlist[i].length - fraglen; from++)
	    {
		to = from + fraglen - 1;
		len = to - from + 1;
		
		type = malloc(len + 1);
		if (type == NULL)
		    fail("mkfragtb2: malloc failed!");
		
		for (k = 0; k < len; k++)
		    type[k] = sscodes[chnlist[i].sstruc[from + k]];
		type[len] = '\0';
		
		if (from >= 0)
		{
		    if (len > 1)
			for (k = 0; k <= seqlen - len; k++)
			{
			    /* Only pick fragments within SVRs */
			    if (mod_mode == MODELM)
			    {
				for (l = 1; l < len; l++)
				    if (!svrflg[k+l])
					break;
				
				if (l != len)
				    continue;
			    }
#if 1
			    for (l = 0; l < len; l++)
				if (!(tplt_ss[k + l] & chnlist[i].sstruc[from + l]))
				    break;
			    
			    if (l != len)
				continue;
#endif			
			    for (rmsd = pair = seqsim = 0.0F, npair = nrms = l = 0; l < len; l++)
			    {
				seqsim += dmmmat[chnlist[i].seq[from + l]][seq[0][k + l]];
				
				for (m = l + 1; m < len; m++)
				{
				    if (cb_targ != NULL)
				    {
					if (cb_targ[k + l][k + m] > 0.0F)
					{
					    rmsd += SQR(cb_d[from + l][from + m] - cb_targ[k + l][k + m]);
					    nrms++;
					}
				    }
				    
				    for (ns = 0; ns < n_ds; ns++)
					if (ds_from[ns] == k+l && ds_to[ns] == k+m && cb_d[from+l][from+m] > 6.0F)
					    pair -= 1000.0F;
				    
				    if (cb_d[from + l][from + m] < 15.0F)
				    {
					pair += cbcb_pairpot(k + l, k + m, cb_d[from + l][from + m]);
					npair++;
				    }
				}
			    }
			    
			    
/*		    tot = ran0(); */
/*		    tot = pair; */
			    
			    tot = 0.05F * seqsim / len + pair / MAX(1, npair);

			    /* Try to use target structure for selecting fragments if modelling or refining */
			    if (mod_mode == REFINEM && nrms > 1 && !svrflg[k] && !svrflg[k+len-1] && segidx[k] == segidx[k+len-1])
			    {
				    /* Find fragments which will cause small global RMSD changes */
				    calcxf(targchn[k].n, targchn[k].ca, targchn[k].c, xf1);
				    calcxf(chnlist[i].chain[from].n, chnlist[i].chain[from].ca, chnlist[i].chain[from].c, xf2);
				    
				    /* See where fragment will "land" when spliced */
				    transform_point(xf1, targchn[k+len-1].n, point1);
				    transform_point(xf2, chnlist[i].chain[from+len-1].n, point2);
				    tot = -distsq(point1, point2);
				    transform_point(xf1, targchn[k+len-1].ca, point1);
				    transform_point(xf2, chnlist[i].chain[from+len-1].ca, point2);
				    tot -= distsq(point1, point2);
				    transform_point(xf1, targchn[k+len-1].c, point1);
				    transform_point(xf2, chnlist[i].chain[from+len-1].c, point2);
				    tot -= distsq(point1, point2);
			    }
			    else if (fragsel && nrms > 1)
				tot = -rmsd / nrms;

#if 0
			    printf("%g %g\n", pair, tot);
#endif
			    
			    if (tot > fraglist2[k].totbest)
			    {
				fraglist2[k].totbest = tot;
				fraglist2[k].bestlen = len;
				fraglist2[k].nchn = i;
				fraglist2[k].from = from;
				fraglist2[k].type = type;
				fraglist2[k].rmsd = rmsd > 0.0F ? sqrt(rmsd / nrms) : 999.0F;
			    }
			    
			    if (fraglist2[k].nfrags < MAXFRAGS2)
			    {
				fraglist2[k].frags[fraglist2[k].nfrags].chain = i;
				fraglist2[k].frags[fraglist2[k].nfrags].pos = from;
				fraglist2[k].frags[fraglist2[k].nfrags].weight = tot;
				fraglist2[k].frags[fraglist2[k].nfrags].rmsd = rmsd > 0.0F ? sqrt(rmsd / nrms) : 999.0F;
				fraglist2[k].frags[fraglist2[k].nfrags].length = len;
				fraglist2[k].totwt += tot;
				fraglist2[k].totsum += tot;
				fraglist2[k].totsumsq += SQR(tot);
				fraglist2[k].nfrags++;
				fraglist2[k].ntot++;
			    }
			    else
			    {
				minwt = fraglist2[k].frags[0].weight;
				j = 0;
				for (l = 1; l < MAXFRAGS2; l++)
				    if (fraglist2[k].frags[l].weight < minwt)
				    {
					minwt = fraglist2[k].frags[l].weight;
					j = l;
				    }
				if (minwt < tot)
				{
				    fraglist2[k].frags[j].chain = i;
				    fraglist2[k].frags[j].pos = from;
				    fraglist2[k].frags[j].length = len;
				    fraglist2[k].totwt -= minwt;
				    fraglist2[k].frags[j].rmsd = rmsd > 0.0F ? sqrt(rmsd / nrms) : 999.0F;
				    fraglist2[k].frags[j].weight = tot;
				    fraglist2[k].totwt += tot;
				    fraglist2[k].totsum += tot;
				    fraglist2[k].totsumsq += SQR(tot);
				    fraglist2[k].ntot++;
				}
			    }
			}
		}
	    }
	freemat(cb_d, chnlist[i].length);
    }
    
    for (rmsdtot = nrms = av_n = nfrag = nclose = i = 0; i < seqlen; i++)
	if (fraglist2[i].nfrags)
	{
	    for (len = j = 0; j < fraglist2[i].nfrags; j++)
		len = MAX(fraglist2[i].frags[j].length, len);
	    av = fraglist2[i].totsum / fraglist2[i].ntot;
	    sd = sqrt(fraglist2[i].totsumsq / fraglist2[i].ntot - SQR(av));
	    if (verboseflg)
		printf("%4d %2d %4d %4d %4d %8.3f %8.3f %5.2f %-17s ", i, fraglist2[i].bestlen, fraglist2[i].nchn, fraglist2[i].from, fraglist2[i].nfrags, fraglist2[i].totbest, (fraglist2[i].totbest - av) / sd, fraglist2[i].rmsd, fraglist2[i].type);
	    minrms = 1000.0F;
	    maxrms = 0.0F;
	    for (j = 0; j < fraglist2[i].nfrags; j++)
	    {
		minrms = MIN(minrms, fraglist2[i].frags[j].rmsd);
		maxrms = MAX(maxrms, fraglist2[i].frags[j].rmsd);
	    }

	    if (minrms < 1.0F)
		nclose++;
	    nfrag++;

	    if (verboseflg)
	    {
		printf(" %4.1f %4.1f", minrms, maxrms);
		putchar('\n');
	    }
	    av_n += fraglist2[i].nfrags;
	    rmsdtot += minrms;
	    nrms++;
	}
    
    if (verboseflg)
    {
	printf("Average Min RMSD = %f\n", rmsdtot / nrms);
	printf("p(RMSD < 1.0) = %f\n", (float) nclose / nfrag);
	printf("Average frags per position = %f\n", av_n / seqlen);
    }
}


unsigned int movefreq[3] = { 5000, 5000, 5000 };
unsigned int lastmove;
float moveprob[3] = { 0.333333, 0.333333, 0.333333 };

/* Choose large fragment from supersec. or fixed-length list */
void            fragsamp(int *to, int *src, int *from, int *len)
{
    int             fn = -1;

    do
    {
	if (mod_mode == FOLDM)
	    *to = randint(0, seqlen-1);
	else
	    *to = randint(REFSTART, MIN(seqlen, REFEND)-1);

	if (MAXFRAGS > 0 && ran0() < moveprob[0] / (moveprob[0] + moveprob[1]))
	{
	    if (fraglist[*to].nfrags)
	    {
		fn = randint(0, fraglist[*to].nfrags - 1);
		*len = fraglist[*to].frags[fn].length;
		*src = fraglist[*to].frags[fn].chain;
		*from = fraglist[*to].frags[fn].pos;
		lastmove = 0;
	    }
	}
	else
	{
	    if (fraglist2[*to].nfrags)
	    {
		fn = randint(0, fraglist2[*to].nfrags - 1);
		*len = fraglist2[*to].frags[fn].length;
		*src = fraglist2[*to].frags[fn].chain;
		*from = fraglist2[*to].frags[fn].pos;
		lastmove = 1;
	    }
	}
    }
    while (fn < 0);
}


/* Initialize the 'world' */
void
                ga_init(void)
{
    int             i, j, src, from, fn, len, to;
    static short    allocd;
    static RESATM  *save;

    /* Initialize population with fully extended chain */
    fullyextended(curchn);
    
    /* population arrays */
    if (!allocd)
    {
	curpool = calloc(poolsize, sizeof(Schema));
	newpool = calloc(poolsize, sizeof(Schema));
	save = calloc(seqlen, sizeof(RESATM));
    
	samparr = calloc(poolsize, sizeof(int));

	if (!curpool || !newpool || !save || !samparr)
	    fail("ga_init: cannot create population arrays!");
    }

    for (i = 0; i < poolsize; i++)
    {
	if (!allocd)
	{
	    curpool[i].genome = calloc(seqlen, sizeof(RESATM));
	    newpool[i].genome = calloc(seqlen, sizeof(RESATM));

	    if (curpool[i].genome == NULL || newpool[i].genome == NULL)
		fail("ga_init: cannot create schema!");
	}

	if (mod_mode == REFINEM)
	    chaincpy(curpool[i].genome, targchn, seqlen);
	else
	{
	    chaincpy(curpool[i].genome, curchn, seqlen);
	    
#if 0
	    for (;;)
	    {
		src = randint(0, nchn - 1);
		len = MIN(chnlist[src].length, seqlen);
		to = randint(0, seqlen - len);
		from = randint(0, chnlist[src].length - len);
		splice(curpool[i].genome, chnlist[src].chain + from, to, len);
		if (!calcdist(curpool[i].genome, TRUE))
		    break;
		chaincpy(curpool[i].genome, curchn, seqlen);
	    }
#endif
	    
	    /* Initialize population with fragments */
	    for (j = 0; j < seqlen; j++)
	    {
		if (fraglist2[j].nfrags)
		{
		    fn = randint(0, fraglist2[j].nfrags - 1);
		    src = fraglist2[j].frags[fn].chain;
		    from = fraglist2[j].frags[fn].pos;
		    len = fraglist2[j].frags[fn].length;
		    chaincpy(save, curpool[i].genome, seqlen);
		    splice(curpool[i].genome, chnlist[src].chain + from, j, len);
		    if (calcdist(curpool[i].genome, TRUE))
		    {
			chaincpy(curpool[i].genome, save, seqlen);
			len = 1;
		    }
		}
		j += len - 1;
	    }
	}

	(void) calcdist(curpool[i].genome, FALSE);
	curpool[i].perfval = eval(curpool[i].genome, FALSE);
	curpool[i].evalflg = FALSE;
    }

    allocd = TRUE;
}

int             schcmp(const void *sch1, const void *sch2)
{
    if (((Schema *) sch1)->perfval < ((Schema *) sch2)->perfval)
	return -1;
    else if (((Schema *) sch1)->perfval > ((Schema *) sch2)->perfval)
	return 1;

    return 0;
}

/* Sort pool into ascending order of perfval */
void            sortpool(Schema * pool)
{
    qsort((void *) pool, poolsize, sizeof(Schema), schcmp);
}

/* Select new population from old */
void
                gaselect(void)
{
    float           expected;	/* expected number of offspring          */
    float           factor;	/* normalizer for expected value        */
    float           ptr;	/* determines fractional selection       */
    float           sum;	/* control for selection loop           */
    float           fitsum;	/* sum of fitness values */
    int             i, j, k, ncopy, temp, sortkey, rank;

#if 0
    sortpool(curpool);

    /* denominator for ordinal selection probabilities */

    rank = poolsize;

    for (fitsum = i = 0; i < poolsize; i++)
    {
	if (i && curpool[i].perfval > curpool[i-1].perfval)
	    rank--;
	curpool[i].selval = rank;
	fitsum += curpool[i].selval;
	if (curpool[i].perfval == best)
	    besti = i;
    }
#else
    /* denominator for selection probabilities */
    for (fitsum = 0.0F, i = 0; i < poolsize; i++)
    {
	curpool[i].selval = worst - curpool[i].perfval;
	fitsum += curpool[i].selval;
	if (curpool[i].perfval == best)
	    besti = i;
    }
#endif

    for (i = 0; i < poolsize; i++)
    {
	sum = 0.0F;
	k = -1;			/* index of next Selected structure */

	ptr = fitsum * ran0();	/* spin the wheel one time */

	do
	{
	    k++;
	    sum += curpool[k].selval;
	}
	while (sum < ptr && k < poolsize - 1);

	samparr[i] = k;
    }

#ifndef ELITIST
    besti = -1;
#endif

#if 0
    /* randomly shuffle indices to new structures */
    for (i = 0; i < poolsize; i++)
    {
	j = randint(i, poolsize - 1);
	temp = samparr[j];
	samparr[j] = samparr[i];
	samparr[i] = temp;
    }
#endif

#if 1
    /* n-way tournament selection */
    for (i = 0; i < poolsize; i++)
    {
	for (k=0; k<1; k++)
	{
	    do
		j = randint(0, poolsize - 1);
	    while (j == i);
#ifdef ELITIST
	    if (i == besti)
		j = i;
#endif
	    if (curpool[j].perfval < curpool[i].perfval)
		break;
	}
	
	if (curpool[i].perfval < curpool[j].perfval)
	    samparr[i] = i;
	else
	    samparr[i] = j;
    }
#endif

    /* Form the new population */
    for (i = 0; i < poolsize; i++)
    {
	k = samparr[i];
	chaincpy(newpool[i].genome, curpool[k].genome, seqlen);
	newpool[i].perfval = curpool[k].perfval;
	newpool[i].evalflg = FALSE;
    }
}

void
                ga_randmove(void)
{
    int             i, nr, len, src, dest, from, to;
    float rotasel;

    do
	dest = randint(0, poolsize - 1);
    while (dest == besti);

    for (;;)
    {
	if (ran0() < moveprob[0] + moveprob[1])
	{
	    fragsamp(&to, &src, &from, &len);
	}
	else
	{
	    do
	    {
		len = 2;
		to = randint(0, seqlen - len);
		src = randint(0, nchn - 1);
		from = randint(0, chnlist[src].length - len);
	    } while (!chksst(src, from, to, len));
	}
	chaincpy(oldchn, newpool[dest].genome, seqlen);
	splice(newpool[dest].genome, chnlist[src].chain + from, to, len);
	if (!calcdist(newpool[dest].genome, TRUE))
	    break;
	chaincpy(newpool[dest].genome, oldchn, seqlen);
    }

    /* Set a rotamer to a random state (biased to most common rotamers) */
    i = randint(0, seqlen-1);

    if (nrots[seq[0][i]] > 1)
    {
	nr = 0;
	rotasel = ran0() * 100.0;
	while ((rotasel -= rotafreq[seq[0][i]][nr]) > 0.001)
	    nr++;
	
	if (nr >= nrots[seq[0][i]])
	    fail("nr > nrots!");
	
	newpool[dest].genome[i].rotnum = nr;
    }
    
    newpool[dest].evalflg = TRUE;
}

void
                mutate(float prob)
{
    int             n;

    if (prob > 0.0F)
	for (n = poolsize; n--;)
	    if (ran0() < prob)
		ga_randmove();
}

/* Randomly crossover pool */
void
                crossovr(void)
{
    int             i, j, p1;
    static RESATM          *temp, *save;

    if (!temp)
    {
	temp = calloc(seqlen, sizeof(RESATM));
	save = calloc(seqlen, sizeof(RESATM));
	
	if (temp == NULL || save == NULL)
	    fail("Out of memory in crossover!");
    }

    for (i = 0; i < poolsize - 1; i++)
	if (ran0() < crosrate)
	{
	    /* Attempt 1 single point crossover(s) */
	    for (j=0; j<1; j++)
	    {
		/* Single point crossover */
		p1 = randint(4, seqlen - 4);
		
		chaincpy(temp, newpool[i].genome, seqlen);
		chaincpy(save, newpool[i + 1].genome, seqlen);
		
		if (i != besti)
		{
		    splice(newpool[i].genome, newpool[i + 1].genome + p1, p1, seqlen - p1);
		    if (calcdist(newpool[i].genome, TRUE))
			chaincpy(newpool[i].genome, temp, seqlen);
		    else
			newpool[i].evalflg = TRUE;
		}
		
		if (i+1 != besti)
		{
		    splice(newpool[i + 1].genome, temp + p1, p1, seqlen - p1);
		    if (calcdist(newpool[i + 1].genome, TRUE))
			chaincpy(newpool[i + 1].genome, save, seqlen);
		    else
			newpool[i + 1].evalflg = TRUE;
		}
	    }
	    
	    i++;
	}
}

void
                statistics(Schema * pool)
{
    int             i;

    for (i = 0; i < poolsize; i++)
	if (pool[i].evalflg)
	{
	    (void) calcdist(pool[i].genome, FALSE);
	    pool[i].perfval = eval(pool[i].genome, FALSE);
	    pool[i].evalflg = FALSE;
	}
    avc_perf = best = worst = pool[0].perfval;

    for (i = 1; i < poolsize; i++)
    {
	avc_perf += pool[i].perfval;
	if (worst < pool[i].perfval)
	    worst = pool[i].perfval;
	if (best > pool[i].perfval)
	    best = pool[i].perfval;
    }
    avc_perf /= (float) poolsize;
}

void
                run_ga(void)
{
    int             i, gen, prevgen, opt_flag;
    float           prevbest, sigma;
    Schema         *temp;

    if (verboseflg)
	printf("Pool size = %d\n", poolsize);
    ga_init();

    prevbest = VBIG;

    statistics(curpool);

    if (verboseflg)
	printf("Initial: %g %g %g\n\n", worst, avc_perf, best);

    opt_flag = FALSE;

    for (gen = 1; gen <= maxngen; gen++)
	if (walltime() < MAXTIME)
	{
	    statistics(curpool);
	    
	    gaselect();
	    
	    statistics(newpool);
	    
	    if (verboseflg)
		printf("%d %g %g %g\n", gen, worst, avc_perf, best);
	    fflush(stdout);
	    
	    if (best < prevbest)
	    {
		prevbest = best;
		prevgen = gen;
		opt_flag = TRUE;
	    }
	    
	    crossovr();
	    mutate(mutrate);
	    
	    temp = newpool;
	    newpool = curpool;
	    curpool = temp;
	}
    
#if 0
    for (int i=0; i<poolsize; i++)
	writepdbensemble(newpool[i].genome, newpool[i].perfval);
#endif
}


/* Simulated Annealing */

void
                sa_randmove()
{
    int             i, len, src, from, to;

    do {
	if (ran0() < moveprob[0] + moveprob[1])
	{
	    fragsamp(&to, &src, &from, &len);
	}
	else
	{
	    do
	    {
		len = 2;
		if (mod_mode == FOLDM)
		    to = randint(0, seqlen-len);
		else
		    to = randint(REFSTART, MIN(seqlen, REFEND)-len);
		src = randint(0, nchn - 1);
		from = randint(0, chnlist[src].length - len);
	    }
	    while (!chksst(src, from, to, len));
	    lastmove = 2;
	}

	/* Only pick fragments within SVRs */
	if (mod_mode == MODELM)
	{
	    for (i = 0; i < len; i++)
		if (!svrflg[to+i] || (i > 0 && svrflg[to+i] != svrflg[to+i-1]))
		    break;
	    
	    if (i == len)
		break;
	}
    } while (mod_mode == MODELM);

    frag_combine(curchn, chnlist[src].chain + from, to, len);
}

void
                sc_randmove()
{
    int             i, nr;
    float rotasel;
    
    /* Set a rotamer to a random state (biased to most common rotamers) */
    i = randint(0, seqlen-1);
    
    if (nrots[seq[0][i]] > 1)
    {
	rotasel = ran0() * 100.0;
	for (nr = 0; nr < nrots[seq[0][i]]-1; nr++)
	    if ((rotasel -= rotafreq[seq[0][i]][nr]) <= 0.0F)
		break;
	
	curchn[i].rotnum = nr;
    }
}


/* Generate initial model by splicing gaps in template */
void
                initmodel(void)
{
    int             i, j, src, len, from, occ_start, occ_end;
    
    chaincpy(curchn, targchn, seqlen);

    if (svrflg[0])
    {
	/* Build N-terminus */
	for (i=1; i<seqlen; i++)
	    if (!svrflg[i])
		break;

	occ_start = i;

	len = i+1;

	do
	    src = randint(0, nchn - 1);
	while (chnlist[src].length < len);

	chaincpy(curchn, chnlist[src].chain, len);
        splice(curchn, targchn+occ_start, occ_start, seqlen-occ_start);
    }
    else
	occ_start = 0;

    if (svrflg[seqlen-1])
    {
	/* Build C-terminus */
	for (i=seqlen-2; i>0; i--)
	    if (!svrflg[i])
		break;

	occ_end = i;

	len = seqlen - i;

	do
	    src = randint(0, nchn - 1);
	while (chnlist[src].length < len);

	from = randint(0, chnlist[src].length - len);

        splice(curchn, chnlist[src].chain+from, occ_end, len);
    }
    else
	occ_end = seqlen-1;

    for (i=occ_start; i<occ_end; i++)
	if (!svrflg[i] && svrflg[i+1])
	{
	    for (j=i+1; j<occ_end; j++)
		if (!svrflg[j])
		    break;

	    len = j - i;

/*	    printf("len=%d at %d\n", len, i+1); */

	    do
		src = randint(0, nchn - 1);
	    while (chnlist[src].length < len);

	    from = randint(0, chnlist[src].length - len);

	    fixed_splice(curchn, chnlist[src].chain+from, i, len);
	    
	    i = j;
	}
    writemainpdb("debug1.pdb", curchn, 0, seqlen-1);
}


/* Refine model with FRAGFOLD energy */
void
                refine_model(void)
{
    int             i, nc, src, len, to, from;
    float           old_e, new_e;

    calcdist(curchn, FALSE);
    old_e = eval(curchn, FALSE);

    for (nc = 0; nc < 1000000; nc++)
    {
	chaincpy(oldchn, curchn, seqlen);

#if 1
	if (ran0() < 0.5)
	{
	    src = randint(0, nchn - 1);
	    len = randint(3, MIN(seqlen, 9));
	    to = randint(0, seqlen - len);
	    from = randint(0, chnlist[src].length - len);
	    splice(curchn, chnlist[src].chain + from, to, len);
	}
	else
	{
	    do {
		len = randint(3, MIN(seqlen, 50));
		to = randint(0, seqlen - len);
		for (i=to; i<to+len; i++)
		    if (!svrflg[i] || chnbrk[i])
			break;
	    } while (i != to+len);
	    splice(curchn, targchn + to, to, len);
	}
#endif

	calcdist(curchn, FALSE);

	new_e = eval(curchn, FALSE);
	
	if (new_e < old_e)
	{
	    old_e = new_e;
	    continue;
	}
	else
	    chaincpy(curchn, oldchn, seqlen);
    }

    buildschn(curchn, 0, seqlen-1);
    writepdb(curchn, 0.0);
}

/* Replica Exchange */
void
                repmc(void)
{
    float           ediff, av_ediff = 0.0, old_e, new_e, tmax, chi0,
	e_min = VBIG, delta, rmsd, lastrmsd, *rep_e, etemp, ttemp;
    int             i, j, cycle=0, *nsucc, *nuphill, nswaps, src, totalsteps = 0, minsteps;
    Replica         *ensemble;
    static unsigned optcount[100], swapcount[100];
    FILE *ifp;
    
    if (mod_mode == REFINEM || mod_mode == MODELM)
	chaincpy(curchn, targchn, seqlen);

    if (calcdist(curchn, TRUE))
    {
	if (verboseflg)
	    puts("WARNING: Clashes in starting conformation!");
	calcdist(curchn, FALSE);
    }
	
    /* Evaluate average energy change for random moves */
    for (i = 0; i < 500; i++)
    {
	chaincpy(oldchn, curchn, seqlen);
	for (;;)
	{
	    sa_randmove();
	    if (!calcdist(curchn, TRUE))
		break;
	    chaincpy(curchn, oldchn, seqlen);
	}

	new_e = eval(curchn, FALSE);
	
	if (i)
	{
	    ediff = fabs(new_e - old_e);
	    av_ediff += ediff;
	}

	if (new_e < e_min)
	    e_min = new_e;

	old_e = new_e;
    }

    av_ediff /= i+1;

    ensemble = (Replica *) calloc(poolsize, sizeof(Replica));
    nsucc = (int *) calloc(poolsize, sizeof(int));
    nuphill = (int *) calloc(poolsize, sizeof(int));

    if (mod_mode == 1)
	chaincpy(curchn, targchn, seqlen);

    if (opt_mode == 4)
	vtrange = 6;

    minsteps = 0;

    for (i=0; i<poolsize; i++)
    {
	ensemble[i].chn = (RESATM *) calloc(seqlen, sizeof(RESATM));

	chaincpy(ensemble[i].chn, curchn, seqlen);

	ensemble[i].energy = eval(curchn, FALSE);
//	ensemble[i].t = av_ediff * (poolsize - i) / poolsize;
	ensemble[i].t = pow(TRATIO, (float)i) * av_ediff * INITEMP;
    }

    if (verboseflg)
	puts("Replica Exchange Annealing...");

    while (totalsteps < MAXSTEPS && walltime() < MAXTIME)
    {
	if (opt_mode == 4 && minsteps > MAXSTEPS / seqlen / 2)
	{
	    printf("New Target Range = %d\n", vtrange = 6 + seqlen * 1.1 * totalsteps / MAXSTEPS);

	    for (i=0; i<poolsize; i++)
	    {
		chaincpy(curchn, ensemble[i].chn, seqlen);

		calcdist(curchn, TRUE);
		ensemble[i].energy = eval(curchn, FALSE);
	    }
	    
	    minsteps = 0;
	}

	for (i=0; i<poolsize; i++)
	{
	    chaincpy(curchn, ensemble[i].chn, seqlen);

	    for (;;)
	    {
		sa_randmove();
		if (!calcdist(curchn, TRUE))
		    break;
		chaincpy(curchn, ensemble[i].chn, seqlen);
	    }

	    sc_randmove();
	    
	    totalsteps++;

	    new_e = eval(curchn, FALSE);

	    if (new_e > ensemble[i].energy)
		nuphill[i]++;

	    minsteps++;

	    if (boltzmann(new_e - ensemble[i].energy, ensemble[i].t))
	    {
		if (new_e < e_min)
		    optcount[i]++;

		if (new_e < e_min)
		    e_min = new_e;

		if (new_e > ensemble[i].energy)
		    nsucc[i]++;

		/* printf("new_e : %d %f\n", nsucc+1, new_e); */
		ensemble[i].energy = new_e;
		chaincpy(ensemble[i].chn, curchn, seqlen);
		movefreq[lastmove]++;
	    }
	}

	/* Replica exchange... */
	nswaps = 0;

//	i = randint(0, poolsize-2);
//	j = randint(i+1, poolsize-1);

	for (i=0; i<poolsize; i++)
	    for (j=i+1; j<poolsize; j++)
	    {
		delta = (1.0/ensemble[i].t - 1.0/ensemble[j].t) * (ensemble[i].energy - ensemble[j].energy);
		if (delta >= 0.0 || ran0() <= exp(delta))
		{
		    ttemp = ensemble[i].t;
		    ensemble[i].t = ensemble[j].t;
		    ensemble[j].t = ttemp;
		    nswaps++;
		    swapcount[i]++;
		    swapcount[j]++;
		}
	    }
	

	if (verboseflg && !(++cycle % 100))
	{
	    printf("\nEner:");
	    for (i=0; i<poolsize; i++)
		printf(" %f", ensemble[i].energy);
	    printf("\nTemp:");
	    for (i=0; i<poolsize; i++)
		printf(" %f", ensemble[i].t / av_ediff);
	    printf("\nchi0:");
	    for (i=0; i<poolsize; i++)
		printf(" %f", (float) nsucc[i] / (float) nuphill[i]);
	    printf("\nOptC:");
	    for (i=0; i<poolsize; i++)
		printf(" %d", optcount[i]);
	    printf("\nSwap:");
	    for (i=0; i<poolsize; i++)
		printf(" %d", swapcount[i]);
	    printf("\nCycle %3d : Swaps = %d Steps = %d Emin = %g\n", cycle, nswaps, totalsteps, e_min);
	}
	
	/* Calculate current move accept probabilities */

#if 0
	moveprob[0] = (float) movefreq[0] / (movefreq[0] + movefreq[1] + movefreq[2]);
	moveprob[1] = (float) movefreq[1] / (movefreq[0] + movefreq[1] + movefreq[2]);
	moveprob[2] = (float) movefreq[2] / (movefreq[0] + movefreq[1] + movefreq[2]);

/*	printf("New move probabilities: %f %f %f\n", moveprob[0], moveprob[1], moveprob[2]); */
#endif

	fflush(stdout);
    }
}

/* Thermodynamic simulated annealing */
void
                thermoanneal(void)
{
    float           ediff, old_e, new_e, t, t0, e_min = VBIG, dst = 0.0, dct = 0.0;
    int             nsucc = 0, totalsteps = 0;

    if (mod_mode == REFINEM)
	chaincpy(curchn, targchn, seqlen);

    if (calcdist(curchn, TRUE))
    {
	if (verboseflg)
	    puts("WARNING: Clashes in starting conformation!");
	calcdist(curchn, FALSE);
    }

    old_e = eval(curchn, FALSE);

    t = t0 = 0.1;

    printf("T0 = %f\n", t0);

    if (verboseflg)
	puts("Thermodynamic Annealing...");

    while ((t > 1e-3 || totalsteps < MAXSTEPS) && walltime() < MAXTIME)
    {
	chaincpy(oldchn, curchn, seqlen);
	
	for (;;)
	{
	    sa_randmove();
	    if (!calcdist(curchn, TRUE))
		break;
	    chaincpy(curchn, oldchn, seqlen);
	}
	
	new_e = eval(curchn, FALSE);

	ediff = new_e - old_e;

	if (ran0() < exp(-ediff / t))
	{
	    nsucc++;
	    dct += ediff;
	    movefreq[lastmove]++;
	    if (new_e < e_min)
		e_min = new_e;
	    old_e = new_e;
	}
	else
	    chaincpy(curchn, oldchn, seqlen);

	if (ediff > 0.0)
	    dst -= ediff / t;

	if (dct >= 0.0 || dst == 0.0)
	    t = t0;
	else
	{
	    t = KAVALUE * (dct / dst);
	}

	if (verboseflg && totalsteps && totalsteps % 10000 == 0)
	{
	    printf("dct = %f dst = %f dct/dst = %f\n", dct, dst, dct/dst);
	    printf("T = %g Steps = %d Nsucc = %d Emin = %g\n", t, totalsteps, nsucc, e_min);
	    fflush(stdout);
	}

#if 0
	moveprob[0] = (float) movefreq[0] / (movefreq[0] + movefreq[1] + movefreq[2]);
	moveprob[1] = (float) movefreq[1] / (movefreq[0] + movefreq[1] + movefreq[2]);
	moveprob[2] = (float) movefreq[2] / (movefreq[0] + movefreq[1] + movefreq[2]);
	
	if (verboseflg)
	    printf("New move probabilities: %f %f %f\n", moveprob[0], moveprob[1], moveprob[2]);
#endif


	totalsteps++;
    }
}


/* Adaptive simulated annealing with reheating */
void
                anneal(void)
{
    float           ediff, av_ediff = 0.0, max_ediff = 0.0, old_e, new_e, t, t0, e_min = VBIG;
    int             i, nsucc = 0, totalsteps = 0, nuphill, nupmax = 0;

    if (mod_mode == REFINEM || mod_mode == MODELM)
	chaincpy(curchn, targchn, seqlen);

    if (calcdist(curchn, TRUE))
    {
	if (verboseflg)
	    puts("WARNING: Clashes in starting conformation!");
	calcdist(curchn, FALSE);
    }

    /* Evaluate energy of random configurations */
    for (i = 0; i < 500; i++)
    {
	chaincpy(oldchn, curchn, seqlen);
	for (;;)
	{
	    sa_randmove();
	    if (!calcdist(curchn, TRUE))
		break;
	    chaincpy(curchn, oldchn, seqlen);
	}

	new_e = eval(curchn, FALSE);

	if (mod_mode == REFINEM)
	    chaincpy(curchn, targchn, seqlen);

	if (i)
	{
	    ediff = fabsf(new_e - old_e);

	    if (ediff > max_ediff)
		max_ediff = ediff;

	    av_ediff += ediff;
	}

	old_e = new_e;
    }

    av_ediff /= i;

    t0 = av_ediff * INITEMP;

    if (verboseflg)
	puts("Adaptive Simulated Annealing with Reheating ...");

    nuphill = 0;
    while (totalsteps < MAXSTEPS && walltime() < MAXTIME)
    {
	chaincpy(oldchn, curchn, seqlen);
	
	for (;;)
	{
	    sa_randmove();
	    if (!calcdist(curchn, TRUE))
		break;
	    chaincpy(curchn, oldchn, seqlen);
	}

	sc_randmove();
	
	new_e = eval(curchn, FALSE);

	ediff = new_e - old_e;

	if (ediff > 0.0)
	{
	    nuphill++;
	    if (nuphill > nupmax)
		nupmax = nuphill;
	}
	else if (ediff < 0.0)
	    nuphill = 0;

	/* Compute new temperature based on uphill move history */
	t = 50.0 + (t0 - 50.0) * log(1.0 + nuphill) / 5.0;

	if (ran0() < exp(-ediff / t))
	{
	    nsucc++;
	    movefreq[lastmove]++;
	    if (new_e < e_min)
		e_min = new_e;
	    old_e = new_e;
	}
	else
	    chaincpy(curchn, oldchn, seqlen);

	if (verboseflg && totalsteps && totalsteps % 10000 == 0)
	{
	    printf("T = %g Steps = %d Nsucc = %d Nupmax = %d Emin = %g\n", t, totalsteps, nsucc, nupmax, e_min);
	    fflush(stdout);
	}

#if 0
	    moveprob[0] = (float) movefreq[0] / (movefreq[0] + movefreq[1] + movefreq[2]);
	    moveprob[1] = (float) movefreq[1] / (movefreq[0] + movefreq[1] + movefreq[2]);
	    moveprob[2] = (float) movefreq[2] / (movefreq[0] + movefreq[1] + movefreq[2]);
	    
	    if (verboseflg)
		printf("New move probabilities: %f %f %f\n", moveprob[0], moveprob[1], moveprob[2]);
#endif

	totalsteps++;
    }
}


/* Sample initial conformations from chain list & calculate potential component weights */
void
                estweights(void)
{
    int             nc, src, len, to, from, nsamp;
    float           refsd;

    if (mod_mode != REFINEM)
    {
	/* Generate fully extended chain to start without clashes */
	fullyextended(curchn);
    }
    else
	chaincpy(curchn, targchn, seqlen);

    nsamp = (wt_mode == 2) ? 50000 : 5000;

    /* Ignore first 100 conformations */
    for (nc = -100; nc < nsamp; nc++)
    {
	chaincpy(oldchn, curchn, seqlen);
	for (;;)
	{
	    src = randint(0, nchn - 1);
	    len = MIN(chnlist[src].length, seqlen);
	    to = randint(0, seqlen - len);
	    from = randint(0, chnlist[src].length - len);
	    splice(curchn, chnlist[src].chain + from, to, len);
	    if (!calcdist(curchn, TRUE))
		break;
	    chaincpy(curchn, oldchn, seqlen);
	}

	if (nc >= 0)
	{
/*	    writepdbensemble(curchn, last_steric); */
	    if (mqapcmd != NULL || !wt_mode)
	    {
		buildschn(curchn, 0, seqlen-1);
		return; /* Leave weights as they are! */
	    }
	    (void) eval(curchn, TRUE);
	}
    }

    /* Calculate weights relative to steric component */
    refsd = sqrt(steric_sumsq / nc - SQR(steric_sum / nc));

    if (refsd <= 1e-6F)
	return;

    if (srwt > 0.0F && sr_sumsq > 0.0)
	srwt *= refsd / sqrt(sr_sumsq / nc - SQR(sr_sum / nc));

    if (lrwt > 0.0F && lr_sumsq > 0.0)
	lrwt *= refsd / sqrt(lr_sumsq / nc - SQR(lr_sum / nc));

    if (solvwt > 0.0F && solv_sumsq > 0.0)
	solvwt *= refsd / sqrt(solv_sumsq / nc - SQR(solv_sum / nc));

    if (hbwt > 0.0F && hb_sumsq > 0.0)
	hbwt *= refsd / sqrt(hb_sumsq / nc - SQR(hb_sum / nc));

    if (compactwt > 0.0F && comp_sumsq > 0.0)
	compactwt *= refsd / sqrt(comp_sumsq / nc - SQR(comp_sum / nc));

    if (dswt > 0.0F && ds_sumsq > 0.0)
	dswt *= refsd / sqrt(ds_sumsq / nc - SQR(ds_sum / nc));

    if (fabs(srrwt) > 0.0F && srr_sumsq > 0.0)
	srrwt *= refsd / sqrt(srr_sumsq / nc - SQR(srr_sum / nc));
	
	if (fabs(lrrwt) > 0.0F && lrr_sumsq > 0.0)
	lrrwt *= refsd / sqrt(lrr_sumsq / nc - SQR(lrr_sum / nc));
}



void
                readxyz(char *buf, float *x, float *y, float *z)
{
    char            temp[9];

    temp[8] = '\0';
    strncpy(temp, buf, 8);
    *x = atof(temp);
    strncpy(temp, buf + 8, 8);
    *y = atof(temp);
    strncpy(temp, buf + 16, 8);
    *z = atof(temp);
}

/* Read target structure coords - return no. of zero occupancy residues */
int
                maketemplate(char *pdbname, char chainid, int model, RESATM *chn)
{
    FILE           *pfp;
    char            buf[160], whichatm[MAXSEQLEN], inscode;
    int             atc, i, j, k, nres, resid, lastid, nseg;
    float           dv, maxd, sx, sy;
    Vector          cca, nca, xx, yy;
    short           at1, at2;
    float           x[MAXSEQLEN][5], y[MAXSEQLEN][5], z[MAXSEQLEN][5];
    static const float caca_maxdist[5] = { 3.90, 7.34, 10.81, 14.19, 17.58 };

    if (verboseflg)
	printf("Reading %s chain %c ...\n", pdbname, chainid);
    pfp = fopen(pdbname, "r");
    if (!pfp)
	fail("maketemplate: Cannot open PDB file!");

    for (i = 0; i < MAXSEQLEN; i++)
	whichatm[i] = 0;

    if (model >= 1)
    {
	while (!feof(pfp))
	{
	    if (fgets(buf, 160, pfp) == NULL)
		fail("maketemplate: cannot find MODEL!");
	    if (strncmp(buf, "MODEL", 5))
		continue;
	    sscanf(buf + 5, "%d", &i);
	    if (i == model)
		break;
	}
    }

    i = -1;
    nsvr = nseg = 0;
    while (!feof(pfp))
    {
	if (fgets(buf, 160, pfp) == NULL)
	    break;
	if (!strncmp(buf, "END", 3))
	    break;
	if (!strncmp(buf, "TER", 3))
	{
	    nseg++;
	    continue;
	}
	/* printf("%d %s\n",i,buf); */
	if (strncmp(buf, "ATOM", 4) || (buf[21] != chainid && !(buf[21] == 'A' && chainid == ' ')) || (buf[16] != ' ' && buf[16] != 'A'))
	    continue;
	for (atc = 0; atc <= CATOM; atc++)
	    if (!strncmp(buf + 13, atmnames[atc], 2))
	    {
		inscode = buf[26];
		buf[26] = ' ';
		sscanf(buf + 22, "%d", &resid);
		if (atc == NATOM)
		{
		    ++i;
		    readxyz(buf + 30, &x[i][atc], &y[i][atc], &z[i][atc]);
		    if (strlen(buf) > 56)
			svrflg[i] = buf[56] != '1';
		    else
			svrflg[i] = 0;
		    if (svrflg[i])
			nsvr++;
		    segidx[i] = nseg;
		    if (!i)
			firstid = resid;
		    else if (inscode == ' ' && resid - lastid > 1)
		    {
			if (SQR(x[i][NATOM] - x[i - 1][NATOM]) + SQR(y[i][NATOM] - y[i - 1][NATOM]) + SQR(z[i][NATOM] - z[i - 1][NATOM]) > 15.0)
			{
			    printf("WARNING: Skipping %d missing residues!\n", resid - lastid - 1);
			    for (k = 0; k < resid - lastid - 1; k++)
				i++;
			}
		    }
		    lastid = resid;
		}
		else if (i >= 0)
		    readxyz(buf + 30, &x[i][atc], &y[i][atc], &z[i][atc]);
		whichatm[i] |= 1 << atc;
	    }
    }
    fclose(pfp);

    nres = i + 1;

    /* Check atoms */
    for (i = 0; i < nres; i++)
	if (!svrflg[i])
	{
	    if (!(whichatm[i] & (1 << NATOM)))
	    {
		printf("FATAL: Missing N atom in %d!\n", i + 1);
		exit(1);
	    }
	    if (!(whichatm[i] & (1 << CAATOM)))
	    {
		printf("WARNING: Missing CA atom in %d!\n", i + 1);
	    }
	    if (!(whichatm[i] & (1 << CBATOM)))
	    {
		if (!(whichatm[i] & (1 << CAATOM)) || !(whichatm[i] & (1 << CATOM)) || !(whichatm[i] & (1 << NATOM)))
		{
		    /* Not much left of this residue! */
		    printf("WARNING: Missing main-chain atom in %d!\n", i + 1);
		    continue;
		}
		
		/* Reconstruct CB atom */
		nca[0] = x[i][CAATOM] - x[i][NATOM];
		nca[1] = y[i][CAATOM] - y[i][NATOM];
		nca[2] = z[i][CAATOM] - z[i][NATOM];
		cca[0] = x[i][CAATOM] - x[i][CATOM];
		cca[1] = y[i][CAATOM] - y[i][CATOM];
		cca[2] = z[i][CAATOM] - z[i][CATOM];
		vecadd(xx, nca, cca);
		vecprod(yy, nca, cca);
		sx = CACBDIST * cosf(TETH_ANG) / sqrtf(dotprod(xx, xx));
		sy = CACBDIST * sinf(TETH_ANG) / sqrtf(dotprod(yy, yy));
		x[i][CBATOM] = x[i][CAATOM] + xx[0] * sx + yy[0] * sy;
		y[i][CBATOM] = y[i][CAATOM] + xx[1] * sx + yy[1] * sy;
		z[i][CBATOM] = z[i][CAATOM] + xx[2] * sx + yy[2] * sy;
		whichatm[i] |= 1 << CBATOM;
	    }
	}

    if (nres != seqlen)
	fail("Sequence length mismatch in PDB file!");
    
    for (i = 0; i < nres; i++)
    {
	chn[i].n[0] = x[i][NATOM];
	chn[i].n[1] = y[i][NATOM];
	chn[i].n[2] = z[i][NATOM];
	chn[i].ca[0] = x[i][CAATOM];
	chn[i].ca[1] = y[i][CAATOM];
	chn[i].ca[2] = z[i][CAATOM];
	chn[i].c[0] = x[i][CATOM];
	chn[i].c[1] = y[i][CATOM];
	chn[i].c[2] = z[i][CATOM];
	chn[i].o[0] = x[i][OATOM];
	chn[i].o[1] = y[i][OATOM];
	chn[i].o[2] = z[i][OATOM];
	chn[i].cb[0] = x[i][CBATOM];
	chn[i].cb[1] = y[i][CBATOM];
	chn[i].cb[2] = z[i][CBATOM];
    }

    for (i = 0; i < nres; i++)
	chn[i].rotnum = 0;

    buildschn(chn, 0, seqlen-1);

    if (verboseflg)
	printf("Target NRES = %d\n", nres);
    
    /* Calculate interatomic distance templates */
    at1 = at2 = CBATOM;
    for (i = 0; i < nres; i++)
	for (j = i + 1; j < nres; j++)
	{
	    dv = SQR(x[i][at1] - x[j][at2]) + SQR(y[i][at1] - y[j][at2]) + SQR(z[i][at1] - z[j][at2]);
	    dsqmat[i][j][CB_CB] = dsqmat[j][i][CB_CB] = dv;
	}

    at1 = at2 = CAATOM;
    for (i = 0; i < nres; i++)
	for (j = i + 1; j < nres; j++)
	{
	    dv = SQR(x[i][at1] - x[j][at2]) + SQR(y[i][at1] - y[j][at2]) + SQR(z[i][at1] - z[j][at2]);
	    dsqmat[i][j][CA_CA] = dsqmat[j][i][CA_CA] = dv;
	}

    for (i = 0; i < nres; i++)
	tcbooi[i] = 0;

    /* Compute CB-Ooi numbers */
    for (i = 0; i < nres; i++)
	for (j = i + 1; j < nres; j++)
	    if (dsqmat[i][j][CB_CB] > 0.0F && dsqmat[i][j][CB_CB] <= SQR(OOICUTOFF))
	    {
		tcbooi[i]++;
		tcbooi[j]++;
	    }

#if 0
    /* Check for gaps which cannot be spanned */
    if (nsvr)
	for (i = 0; i < seqlen; i++)
	    for (j = i+2; j<seqlen; j++)
		if (!svrflg[i] && !svrflg[j])
		{
		    if (j-i < 6)
			maxd = caca_maxdist[j-i-1];
		    else
			maxd = 3.55F * (j - i);
		    dv = dist(chn[i].ca, chn[j].ca);
		    if (dv > maxd)
		    {
			for (k=i; k<=j; k++)
			    svrflg[k] = 1;
			printf("Residues too far apart (unbuilding) : %d - %d\n", i+1, j+1);
		    }
		}
#endif

    for (nsvr = i = 0; i < seqlen; i++)
	if (svrflg[i])
	    nsvr++;

    if (mod_mode != MODELM)
	return nsvr;
    
    for (i = 0; i < seqlen; i++)
	if (svrflg[i])
	    break;
    
    if (i < 9)
	while (--i >= 0)
	    svrflg[i] = 1;

    for (i = seqlen-1; i; i--)
	if (svrflg[i])
	    break;

    if (i > seqlen-9)
	while (++i < seqlen)
	    svrflg[i] = 1;
    
    for (i = 0; i < seqlen; i++)
	if (svrflg[i])
	    svrflg[i] = 2;
	else
	    break;

    for (i = seqlen-1; i; i--)
	if (svrflg[i])
	    svrflg[i] = 3;
	else
	    break;

    if (nsvr > 0 && verboseflg)
    {
	puts("SVR MAP:");
	for (i=0; i<seqlen; i++)
	    putchar('0' + svrflg[i]);
	putchar('\n');
    }
    
    return nsvr;
}

char *getnextp(FILE *pfp)
{
    static char pv[256];

    while (!feof(pfp))
    {
	if (fscanf(pfp, "%s", pv) != 1)
	    return NULL;
	if (pv[0] == '#')
	{
	    fgets(pv, 256, pfp);
	    continue;
	}

	if (verboseflg)
	    printf("Read parameter: %s\n", pv);
	return pv;
    }

    return NULL;
}

/* Read parameter file if present */
void readparams(char *fname)
{
    char *vname, *value;
    FILE *pfp;

    pfp = fopen(fname, "r");

    if (!pfp)
	fail("readparams: cannot open parameter file!");

    if (verboseflg)
	puts("\nReading parameter file...");

    while ((vname = getnextp(pfp)))
    {
	value = "N/A";

	if (!strcasecmp(vname, "ALNFILE"))
	{
	    value = getnextp(pfp);
	    if (value != NULL)
		strcpy(alnfname, value);
	    if (verboseflg)
		puts(value);
	}
	else if (!strcasecmp(vname, "CONFILE"))
	{
	    value = getnextp(pfp);
	    if (value != NULL)
		strcpy(confname, value);
	    if (verboseflg)
		puts(value);
	}
	else if (!strcasecmp(vname, "HBFILE"))
	{
	    value = getnextp(pfp);
	    if (value != NULL)
		strcpy(hbfname, value);
	    if (verboseflg)
		puts(value);
	}
	else if (!strcasecmp(vname, "MQAP"))
	{
	    value = getnextp(pfp);
	    if (value != NULL)
		mqapcmd = strdup(value);
	    puts(value);
	}
	else if (!strcasecmp(vname, "SRWT"))
	{
	    value = getnextp(pfp);
	    if (value != NULL)
		srwt = atof(value);
	}
	else if (!strcasecmp(vname, "LRWT"))
	{
	    value = getnextp(pfp);
	    if (value != NULL)
		lrwt = atof(value);
	}
	else if (!strcasecmp(vname, "SOLVWT"))
	{
	    value = getnextp(pfp);
	    if (value != NULL)
		solvwt = atof(value);
	}
	else if (!strcasecmp(vname, "HBWT"))
	{
	    value = getnextp(pfp);
	    if (value != NULL)
		hbwt = atof(value);
	}
	else if (!strcasecmp(vname, "COMPACTWT"))
	{
	    value = getnextp(pfp);
	    if (value != NULL)
		compactwt = atof(value);
	}
	else if (!strcasecmp(vname, "STERICWT"))
	{
	    value = getnextp(pfp);
	    if (value != NULL)
		stericwt = atof(value);
	}
	else if (!strcasecmp(vname, "DSWT"))
	{
	    value = getnextp(pfp);
	    if (value != NULL)
		dswt = atof(value);
	}
	else if (!strcasecmp(vname, "RRWT"))
	{
	    value = getnextp(pfp);
	    if (value != NULL)
		srrwt = atof(value);
		lrrwt = atof(value);
	}
	else if (!strcasecmp(vname, "SRRWT"))
	{
	    value = getnextp(pfp);
	    if (value != NULL)
		srrwt = atof(value);
	}else if (!strcasecmp(vname, "LRRWT"))
	{
	    value = getnextp(pfp);
	    if (value != NULL)
		lrrwt = atof(value);
	}
	else if (!strcasecmp(vname, "TARGWT"))
	{
	    value = getnextp(pfp);
	    if (value != NULL)
		targwt = atof(value);
	}
	else if (!strcasecmp(vname, "LRWT"))
	{
	    value = getnextp(pfp);
	    if (value != NULL)
		lrwt = atof(value);
	}
	else if (!strcasecmp(vname, "STERMODE"))
	{
	    value = getnextp(pfp);
	    if (value != NULL)
	    {
		if (!strncasecmp(value, "CBONLY", 2))
		    global_ster_mode = 0;
		else if (!strncasecmp(value, "CENTROID", 3))
		    global_ster_mode = 1;
		else if (!strncasecmp(value, "ALLATOM", 3))
		    global_ster_mode = 2;
		else if (!strncasecmp(value, "DFIRE", 3))
		    global_ster_mode = 3;
		else
		    fail("Unknown steric mode specified in parameter file!");
	    }
	}
	else if (!strcasecmp(vname, "MAXFRAGS"))
	{
	    value = getnextp(pfp);
	    if (value != NULL)
		MAXFRAGS = atoi(value);
	}
	else if (!strcasecmp(vname, "MAXFRAGS2"))
	{
	    value = getnextp(pfp);
	    if (value != NULL)
		MAXFRAGS2 = atoi(value);
	}
	else if (!strcasecmp(vname, "MINFRAGS2LEN"))
	{
	    value = getnextp(pfp);
	    if (value != NULL)
		MINFRAGS2LEN = atoi(value);
	}
	else if (!strcasecmp(vname, "MAXFRAGS2LEN"))
	{
	    value = getnextp(pfp);
	    if (value != NULL)
		MAXFRAGS2LEN = atoi(value);
	}
	else if (!strcasecmp(vname, "POOLSIZE"))
	{
	    value = getnextp(pfp);
	    if (value != NULL)
		poolsize = atoi(value);
	}
	else if (!strcasecmp(vname, "MAXNGEN"))
	{
	    value = getnextp(pfp);
	    if (value != NULL)
		maxngen = atoi(value);
	}
	else if (!strcasecmp(vname, "MUTRATE"))
	{
	    value = getnextp(pfp);
	    if (value != NULL)
		mutrate = atof(value);
	}
	else if (!strcasecmp(vname, "TRATIO"))
	{
	    value = getnextp(pfp);
	    if (value != NULL)
		TRATIO = atof(value);
	}
	else if (!strcasecmp(vname, "CROSRATE"))
	{
	    value = getnextp(pfp);
	    if (value != NULL)
		crosrate = atof(value);
	}
	else if (!strcasecmp(vname, "MAXSTEPS"))
	{
	    value = getnextp(pfp);
	    if (value != NULL)
		MAXSTEPS = atoi(value);
	}
	else if (!strcasecmp(vname, "REFSTART"))
	{
	    value = getnextp(pfp);
	    if (value != NULL)
		REFSTART = atoi(value);
	}
	else if (!strcasecmp(vname, "REFEND"))
	{
	    value = getnextp(pfp);
	    if (value != NULL)
		REFEND = atoi(value);
	}
	else if (!strcasecmp(vname, "MAXTIME"))
	{
	    value = getnextp(pfp);
	    if (value != NULL)
		MAXTIME = atoi(value);
	}
	else if (!strcasecmp(vname, "INITEMP"))
	{
	    value = getnextp(pfp);
	    if (value != NULL)
		INITEMP = atof(value);
	}
	else if (!strcasecmp(vname, "KAVALUE"))
	{
	    value = getnextp(pfp);
	    if (value != NULL)
		KAVALUE = atof(value);
	}
	else if (!strcasecmp(vname, "REFRANGE"))
	{
	    value = getnextp(pfp);
	    if (value != NULL)
		REFRANGE = atof(value);
	}
	else if (!strcasecmp(vname, "WTMODE"))
	{
	    value = getnextp(pfp);
	    if (value != NULL)
	    {
		if (!strncasecmp(value, "ABS", 3))
		    wt_mode = 0;
		else if (!strncasecmp(value, "STD", 3))
		    wt_mode = 1;
		else if (!strncasecmp(value, "EST", 3))
		    wt_mode = 2;
		else
		    fail("Unknown weighting mode specified in parameter file!");
	    }
	}
	else if (!strcasecmp(vname, "OPTMETHOD"))
	{
	    value = getnextp(pfp);
	    if (value != NULL)
	    {
		if (!strcasecmp(value, "SA"))
		    opt_mode = 0;
		else if (!strcasecmp(value, "GA"))
		    opt_mode = 1;
		else if (!strcasecmp(value, "TSA"))
		    opt_mode = 2;
		else if (!strcasecmp(value, "REPMC"))
		    opt_mode = 3;
		else if (!strcasecmp(value, "VTREPMC"))
		    opt_mode = 4;
		else
		    fail("Unknown optimization method specified in parameter file!");
	    }
	}
	else if (!strcasecmp(vname, "FRAGSEL"))
	{
	    value = getnextp(pfp);
	    if (value != NULL)
	    {
		if (!strcasecmp(value, "PAIR"))
		    fragsel = 0;
		else if (!strcasecmp(value, "RMSD"))
		    fragsel = 1;
		else
		    fail("Unknown fragment selection method specified in parameter file!");
	    }
	}
	else if (!strcasecmp(vname, "MODE"))
	{
	    value = getnextp(pfp);
	    if (value != NULL)
	    {
		if (!strncasecmp(value, "FOLD", 3))
		    mod_mode = FOLDM;
		else if (!strncasecmp(value, "REFINE", 3))
		    mod_mode = REFINEM;
		else if (!strncasecmp(value, "MODEL", 3))
		    mod_mode = MODELM;
		else if (!strncasecmp(value, "EVAL", 4))
		    mod_mode = EVALM;
		else
		    fail("Unknown mode specified in parameter file!");
	    }
	}
	else
	    fail("Unknown keyword (%s) in parameter file!", vname);

	if (value == NULL)
	{
	    fprintf(stderr, "%s - missing parameter value in parameter file!\n", vname);
	    exit(-1);
	}
    }

    fclose(pfp);
}

	
int main(int argc, char **argv)
{
    int             a, aa, b, i, j, k, optstart = 0;
    char            fpotname[160], dfirename[160], adistname[160];
    float eqd, dsd, maxd, prob, e;
    const int       bigendv=1;
    FILE           *ifp;
    
    printf("FRAGFOLD Version %s (Last edit %s)\n", FragfoldVersion, Last_Edit_Date);
    printf("Motif-constrained Protein Folding Program\n");
    printf("Build date : %s\n",__DATE__);
    printf("Copyright (C) 1994/2000 David T. Jones\n\n");

    if (argc < 2)
	fail("usage : %s paramfile [outpdbfn] [pdbfname] [chainid] [options]\n", argv[0]);

    randomise();

    /* Parse options at end of command line */

    for (i=2; i<argc; i++)
    {
	if (argv[i][0] == '-')
	{
	    if (!optstart)
		optstart = i;
	    switch(argv[i][1])
	    {
	    case 'q':
		verboseflg = FALSE;
		break;
	    }
	}
    }

    if (optstart)
	argc = optstart;

    readparams(argv[1]);

    if (verboseflg)
	puts("\nReading mean force potential tables...");

    if (getenv("FDATA_DIR"))
    {
	strcpy(fpotname, getenv("FDATA_DIR"));
	if (fpotname[strlen(fpotname)-1] != '/')
	    strcat(fpotname, "/");
	strcpy(dfirename, fpotname);
	strcpy(adistname, fpotname);
    }
    else
	fpotname[0] = dfirename[0] = adistname[0] = '\0';

    strcat(fpotname, "fragfpot.dat");
	
    ifp = fopen(fpotname, "rb");
    if (!ifp)
	fail("main: cannot open fragfpot.dat!");
    fread(sr_de, 1, sizeof(sr_de), ifp);
    fread(lr_de, 1, sizeof(lr_de), ifp);
    fread(acc_de, 1, sizeof(acc_de), ifp);
    fread(ooi_de, 1, sizeof(ooi_de), ifp);
    fclose(ifp);

    strcat(dfirename, "dfirepot.dat");
    ifp = fopen(dfirename, "rb");
    if (!ifp)
	fail("main: cannot open dfirepot.dat!");
    fread(dfire, 1, sizeof(dfire), ifp);
    fclose(ifp);

    strcat(adistname, "atomdist.dat");
    ifp = fopen(adistname, "rb");
    if (!ifp)
	fail("main: cannot open atomdist.dat!");
    fread(atomdistsq, 1, sizeof(atomdistsq), ifp);
    fclose(ifp);

    /* Swap bytes if running on bigendian machine */
    if (*(char *)&bigendv != 1)
    {
	byteswap(sr_de, sizeof(sr_de)/4);
	byteswap(lr_de, sizeof(lr_de)/4);
	byteswap(acc_de, sizeof(acc_de)/4);
	byteswap(ooi_de, sizeof(ooi_de)/4);
	byteswap(dfire, sizeof(dfire)/4);
	byteswap(atomdistsq, sizeof(atomdistsq)/4);
    }

    for (i=0; i<167; i++)
	for (j=0; j<167; j++)
	    for (k=0; k<20; k++)
		if (dfire[i][j][k] == 999.999F)
		    dfire[i][j][k] = 0.5F;

    for (i = 0; i < TOPOMAX; i++)
	for (a = 0; a < 4; a++)
	    for (b = 0; b < 4; b++)
		diffdist[i][a][b] = (maxdist[i][a][b] - mindist[i][a][b]);

    for (i=0; i<167; i++)
	for (j=0; j<167; j++)
	    for (k=0; k<3; k++)
		atomdistsq[i][j][k] = SQR(atomdistsq[i][j][k]);

#if 0

    for (i=OOIMIN; i<=OOIMAX; i+=1)
	printf("%d %f\n", i, -ooi_de[TRP][MAX(0, MIN(OOIMAX, i) - OOIMIN) / OOIDIV]);

    putchar('\n');

    for (i=0; i<40; i++)
	printf("%d %f\n", i, -pairpot(CYS, CYS, CBATOM, CBATOM, 50, (float)i));
    putchar('\n');
    for (i=0; i<40; i++)
	printf("%d %f\n", i, -pairpot(LEU, LEU, CBATOM, CBATOM, 50, (float)i));
    putchar('\n');
    for (i=0; i<40; i++)
	printf("%d %f\n", i, -pairpot(GLU, ARG, CBATOM, CBATOM, 50, (float)i));
    putchar('\n');
    for (i=0; i<40; i++)
	printf("%d %f\n", i, -pairpot(GLU, GLU, CBATOM, CBATOM, 50, (float)i));
    putchar('\n');
    exit(0);
#endif

#ifdef RASMOL
    /* Start rasmol if available */
    rfp = popen("rasmol", "w");
#endif

    ifp = fopen(alnfname, "r");
    if (!ifp)
	fail("main: unable to open aln file!");
    fgets(buf, sizeof(buf), ifp);
    fscanf(ifp, "%d\n", &nseqs);
    fgets(tplt_ss, sizeof(tplt_ss), ifp);
    seqlen = strlen(tplt_ss) - 1;

    if (seqlen < 10)
	fail("Target sequence length < 10!!");

    for (i=0; i<seqlen; i++)
    {
	if (tplt_ss[i] == 'C')
	    tplt_ss[i] = COIL;
	else if (tplt_ss[i] == 'H')
	    tplt_ss[i] = HELIX;
	else if (tplt_ss[i] == 'E')
	    tplt_ss[i] = STRAND;
	else if (tplt_ss[i] == 'h')
	    tplt_ss[i] = HELIX|COIL;
	else if (tplt_ss[i] == 'e')
	    tplt_ss[i] = STRAND|COIL;
	else
	    tplt_ss[i] = HELIX|STRAND|COIL;
    }
	
    fgets(buf, sizeof(buf), ifp);

    if (strlen(buf) != seqlen+1)
	fail("main: mismatching line length in aln file!");

    for (i=0; i<9; i++)
	ds_from[i] = -1;
    for (i=0; i<seqlen; i++)
	if (isdigit(buf[i]))
	{
	    j = buf[i] - '1';
	    buf[i] = 'C';
	    if (ds_from[j] < 0)
		ds_from[j] = i;
	    else
		ds_to[j] = i;
	    if (j+1 > n_ds)
		n_ds = j+1;
	}

    seq = (char **) malloc(nseqs * sizeof(char *));

    if (!seq)
	fail("main: out of memory!");
    seq[0] = malloc(seqlen);
    if (!seq[0])
	fail("main: out of memory!");

    memcpy(seq[0], buf, seqlen);

    for (i = 1; i < nseqs; i++)
    {
	seq[i] = malloc(seqlen);
	if (!seq[i])
	    fail("main: out of memory!");
	if (!fgets(buf, MAXSEQLEN, ifp))
	    fail("main: error while reading aln file!");
	if (strlen(buf) != seqlen+1)
	    fail("main: mismatching line length in aln file!");
	memcpy(seq[i], buf, seqlen);
    }

    for (i = 0; i < nseqs; i++)
	for (j = 0; j < seqlen; j++)
	{
	    if (!isalpha(seq[i][j]))
	    {
		seq[i][j] = GAP;
		gaps[j]++;
	    }
	    else
	    {
		aa = aanum(seq[i][j]);
		switch (aa)
		{
		case 20:
		    /* ASP is more common than ASN! */
		    seq[i][j] = ASP;
		    break;
		case 21:
		    /* GLU is more common than GLN! */
		    seq[i][j] = GLU;
		    break;
		case 22:
		    seq[i][j] = UNK;
		    break;
		default:
		    seq[i][j] = aa;
		    break;
		}
	    }
	}

    fclose(ifp);

    if (confname[0])
    {
	conmat = allocmat(seqlen, seqlen, sizeof(float));    
	dcomat = allocmat(seqlen, seqlen, sizeof(float));    

	rrcbflag = 0;
	ifp = fopen(confname, "r");
	if (!ifp)
	    fail("main: unable to open con file!");
	
	while (!feof(ifp))
	{
	    if (!fgets(buf, sizeof(buf), ifp))
		break;

	    if (!strncmp(buf, "PFRMAT RR", 9))
	    {
		rrcbflag = 1;
		continue;
	    }

	    if (!strncmp(buf, "RR CONSTR", 9))
	    {
		rrcbflag = 2;
		continue;
	    }

	    if (isalpha(buf[0]))
		continue;

	    if (rrcbflag == 2)
	    {
		if (sscanf(buf, "%d%d%f%f", &i, &j, &eqd, &dsd) != 4)
		    continue;
		conmat[i-1][j-1] = conmat[j-1][i-1] = eqd;
		dcomat[i-1][j-1] = dcomat[j-1][i-1] = dsd;
	    }
	    else
	    {
		if (sscanf(buf, "%d%d%*s%f%f", &i, &j, &maxd, &prob) != 4)
		    continue;

//		conmat[i-1][j-1] += log((1.0/prob - 1.0) / (0.9/0.1));
		conmat[i-1][j-1] += log(1.0 - prob);
		conmat[j-1][i-1] = conmat[i-1][j-1];
		dcomat[i-1][j-1] = MAX(dcomat[i-1][j-1], maxd);
		dcomat[j-1][i-1] = dcomat[i-1][j-1];
	    }
	}
	
	fclose(ifp);
    }

    if (hbfname[0])
    {
	hbmat = allocmat(seqlen, seqlen, sizeof(float));    

	ifp = fopen(hbfname, "r");
	if (!ifp)
	    fail("main: unable to open HB file!");
	
	while (!feof(ifp))
	{
	    if (!fgets(buf, sizeof(buf), ifp))
		break;

	    if (isalpha(buf[0]) || ispunct(buf[0]))
		continue;

	    if (sscanf(buf, "%d%d%*s%*s%f", &i, &j, &prob) != 3)
		continue;
	    
	    hbmat[i-1][j-1] = prob;
	}
	
	fclose(ifp);
    }
    
    build_emats();

    exprad = pow(3 * seqlen / 12.56637 / 0.0065, 1.0 / 3.0);
    expsurf = 0.0F;
    for (i = 0; i < seqlen; i++)
	expsurf += molwt[seq[0][i]];
    expsurf = 6.32 * pow(expsurf, 0.728);
    expmean = pow(66.6434 * seqlen, 1.0 / 3.0);

    if (verboseflg)
	printf("Expected radius = %f\nExpected mean distance = %f\nExpected surface area = %f\n", exprad, expmean, expsurf);

    dsqmat = allocmat(seqlen, seqlen, sizeof(**dsqmat));

    bestchn = calloc(seqlen, sizeof(RESATM));

    if (bestchn == NULL)
	fail("main: calloc bestchn failed!");

    curchn = calloc(seqlen, sizeof(RESATM));

    for (i=0; i<seqlen; i++)
	curchn[i].sc = NULL;

    if (curchn == NULL)
	fail("main: calloc curchn failed!");

    oldchn = calloc(seqlen, sizeof(RESATM));

    if (oldchn == NULL)
	fail("main: calloc oldchn failed!");

    targchn = calloc(seqlen, sizeof(RESATM));
    
    if (targchn == NULL)
	fail("main: calloc targchn failed!");


    if (verboseflg)
	puts("\nReading fold library...");
    readchn();

    if (argc == 4)
    {
	maketemplate(argv[3], ' ', -1, targchn);

	for (i = 0; i < seqlen; i++)
	    chnbrk[i] = 0;
	
	for (i = 0; i < seqlen-1; i++)
	    if (distsq(targchn[i].ca, targchn[i+1].ca) > 16.0)
		chnbrk[i] = chnbrk[i+1] = 1;

	CenterCoords(targchn, seqlen);

	cb_targ = allocmat(seqlen, seqlen, sizeof(float));

	for (i = 0; i < seqlen; i++)
	    for (j = 0; j < seqlen; j++)
		if (!svrflg[i] && !svrflg[j] && segidx[i] == segidx[j])
		    cb_targ[i][j] = sqrtf(dsqmat[i][j][CB_CB]);
		else
		    cb_targ[i][j] = 0.0;
    }
    else if (argc >= 5)
    {
	if (strlen(argv[4]) != 1)
	    sscanf(argv[4] + 1, "%d", &a);
	else
	    a = -1;

	maketemplate(argv[3], argv[4][0], a, targchn);

	for (i = 0; i < seqlen; i++)
	    chnbrk[i] = 0;
	
	for (i = 0; i < seqlen-1; i++)
	    if (distsq(targchn[i].ca, targchn[i+1].ca) > 16.0)
		chnbrk[i] = chnbrk[i+1] = 1;

	CenterCoords(targchn, seqlen);

	cb_targ = allocmat(seqlen, seqlen, sizeof(float));

	for (i = 0; i < seqlen; i++)
	    for (j = 0; j < seqlen; j++)
		if (!svrflg[i] && !svrflg[j] && segidx[i] == segidx[j])
		    cb_targ[i][j] = sqrtf(dsqmat[i][j][CB_CB]);
		else
		    cb_targ[i][j] = 0.0;
    }

    if (argc > 2)
    {
	outpdbn = argv[2];
	if (fopen(argv[2], "r") != NULL)
	    fail("main: output PDB file exists!");
    }

    if (!cb_targ && mod_mode)
	fail("No target PDB file was specified in modelling/refinement mode!");

    readrots();
        
    estweights();

    if (verboseflg)
	printf("SRWT %f\nLRWT %f\nSOLVWT %f\nHBWT %f\nCOMPACTWT %f\nSTERICWT %f\nDSWT %f\nSRRWT %f\nLRRWT %f\nTARGWT %f\n\n", srwt, lrwt, solvwt, hbwt, compactwt, stericwt, dswt, srrwt, lrrwt, targwt);

    if (wt_mode == 2)
	exit(0);

    if (verboseflg)
	printf("MOD MODE = %d\n", mod_mode);

    /* If we are not modelling then calculate energy for target structure */

    if (cb_targ)
    {
	calcdist(targchn, FALSE);
	sc_opt(dsqmat, targchn, global_ster_mode);
/*	sc_anneal(dsqmat, targchn, global_ster_mode); */
    }

    if (argc >= 4 && mod_mode != MODELM && verboseflg)
    {
	puts("\nEnergy sums for target structure:");

	e_model(dsqmat, targchn, TRUE);

	if (mod_mode == EVALM)
	{
	    printf("%f %f %f %f %f %f %f %f %f\n", srwt, lrwt, solvwt, hbwt, compactwt, stericwt, dswt, srrwt, lrrwt);
	    printf("%f %f %f %f %f %f %f %f %f\n", last_sr, last_lr, last_solv, last_hbond, last_compact, last_steric, last_ds, last_srr, last_lrr);
	    exit(0);
	}
    }

    cb_last = allocmat(seqlen, seqlen, sizeof(float));

    if (MAXFRAGS > 0)
	mkfragtb();

    if (MAXFRAGS2 > 0)
	mkfragtb2();

    (void) cputime();

    if (opt_mode == 1)
	run_ga();
    else
    {
	if (mod_mode == MODELM)
	    initmodel();

	if (opt_mode == 2)
	    thermoanneal();
	else if (opt_mode >= 3)
	    repmc();
	else
	    anneal();
    }
    
    if (rfp)
	pclose(rfp);

    calcdist(bestchn, FALSE);
    sc_opt(dsqmat, bestchn, 3); /* Optimize side chains with DFIRE */

    /* Print final energy report */
    verboseflg = TRUE;

    printf("Final energy report:\n");
    e = e_model(dsqmat, bestchn, TRUE);

    puts("\nFINISHED!");

    /* write out final best structure */
    writepdb(bestchn, e);

    return 0;
}
