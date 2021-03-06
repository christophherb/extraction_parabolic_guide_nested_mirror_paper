/* Automatically generated file. Do not edit. 
 * Format:     ANSI C source code
 * Creator:    McStas <http://www.mcstas.org>
 * Instrument: /home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr (template_simple)
 * Date:       Wed Aug 18 10:01:52 2021
 * File:       ./ess_extraction.c
 * Compile:    cc -o template_simple.out ./ess_extraction.c 
 * CFLAGS=
 */


#define MCCODE_STRING "McStas 2.6.1 - May. 04, 2020"
#define FLAVOR "mcstas"
#define FLAVOR_UPPER "MCSTAS"
#define MC_USE_DEFAULT_MAIN
#define MC_TRACE_ENABLED
#define MC_EMBEDDED_RUNTIME

#line 1 "mccode-r.h"
/*******************************************************************************
*
* McCode, neutron/xray ray-tracing package
*         Copyright (C) 1997-2009, All rights reserved
*         Risoe National Laboratory, Roskilde, Denmark
*         Institut Laue Langevin, Grenoble, France
*
* Runtime: share/mccode-r.h
*
* %Identification
* Written by: KN
* Date:    Aug 29, 1997
* Release: McStas 2.6.1
* Version: $Revision$
*
* Runtime system header for McStas/McXtrace.
*
* In order to use this library as an external library, the following variables
* and macros must be declared (see details in the code)
*
*   struct mcinputtable_struct mcinputtable[];
*   int mcnumipar;
*   char mcinstrument_name[], mcinstrument_source[];
*   int mctraceenabled, mcdefaultmain;
*   extern MCNUM  mccomp_storein[];
*   extern MCNUM  mcAbsorbProp[];
*   extern MCNUM  mcScattered;
*   #define MCCODE_STRING "the McStas/McXtrace version"
*
* Usage: Automatically embbeded in the c code.
*
* $Id$
*
*******************************************************************************/

#ifndef MCCODE_R_H
#define MCCODE_R_H "$Revision$"

#include <math.h>
#include <string.h>
#include <strings.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <limits.h>
#include <errno.h>
#include <time.h>
#include <float.h>
#include <inttypes.h>

/* If the runtime is embedded in the simulation program, some definitions can
   be made static. */

#ifdef MC_EMBEDDED_RUNTIME
#define mcstatic static
#else
#define mcstatic
#endif

#ifdef __dest_os
#if (__dest_os == __mac_os)
#define MAC
#endif
#endif

#ifdef __FreeBSD__
#define NEED_STAT_H
#endif

#if defined(__APPLE__) && defined(__GNUC__)
#define NEED_STAT_H
#endif

#ifdef NEED_STAT_H
#include <sys/stat.h>
#endif

#ifndef MC_PATHSEP_C
#ifdef WIN32
#define MC_PATHSEP_C '\\'
#define MC_PATHSEP_S "\\"
#else  /* !WIN32 */
#define MC_PATHSEP_C '/'
#define MC_PATHSEP_S "/"
#endif /* !WIN32 */
#endif /* MC_PATHSEP_C */

#ifndef WIN32
#ifndef _POSIX_SOURCE
#define _POSIX_SOURCE 1
#endif
#endif

/* the version string is replaced when building distribution with mkdist */
#ifndef MCCODE_STRING
#define MCCODE_STRING "McStas 2.6.1 - May. 04, 2020"
#endif

#ifndef MCCODE_DATE
#define MCCODE_DATE "May. 04, 2020"
#endif

#ifndef MCCODE_VERSION
#define MCCODE_VERSION "2.6.1"
#endif

#ifndef MCCODE_NAME
#define MCCODE_NAME "McStas"
#endif

#ifndef MCCODE_PARTICLE
#define MCCODE_PARTICLE "neutron"
#endif

#ifndef MCCODE_LIBENV
#define MCCODE_LIBENV "MCSTAS"
#endif

#ifndef FLAVOR_UPPER
#define FLAVOR_UPPER MCCODE_NAME
#endif

#ifdef MC_PORTABLE
#ifndef NOSIGNALS
#define NOSIGNALS 1
#endif
#endif

#ifdef MAC
#ifndef NOSIGNALS
#define NOSIGNALS 1
#endif
#endif

#if (USE_MPI == 0)
#undef USE_MPI
#endif

#ifdef USE_MPI  /* default is to disable signals with MPI, as MPICH uses them to communicate */
#ifndef NOSIGNALS
#define NOSIGNALS 1
#endif
#endif

#if (NOSIGNALS == 0)
#undef NOSIGNALS
#endif

/* Note: the enum instr_formal_types definition MUST be kept
   synchronized with the one in mccode.h and with the
   instr_formal_type_names array in cogen.c. */
enum instr_formal_types
  {
    instr_type_double, instr_type_int, instr_type_string
  };
struct mcinputtable_struct { /* defines instrument parameters */
  char *name; /* name of parameter */
  void *par;  /* pointer to instrument parameter (variable) */
  enum instr_formal_types type;
  char *val;  /* default value */
};

typedef double MCNUM;
typedef struct {MCNUM x, y, z;} Coords;
typedef MCNUM Rotation[3][3];

/* the following variables are defined in the McStas generated C code
   but should be defined externally in case of independent library usage */
#ifndef DANSE
extern struct mcinputtable_struct mcinputtable[]; /* list of instrument parameters */
extern int    mcnumipar;                          /* number of instrument parameters */
extern char   mcinstrument_name[], mcinstrument_source[]; /* instrument name and filename */
extern char  *mcinstrument_exe;                           /* executable path = argv[0] or NULL */
extern MCNUM  mccomp_storein[]; /* 11 coords * number of components in instrument */
extern MCNUM  mcAbsorbProp[];
extern MCNUM  mcScattered;      /* number of SCATTER calls in current component */
extern MCNUM  mcRestore;        /* Flag to indicate if neutron needs to be restored */
#ifndef MC_ANCIENT_COMPATIBILITY
extern int mctraceenabled, mcdefaultmain;
#endif
#endif


/* Useful macros ============================================================ */

/* MPI stuff */

#ifdef USE_MPI
#include "mpi.h"

#ifdef OMPI_MPI_H  /* openmpi does not use signals: we may install our sighandler */
#undef NOSIGNALS
#endif

/*
 * MPI_MASTER(i):
 * execution of i only on master node
 */
#define MPI_MASTER(statement) { \
  if(mpi_node_rank == mpi_node_root)\
  { statement; } \
}

#ifndef MPI_REDUCE_BLOCKSIZE
#define MPI_REDUCE_BLOCKSIZE 1000
#endif

int mc_MPI_Sum(double* buf, long count);
int mc_MPI_Send(void *sbuf, long count, MPI_Datatype dtype, int dest);
int mc_MPI_Recv(void *rbuf, long count, MPI_Datatype dtype, int source);

/* MPI_Finalize exits gracefully and should be preferred to MPI_Abort */
#define exit(code) do {                                   \
    MPI_Finalize();                                       \
    exit(code);                                           \
  } while(0)

#else /* !USE_MPI */
#define MPI_MASTER(instr) instr
#endif /* USE_MPI */

#ifdef USE_MPI
static int mpi_node_count;
#endif

#ifdef USE_THREADS  /* user want threads */
#error Threading (USE_THREADS) support has been removed for very poor efficiency. Use MPI/SSH grid instead.
#endif


void   mcset_ncount(unsigned long long count);    /* wrapper to get mcncount */
unsigned long long int mcget_ncount(void);            /* wrapper to set mcncount */
unsigned long long mcget_run_num(void);           /* wrapper to get mcrun_num=0:mcncount */


/* Following part is only embedded when not redundant with mccode.h ========= */

#ifndef MCCODE_H

#ifndef NOSIGNALS
#include <signal.h>
#define SIG_MESSAGE(msg) strcpy(mcsig_message, msg);
#else
#define SIG_MESSAGE(msg)
#endif /* !NOSIGNALS */

/* Useful macros and constants ============================================== */

#ifndef FLT_MAX
#define FLT_MAX         3.40282347E+38F /* max decimal value of a "float" */
#endif

#ifndef MIN
#define MIN(a, b)  (((a) < (b)) ? (a) : (b))
#endif
#ifndef MAX
#define MAX(a, b)  (((a) > (b)) ? (a) : (b))
#endif
#ifndef SQR
#define SQR(x) ( (x) * (x) )
#endif
#ifndef SIGN
#define SIGN(x) (((x)>0.0)?(1):(-1))
#endif

#ifndef PI
# ifdef M_PI
#  define PI M_PI
# else
/* When using c99 in the CFLAGS, some of these consts
   are lost... Perhaps we should in fact include everything from
   https://www.gnu.org/software/libc/manual/html_node/Mathematical-Constants.html
*/
#  define PI 3.14159265358979323846
#  define M_PI PI
#  define M_PI_2 M_PI/2.0
#  define M_PI_4 M_PI/4.0
#  define M_1_PI 1.0/M_PI
#  define M_2_PI 2*M_1_PI
#  define M_2_SQRTPI 2/sqrt(M_PI)
#  define M_SQRT2 sqrt(2)
#  define M_SQRT1_2 sqrt(1/2)
# endif
#endif

#define RAD2MIN  ((180*60)/PI)
#define MIN2RAD  (PI/(180*60))
#define DEG2RAD  (PI/180)
#define RAD2DEG  (180/PI)
#define FWHM2RMS 0.424660900144    /* Convert between full-width-half-max and */
#define RMS2FWHM 2.35482004503     /* root-mean-square (standard deviation) */
#define HBAR     1.05457168e-34    /* [Js] h bar Planck constant CODATA 2002 */
#define MNEUTRON 1.67492728e-27    /* [kg] mass of neutron CODATA 2002 */
#define GRAVITY  9.81              /* [m/s^2] gravitational acceleration */
#define NA       6.02214179e23     /* [#atoms/g .mole] Avogadro's number*/


/* wrapper to get absolute and relative position of comp */
/* mccomp_posa and mccomp_posr are defined in McStas generated C code */
#define POS_A_COMP_INDEX(index) \
    (mccomp_posa[index])
#define POS_R_COMP_INDEX(index) \
    (mccomp_posr[index])
/* number of SCATTER calls in current comp: mcScattered defined in generated C code */
#define SCATTERED mcScattered
/* Flag to indicate if neutron needs to be restored: mcRestore defined in generated C code */
#define RESTORE mcRestore


/* Retrieve component information from the kernel */
/* Name, position and orientation (both absolute and relative)  */
/* Any component: For "redundancy", see comment by KN */
#define tmp_name_comp(comp) #comp
#define NAME_COMP(comp) tmp_name_comp(comp)
#define tmp_pos_a_comp(comp) (mcposa ## comp)
#define POS_A_COMP(comp) tmp_pos_a_comp(comp)
#define tmp_pos_r_comp(comp) (mcposr ## comp)
#define POS_R_COMP(comp) tmp_pos_r_comp(comp)
#define tmp_rot_a_comp(comp) (mcrota ## comp)
#define ROT_A_COMP(comp) tmp_rot_a_comp(comp)
#define tmp_rot_r_comp(comp) (mcrotr ## comp)
#define ROT_R_COMP(comp) tmp_rot_r_comp(comp)

/* Current component name, index, position and orientation */
#define NAME_CURRENT_COMP  NAME_COMP(mccompcurname)
#define INDEX_CURRENT_COMP mccompcurindex
#define POS_A_CURRENT_COMP POS_A_COMP(mccompcurname)
#define POS_R_CURRENT_COMP POS_R_COMP(mccompcurname)
#define ROT_A_CURRENT_COMP ROT_A_COMP(mccompcurname)
#define ROT_R_CURRENT_COMP ROT_R_COMP(mccompcurname)

/* Note: The two-stage approach to MC_GETPAR is NOT redundant; without it,
* after #define C sample, MC_GETPAR(C,x) would refer to component C, not to
* component sample. Such are the joys of ANSI C.

* Anyway the usage of MCGETPAR requires that we use sometimes bare names...
*/
#define MC_GETPAR2(comp, par) (mcc ## comp ## _ ## par)
#define MC_GETPAR(comp, par) MC_GETPAR2(comp,par)

/* MCDISPLAY/trace and debugging message sent to stdout */
#ifdef MC_TRACE_ENABLED
#define DEBUG
#endif

#ifdef DEBUG
#define mcDEBUG_INSTR() if(!mcdotrace); else { printf("\nINSTRUMENT:\n"); printf("Instrument '%s' (%s)\n", mcinstrument_name, mcinstrument_source); }
#define mcDEBUG_COMPONENT(name,c,t) if(!mcdotrace); else {\
  printf("COMPONENT: \"%s\"\n" \
         "POS: %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g\n", \
         name, c.x, c.y, c.z, t[0][0], t[0][1], t[0][2], \
         t[1][0], t[1][1], t[1][2], t[2][0], t[2][1], t[2][2]); \
  mcAccumulatedILength += coords_len(coords_sub(mcLastComp,c)); \
  printf("Component %30s AT (%g,%g,%g)    %g m from origin\n", name, c.x, c.y, c.z, mcAccumulatedILength); \
  mcLastComp=c;\
  }
#define mcDEBUG_INSTR_END() if(!mcdotrace); else printf("INSTRUMENT END:\n");
#define mcDEBUG_ENTER() if(!mcdotrace); else printf("ENTER:\n");
#define mcDEBUG_COMP(c) if(!mcdotrace); else printf("COMP: \"%s\"\n", c);
#define mcDEBUG_LEAVE() if(!mcdotrace); else printf("LEAVE:\n");
#define mcDEBUG_ABSORB() if(!mcdotrace); else printf("ABSORB:\n");
#else
#define mcDEBUG_INSTR()
#define mcDEBUG_COMPONENT(name,c,t)
#define mcDEBUG_INSTR_END()
#define mcDEBUG_ENTER()
#define mcDEBUG_COMP(c)
#define mcDEBUG_LEAVE()
#define mcDEBUG_ABSORB()
#endif

// mcDEBUG_STATE and mcDEBUG_SCATTER are defined by mcstas-r.h and mcxtrace-r.h



#ifdef TEST
#define test_printf printf
#else
#define test_printf while(0) printf
#endif

/* send MCDISPLAY message to stdout to show gemoetry */
void mcdis_magnify(char *what);
void mcdis_line(double x1, double y1, double z1,
                double x2, double y2, double z2);
void mcdis_dashed_linemcdis_dashed_line(double x1, double y1, double z1,
		       double x2, double y2, double z2, int n);
void mcdis_multiline(int count, ...);
void mcdis_rectangle(char* plane, double x, double y, double z,
		     double width, double height);
void mcdis_box(double x, double y, double z,
	       double width, double height, double length);
void mcdis_circle(char *plane, double x, double y, double z, double r);
void mcdis_Circle(double x, double y, double z, double r, double nx, double ny, double nz);
void mcdis_cylinder( double x, double y, double z,
        double r, double height, int N, double nx, double ny, double nz);
void mcdis_sphere(double x, double y, double z, double r, int N);

/* selection of random number generator. default is MT */
#ifndef MC_RAND_ALG
#define MC_RAND_ALG 1
#endif

#if MC_RAND_ALG == 0
   /* Use system random() (not recommended). */
#  define MC_RAND_MAX RAND_MAX
#elif MC_RAND_ALG == 1
   /* "Mersenne Twister", by Makoto Matsumoto and Takuji Nishimura. */
#  define MC_RAND_MAX ((unsigned long)0xffffffff)
#  define random mt_random
#  define srandom mt_srandom
#elif MC_RAND_ALG == 2
   /* Algorithm used in McStas CVS-080208 and earlier (not recommended). */
#  define MC_RAND_MAX 0x7fffffff
#  define random mc_random
#  define srandom mc_srandom
#else
#  error "Bad value for random number generator choice."
#endif

typedef int mc_int32_t;
mc_int32_t mc_random(void);
void mc_srandom (unsigned int x);
unsigned long mt_random(void);
void mt_srandom (unsigned long x);

double rand01();
double randpm1();
double rand0max(double max);
double randminmax(double min, double max);

double randnorm(void);
double randtriangle(void);

#ifndef DANSE
void mcinit(void);
void mcraytrace(void);
void mcsave(FILE *);
void mcfinally(void);
void mcdisplay(void);
#endif

/* simple vector algebra ==================================================== */
#define vec_prod(x, y, z, x1, y1, z1, x2, y2, z2) \
	vec_prod_func(&x, &y, &z, x1, y1, z1, x2, y2, z2)
mcstatic void vec_prod_func(double *x, double *y, double *z,
		double x1, double y1, double z1, double x2, double y2, double z2);

mcstatic double scalar_prod(
		double x1, double y1, double z1, double x2, double y2, double z2);

#define NORM(x,y,z) \
	norm_func(&x, &y, &z)
mcstatic void norm_func(double *x, double *y, double *z) {
	double temp = (*x * *x) + (*y * *y) + (*z * *z);
	if (temp != 0) {
		temp = sqrt(temp);
		*x /= temp;
		*y /= temp;
		*z /= temp;
	}
}
#define normal_vec(nx, ny, nz, x, y, z) \
    normal_vec_func(&(nx), &(ny), &(nz), x, y, z)
mcstatic void normal_vec_func(double *nx, double *ny, double *nz,
    double x, double y, double z);

/**
 * Rotate the vector vx,vy,vz psi radians around the vector ax,ay,az
 * and put the result in x,y,z.
 */
#define rotate(x, y, z, vx, vy, vz, phi, ax, ay, az) \
  do { \
    double mcrt_tmpx = (ax), mcrt_tmpy = (ay), mcrt_tmpz = (az); \
    double mcrt_vp, mcrt_vpx, mcrt_vpy, mcrt_vpz; \
    double mcrt_vnx, mcrt_vny, mcrt_vnz, mcrt_vn1x, mcrt_vn1y, mcrt_vn1z; \
    double mcrt_bx, mcrt_by, mcrt_bz; \
    double mcrt_cos, mcrt_sin; \
    NORM(mcrt_tmpx, mcrt_tmpy, mcrt_tmpz); \
    mcrt_vp = scalar_prod((vx), (vy), (vz), mcrt_tmpx, mcrt_tmpy, mcrt_tmpz); \
    mcrt_vpx = mcrt_vp*mcrt_tmpx; \
    mcrt_vpy = mcrt_vp*mcrt_tmpy; \
    mcrt_vpz = mcrt_vp*mcrt_tmpz; \
    mcrt_vnx = (vx) - mcrt_vpx; \
    mcrt_vny = (vy) - mcrt_vpy; \
    mcrt_vnz = (vz) - mcrt_vpz; \
    vec_prod(mcrt_bx, mcrt_by, mcrt_bz, \
             mcrt_tmpx, mcrt_tmpy, mcrt_tmpz, mcrt_vnx, mcrt_vny, mcrt_vnz); \
    mcrt_cos = cos((phi)); mcrt_sin = sin((phi)); \
    mcrt_vn1x = mcrt_vnx*mcrt_cos + mcrt_bx*mcrt_sin; \
    mcrt_vn1y = mcrt_vny*mcrt_cos + mcrt_by*mcrt_sin; \
    mcrt_vn1z = mcrt_vnz*mcrt_cos + mcrt_bz*mcrt_sin; \
    (x) = mcrt_vpx + mcrt_vn1x; \
    (y) = mcrt_vpy + mcrt_vn1y; \
    (z) = mcrt_vpz + mcrt_vn1z; \
  } while(0)

/**
 * Mirror (xyz) in the plane given by the point (rx,ry,rz) and normal (nx,ny,nz)
 *
 * TODO: This define is seemingly never used...
 */
#define mirror(x,y,z,rx,ry,rz,nx,ny,nz) \
  do { \
    double mcrt_tmpx= (nx), mcrt_tmpy = (ny), mcrt_tmpz = (nz); \
    double mcrt_tmpt; \
    NORM(mcrt_tmpx, mcrt_tmpy, mcrt_tmpz); \
    mcrt_tmpt=scalar_prod((rx),(ry),(rz),mcrt_tmpx,mcrt_tmpy,mcrt_tmpz); \
    (x) = rx -2 * mcrt_tmpt*mcrt_rmpx; \
    (y) = ry -2 * mcrt_tmpt*mcrt_rmpy; \
    (z) = rz -2 * mcrt_tmpt*mcrt_rmpz; \
  } while (0)

Coords coords_set(MCNUM x, MCNUM y, MCNUM z);
Coords coords_get(Coords a, MCNUM *x, MCNUM *y, MCNUM *z);
Coords coords_add(Coords a, Coords b);
Coords coords_sub(Coords a, Coords b);
Coords coords_neg(Coords a);
Coords coords_scale(Coords b, double scale);
double coords_sp(Coords a, Coords b);
Coords coords_xp(Coords b, Coords c);
double coords_len(Coords a);
void   coords_print(Coords a);
mcstatic void coords_norm(Coords* c);

void rot_set_rotation(Rotation t, double phx, double phy, double phz);
int  rot_test_identity(Rotation t);
void rot_mul(Rotation t1, Rotation t2, Rotation t3);
void rot_copy(Rotation dest, Rotation src);
void rot_transpose(Rotation src, Rotation dst);
Coords rot_apply(Rotation t, Coords a);

void mccoordschange(Coords a, Rotation t, double *x, double *y, double *z,
    double *vx, double *vy, double *vz, double *sx, double *sy, double *sz);
void
mccoordschange_polarisation(Rotation t, double *sx, double *sy, double *sz);

double mcestimate_error(double N, double p1, double p2);
void mcreadparams(void);

/* this is now in mcstas-r.h and mcxtrace-r.h as the number of state parameters is no longer equal*/
/* void mcsetstate(double x, double y, double z, double vx, double vy, double vz,
                double t, double sx, double sy, double sz, double p);
*/
void mcgenstate(void);

/* trajectory/shape intersection routines */
int inside_rectangle(double, double, double, double);
int box_intersect(double *dt_in, double *dt_out, double x, double y, double z,
    double vx, double vy, double vz, double dx, double dy, double dz);
int cylinder_intersect(double *t0, double *t1, double x, double y, double z,
    double vx, double vy, double vz, double r, double h);
int sphere_intersect(double *t0, double *t1, double x, double y, double z,
                 double vx, double vy, double vz, double r);
/* second order equation roots */
int solve_2nd_order(double *t1, double *t2,
    double A,  double B,  double C);

/* random vector generation to shape */
void randvec_target_circle(double *xo, double *yo, double *zo,
    double *solid_angle, double xi, double yi, double zi, double radius);
#define randvec_target_sphere randvec_target_circle
void randvec_target_rect_angular(double *xo, double *yo, double *zo,
    double *solid_angle,
               double xi, double yi, double zi, double height, double width, Rotation A);
#define randvec_target_rect(p0,p1,p2,p3,p4,p5,p6,p7,p8,p9)  randvec_target_rect_real(p0,p1,p2,p3,p4,p5,p6,p7,p8,p9,0,0,0,1)
void randvec_target_rect_real(double *xo, double *yo, double *zo,
    double *solid_angle,
	       double xi, double yi, double zi, double height, double width, Rotation A,
			 double lx, double ly, double lz, int order);

/* this is the main() */
int mccode_main(int argc, char *argv[]);


#endif /* !MCCODE_H */

#ifndef MCCODE_R_IO_H
#define MCCODE_R_IO_H "$Revision$"

#if (USE_NEXUS == 0)
#undef USE_NEXUS
#endif

#ifndef CHAR_BUF_LENGTH
#define CHAR_BUF_LENGTH 1024
#endif

/* I/O section part ========================================================= */

/* ========================================================================== */

/*                               MCCODE_R_IO_C                                */

/* ========================================================================== */


/* main DETECTOR structure which stores most information to write to data files */
struct mcdetector_struct {
  char   filename[CHAR_BUF_LENGTH];   /* file name of monitor */
  char   position[CHAR_BUF_LENGTH];   /* position of detector component */
  char   component[CHAR_BUF_LENGTH];  /* component instance name */
  char   instrument[CHAR_BUF_LENGTH]; /* instrument name */
  char   type[CHAR_BUF_LENGTH];       /* data type, e.g. 0d, 1d, 2d, 3d */
  char   user[CHAR_BUF_LENGTH];       /* user name, e.g. HOME */
  char   date[CHAR_BUF_LENGTH];       /* date of simulation end/write time */
  char   title[CHAR_BUF_LENGTH];      /* title of detector */
  char   xlabel[CHAR_BUF_LENGTH];     /* X axis label */
  char   ylabel[CHAR_BUF_LENGTH];     /* Y axis label */
  char   zlabel[CHAR_BUF_LENGTH];     /* Z axis label */
  char   xvar[CHAR_BUF_LENGTH];       /* X variable name */
  char   yvar[CHAR_BUF_LENGTH];       /* Y variable name */
  char   zvar[CHAR_BUF_LENGTH];       /* Z variable name */
  char   ncount[CHAR_BUF_LENGTH];     /* number of events initially generated */
  char   limits[CHAR_BUF_LENGTH];     /* X Y Z limits, e.g. [xmin xmax ymin ymax zmin zmax] */
  char   variables[CHAR_BUF_LENGTH];  /* variables written into data block */
  char   statistics[CHAR_BUF_LENGTH]; /* center, mean and half width along axis */
  char   signal[CHAR_BUF_LENGTH];     /* min max and mean of signal (data block) */
  char   values[CHAR_BUF_LENGTH];     /* integrated values e.g. [I I_err N] */
  double xmin,xmax;                   /* min max of axes */
  double ymin,ymax;
  double zmin,zmax;
  double intensity;                   /* integrated values for data block */
  double error;
  double events;
  double min;                         /* statistics for data block */
  double max;
  double mean;
  double centerX;                     /* statistics for axes */
  double halfwidthX;
  double centerY;
  double halfwidthY;
  int    rank;                        /* dimensionaly of monitor, e.g. 0 1 2 3 */
  char   istransposed;                /* flag to transpose matrix for some formats */

  long   m,n,p;                       /* dimensions of data block and along axes */
  long   date_l;                      /* same as date, but in sec since 1970 */

  double *p0, *p1, *p2;               /* pointers to saved data, NULL when freed */
  char   format[CHAR_BUF_LENGTH];    /* format for file generation */
};

typedef struct mcdetector_struct MCDETECTOR;

static   char *mcdirname             = NULL;      /* name of output directory */
static   char *mcsiminfo_name        = "mccode";  /* default output sim file name */
char    *mcformat                    = NULL;      /* NULL (default) or a specific format */

/* file I/O definitions and function prototypes */

#ifndef MC_EMBEDDED_RUNTIME /* the mcstatic variables (from mccode-r.c) */
extern FILE * mcsiminfo_file;     /* handle to the output siminfo file */
extern int    mcgravitation;      /* flag to enable gravitation */
extern int    mcdotrace;          /* flag to print MCDISPLAY messages */
#else
mcstatic FILE *mcsiminfo_file        = NULL;
#endif

/* I/O function prototypes ================================================== */

/* output functions */
MCDETECTOR mcdetector_out_0D(char *t, double p0, double p1, double p2, char *c, Coords pos);
MCDETECTOR mcdetector_out_1D(char *t, char *xl, char *yl,
                  char *xvar, double x1, double x2, long n,
                  double *p0, double *p1, double *p2, char *f, char *c, Coords pos);
MCDETECTOR mcdetector_out_2D(char *t, char *xl, char *yl,
                  double x1, double x2, double y1, double y2, long m,
                  long n, double *p0, double *p1, double *p2, char *f,
                  char *c, Coords pos);
MCDETECTOR mcdetector_out_list(char *t, char *xl, char *yl,
                  long m, long n,
                  double *p1, char *f,
                  char *c, Coords posa);

/* wrappers to output functions, that automatically set NAME and POSITION */
#define DETECTOR_OUT(p0,p1,p2) mcdetector_out_0D(NAME_CURRENT_COMP,p0,p1,p2,NAME_CURRENT_COMP,POS_A_CURRENT_COMP)
#define DETECTOR_OUT_0D(t,p0,p1,p2) mcdetector_out_0D(t,p0,p1,p2,NAME_CURRENT_COMP,POS_A_CURRENT_COMP)
#define DETECTOR_OUT_1D(t,xl,yl,xvar,x1,x2,n,p0,p1,p2,f) \
     mcdetector_out_1D(t,xl,yl,xvar,x1,x2,n,p0,p1,p2,f,NAME_CURRENT_COMP,POS_A_CURRENT_COMP)
#define DETECTOR_OUT_2D(t,xl,yl,x1,x2,y1,y2,m,n,p0,p1,p2,f) \
     mcdetector_out_2D(t,xl,yl,x1,x2,y1,y2,m,n,p0,p1,p2,f,NAME_CURRENT_COMP,POS_A_CURRENT_COMP)

#ifdef USE_NEXUS
#include "napi.h"
NXhandle nxhandle;
#endif

#endif /* ndef MCCODE_R_IO_H */

#endif /* MCCODE_R_H */
/* End of file "mccode-r.h". */

#line 712 "./ess_extraction.c"

#line 1 "mcstas-r.h"
/*******************************************************************************
*
* McStas, neutron ray-tracing package
*         Copyright (C) 1997-2009, All rights reserved
*         Risoe National Laboratory, Roskilde, Denmark
*         Institut Laue Langevin, Grenoble, France
*
* Runtime: share/mcstas-r.h
*
* %Identification
* Written by: KN
* Date:    Aug 29, 1997
* Release: McStas X.Y
* Version: $Revision$
*
* Runtime system header for McStas.
*
* In order to use this library as an external library, the following variables
* and macros must be declared (see details in the code)
*
*   struct mcinputtable_struct mcinputtable[];
*   int mcnumipar;
*   char mcinstrument_name[], mcinstrument_source[];
*   int mctraceenabled, mcdefaultmain;
*   extern MCNUM  mccomp_storein[];
*   extern MCNUM  mcAbsorbProp[];
*   extern MCNUM  mcScattered;
*   #define MCCODE_STRING "the McStas version"
*
* Usage: Automatically embbeded in the c code.
*
* $Id$
*
*******************************************************************************/

#ifndef MCSTAS_R_H
#define MCSTAS_R_H "$Revision$"

/* Following part is only embedded when not redundent with mcstas.h ========= */

#ifndef MCCODE_H

#define AA2MS    629.622368        /* Convert k[1/AA] to v[m/s] */
#define MS2AA    1.58825361e-3     /* Convert v[m/s] to k[1/AA] */
#define K2V      AA2MS
#define V2K      MS2AA
#define Q2V      AA2MS
#define V2Q      MS2AA
#define SE2V     437.393377        /* Convert sqrt(E)[meV] to v[m/s] */
#define VS2E     5.22703725e-6     /* Convert (v[m/s])**2 to E[meV] */

#define SCATTER do {mcDEBUG_SCATTER(mcnlx, mcnly, mcnlz, mcnlvx, mcnlvy, mcnlvz, \
        mcnlt,mcnlsx,mcnlsy,mcnlsz, mcnlp); mcScattered++;} while(0)
#define ABSORB do {mcDEBUG_STATE(mcnlx, mcnly, mcnlz, mcnlvx, mcnlvy, mcnlvz, \
        mcnlt,mcnlsx,mcnlsy,mcnlsz, mcnlp); mcDEBUG_ABSORB(); MAGNET_OFF; goto mcabsorb;} while(0)

#define STORE_NEUTRON(index, x, y, z, vx, vy, vz, t, sx, sy, sz, p) \
  mcstore_neutron(mccomp_storein,index, x, y, z, vx, vy, vz, t, sx, sy, sz, p);
#define RESTORE_NEUTRON(index, x, y, z, vx, vy, vz, t, sx, sy, sz, p) \
  mcrestore_neutron(mccomp_storein,index, &x, &y, &z, &vx, &vy, &vz, &t, &sx, &sy, &sz, &p);

#define MAGNET_ON \
  do { \
    mcMagnet = 1; \
  } while(0)

#define MAGNET_OFF \
  do { \
    mcMagnet = 0; \
  } while(0)

#define ALLOW_BACKPROP \
  do { \
    mcallowbackprop = 1; \
  } while(0)

#define DISALLOW_BACKPROP \
  do { \
    mcallowbackprop = 0; \
  } while(0)

#define PROP_MAGNET(dt) \
  do { \
  }while (0)
    /* change coordinates from local system to magnet system */
/*    Rotation rotLM, rotTemp; \
      Coords   posLM = coords_sub(POS_A_CURRENT_COMP, mcMagnetPos); \
      rot_transpose(ROT_A_CURRENT_COMP, rotTemp); \
      rot_mul(rotTemp, mcMagnetRot, rotLM); \
      mcMagnetPrecession(mcnlx, mcnly, mcnlz, mcnlt, mcnlvx, mcnlvy, mcnlvz, \
               &mcnlsx, &mcnlsy, &mcnlsz, dt, posLM, rotLM); \
      } while(0)
*/

#define mcPROP_DT(dt) \
  do { \
    if (mcMagnet && dt > 0) PROP_MAGNET(dt);\
    mcnlx += mcnlvx*(dt); \
    mcnly += mcnlvy*(dt); \
    mcnlz += mcnlvz*(dt); \
    mcnlt += (dt); \
    if (isnan(p) || isinf(p)) { mcAbsorbProp[INDEX_CURRENT_COMP]++; ABSORB; }\
  } while(0)

/* ADD: E. Farhi, Aug 6th, 2001 PROP_GRAV_DT propagation with acceleration */
#define PROP_GRAV_DT(dt, Ax, Ay, Az) \
  do { \
    if(dt < 0 && mcallowbackprop == 0) { mcAbsorbProp[INDEX_CURRENT_COMP]++; ABSORB; }\
    if (mcMagnet) printf("Spin precession gravity\n"); \
    mcnlx  += mcnlvx*(dt) + (Ax)*(dt)*(dt)/2; \
    mcnly  += mcnlvy*(dt) + (Ay)*(dt)*(dt)/2; \
    mcnlz  += mcnlvz*(dt) + (Az)*(dt)*(dt)/2; \
    mcnlvx += (Ax)*(dt); \
    mcnlvy += (Ay)*(dt); \
    mcnlvz += (Az)*(dt); \
    mcnlt  += (dt); \
    DISALLOW_BACKPROP;\
  } while(0)


#define PROP_DT(dt) \
  do { \
    if(dt < 0) { RESTORE=1; goto mcabsorbComp; }; \
    if (mcgravitation) { Coords mcLocG; double mc_gx, mc_gy, mc_gz; \
    mcLocG = rot_apply(ROT_A_CURRENT_COMP, coords_set(0,-GRAVITY,0)); \
    coords_get(mcLocG, &mc_gx, &mc_gy, &mc_gz); \
    PROP_GRAV_DT(dt, mc_gx, mc_gy, mc_gz); } \
    else mcPROP_DT(dt); \
    DISALLOW_BACKPROP;\
  } while(0)


#define PROP_Z0 \
  do { \
    if (mcgravitation) { Coords mcLocG; int mc_ret; \
    double mc_dt, mc_gx, mc_gy, mc_gz; \
    mcLocG = rot_apply(ROT_A_CURRENT_COMP, coords_set(0,-GRAVITY,0)); \
    coords_get(mcLocG, &mc_gx, &mc_gy, &mc_gz); \
    mc_ret = solve_2nd_order(&mc_dt, NULL, -mc_gz/2, -mcnlvz, -mcnlz); \
    if (mc_ret && mc_dt>=0) {PROP_GRAV_DT(mc_dt, mc_gx, mc_gy, mc_gz); mcnlz=0;}\
    else { if (mcallowbackprop ==0) {mcAbsorbProp[INDEX_CURRENT_COMP]++; ABSORB; }}; }\
    else mcPROP_Z0; \
    DISALLOW_BACKPROP;\
  } while(0)

#define mcPROP_Z0 \
  do { \
    double mc_dt; \
    if(mcnlvz == 0) { mcAbsorbProp[INDEX_CURRENT_COMP]++; ABSORB; }; \
    mc_dt = -mcnlz/mcnlvz; \
    if(mc_dt < 0 && mcallowbackprop == 0) { mcAbsorbProp[INDEX_CURRENT_COMP]++; ABSORB; }; \
    mcPROP_DT(mc_dt); \
    mcnlz = 0; \
    DISALLOW_BACKPROP;\
  } while(0)

#define PROP_X0 \
  do { \
    if (mcgravitation) { Coords mcLocG; int mc_ret; \
    double mc_dt, mc_gx, mc_gy, mc_gz; \
    mcLocG = rot_apply(ROT_A_CURRENT_COMP, coords_set(0,-GRAVITY,0)); \
    coords_get(mcLocG, &mc_gx, &mc_gy, &mc_gz); \
    mc_ret = solve_2nd_order(&mc_dt, NULL, -mc_gx/2, -mcnlvx, -mcnlx); \
    if (mc_ret && mc_dt>=0) PROP_GRAV_DT(mc_dt, mc_gx, mc_gy, mc_gz); \
    else { if (mcallowbackprop ==0) {mcAbsorbProp[INDEX_CURRENT_COMP]++; ABSORB; }}; }\
    else mcPROP_X0; \
    DISALLOW_BACKPROP;\
  } while(0)

#define mcPROP_X0 \
  do { \
    double mc_dt; \
    if(mcnlvx == 0) { mcAbsorbProp[INDEX_CURRENT_COMP]++; ABSORB; }; \
    mc_dt = -mcnlx/mcnlvx; \
    if(mc_dt < 0 && mcallowbackprop == 0) { mcAbsorbProp[INDEX_CURRENT_COMP]++; ABSORB; }; \
    mcPROP_DT(mc_dt); \
    mcnlx = 0; \
    DISALLOW_BACKPROP;\
  } while(0)

#define PROP_Y0 \
  do { \
    if (mcgravitation) { Coords mcLocG; int mc_ret; \
    double mc_dt, mc_gx, mc_gy, mc_gz; \
    mcLocG = rot_apply(ROT_A_CURRENT_COMP, coords_set(0,-GRAVITY,0)); \
    coords_get(mcLocG, &mc_gx, &mc_gy, &mc_gz); \
    mc_ret = solve_2nd_order(&mc_dt, NULL, -mc_gy/2, -mcnlvy, -mcnly); \
    if (mc_ret && mc_dt>=0) PROP_GRAV_DT(mc_dt, mc_gx, mc_gy, mc_gz); \
    else { if (mcallowbackprop ==0) {mcAbsorbProp[INDEX_CURRENT_COMP]++; ABSORB; }}; }\
    else mcPROP_Y0; \
    DISALLOW_BACKPROP;\
  } while(0)


#define mcPROP_Y0 \
  do { \
    double mc_dt; \
    if(mcnlvy == 0) { mcAbsorbProp[INDEX_CURRENT_COMP]++; ABSORB; }; \
    mc_dt = -mcnly/mcnlvy; \
    if(mc_dt < 0 && mcallowbackprop == 0) { mcAbsorbProp[INDEX_CURRENT_COMP]++; ABSORB; }; \
    mcPROP_DT(mc_dt); \
    mcnly = 0; \
    DISALLOW_BACKPROP; \
  } while(0)

/*moved from mccode-r.h*/
void mcsetstate(double x, double y, double z, double vx, double vy, double vz,
                double t, double sx, double sy, double sz, double p);

#ifdef DEBUG

#define mcDEBUG_STATE(x,y,z,vx,vy,vz,t,sx,sy,sz,p) if(!mcdotrace); else \
  printf("STATE: %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g\n", \
         x,y,z,vx,vy,vz,t,sx,sy,sz,p);
#define mcDEBUG_SCATTER(x,y,z,vx,vy,vz,t,sx,sy,sz,p) if(!mcdotrace); else \
  printf("SCATTER: %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g\n", \
         x,y,z,vx,vy,vz,t,sx,sy,sz,p);

#else

#define mcDEBUG_STATE(x,y,z,vx,vy,vz,t,sx,sy,sz,p)
#define mcDEBUG_SCATTER(x,y,z,vx,vy,vz,t,sx,sy,sz,p)

#endif

#endif /* !MCCODE_H */

#endif /* MCSTAS_R_H */
/* End of file "mcstas-r.h". */

#line 945 "./ess_extraction.c"

#line 1 "mccode-r.c"
/*******************************************************************************
*
* McCode, neutron/xray ray-tracing package
*         Copyright (C) 1997-2009, All rights reserved
*         Risoe National Laboratory, Roskilde, Denmark
*         Institut Laue Langevin, Grenoble, France
*
* Runtime: share/mccode-r.c
*
* %Identification
* Written by: KN
* Date:    Aug 29, 1997
* Release: McStas X.Y/McXtrace X.Y
* Version: $Revision$
*
* Runtime system for McStas and McXtrace.
* Embedded within instrument in runtime mode.
* Contains SECTIONS:
*   MPI handling (sum, send, recv)
*   format definitions
*   I/O
*   mcdisplay support
*   random numbers
*   coordinates handling
*   vectors math (solve 2nd order, normals, randvec...)
*   parameter handling
*   signal and main handlers
*
* Usage: Automatically embbeded in the c code whenever required.
*
* $Id$
*
*******************************************************************************/

/*******************************************************************************
* The I/O format definitions and functions
*******************************************************************************/


/** Include header files to avoid implicit declarations (not allowed on LLVM) */
#include <ctype.h>
#include <sys/types.h>

// UNIX specific headers (non-Windows)
#if defined(__unix__) || defined(__APPLE__)
#include <unistd.h>
#endif

#include <sys/stat.h>

#ifdef _WIN32 
#include <direct.h>
# define  mkdir( D, M )   _mkdir( D ) 
#endif 

#ifndef DANSE
#ifdef MC_ANCIENT_COMPATIBILITY
int mctraceenabled = 0;
int mcdefaultmain  = 0;
#endif
/* else defined directly in the McCode generated C code */

static   long mcseed                 = 0; /* seed for random generator */
static   long mcstartdate            = 0; /* start simulation time */
static   int  mcdisable_output_files = 0; /* --no-output-files */
mcstatic int  mcgravitation          = 0; /* use gravitation flag, for PROP macros */
int      mcMagnet                    = 0; /* magnet stack flag */
mcstatic int  mcdotrace              = 0; /* flag for --trace and messages for DISPLAY */
int      mcallowbackprop             = 0;         /* flag to enable negative/backprop */

/* Number of particle histories to simulate. */
#ifdef NEUTRONICS
mcstatic unsigned long long int mcncount             = 1;
mcstatic unsigned long long int mcrun_num            = 0;
#else
mcstatic unsigned long long int mcncount             = 1000000;
mcstatic unsigned long long int mcrun_num            = 0;
#endif /* NEUTRONICS */

#else
#include "mcstas-globals.h"
#endif /* !DANSE */

/* SECTION: MPI handling ==================================================== */

#ifdef USE_MPI
/* MPI rank */
static int mpi_node_rank;
static int mpi_node_root = 0;


/*******************************************************************************
* mc_MPI_Reduce: Gathers arrays from MPI nodes using Reduce function.
*******************************************************************************/
int mc_MPI_Sum(double *sbuf, long count)
{
  if (!sbuf || count <= 0) return(MPI_SUCCESS); /* nothing to reduce */
  else {
    /* we must cut the buffer into blocks not exceeding the MPI max buffer size of 32000 */
    long   offset=0;
    double *rbuf=NULL;
    int    length=MPI_REDUCE_BLOCKSIZE; /* defined in mccode-r.h */
    int    i=0;
    rbuf = calloc(count, sizeof(double));
    if (!rbuf)
      exit(-fprintf(stderr, "Error: Out of memory %li (mc_MPI_Sum)\n", count*sizeof(double)));
    while (offset < count) {
      if (!length || offset+length > count-1) length=count-offset;
      else length=MPI_REDUCE_BLOCKSIZE;
      if (MPI_Reduce((double*)(sbuf+offset), (double*)(rbuf+offset),
              length, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD) != MPI_SUCCESS)
        return MPI_ERR_COUNT;
      offset += length;
    }

    for (i=0; i<count; i++) sbuf[i] = rbuf[i];
    free(rbuf);
  }
  return MPI_SUCCESS;
} /* mc_MPI_Sum */

/*******************************************************************************
* mc_MPI_Send: Send array to MPI node by blocks to avoid buffer limit
*******************************************************************************/
int mc_MPI_Send(void *sbuf,
                  long count, MPI_Datatype dtype,
                  int dest)
{
  int dsize;
  long offset=0;
  int  tag=1;
  int  length=MPI_REDUCE_BLOCKSIZE; /* defined in mccode-r.h */

  if (!sbuf || count <= 0) return(MPI_SUCCESS); /* nothing to send */
  MPI_Type_size(dtype, &dsize);

  while (offset < count) {
    if (offset+length > count-1) length=count-offset;
    else length=MPI_REDUCE_BLOCKSIZE;
    if (MPI_Send((void*)(sbuf+offset*dsize), length, dtype, dest, tag++, MPI_COMM_WORLD) != MPI_SUCCESS)
      return MPI_ERR_COUNT;
    offset += length;
  }

  return MPI_SUCCESS;
} /* mc_MPI_Send */

/*******************************************************************************
* mc_MPI_Recv: Receives arrays from MPI nodes by blocks to avoid buffer limit
*             the buffer must have been allocated previously.
*******************************************************************************/
int mc_MPI_Recv(void *sbuf,
                  long count, MPI_Datatype dtype,
                  int source)
{
  int dsize;
  long offset=0;
  int  tag=1;
  int  length=MPI_REDUCE_BLOCKSIZE; /* defined in mccode-r.h */

  if (!sbuf || count <= 0) return(MPI_SUCCESS); /* nothing to recv */
  MPI_Type_size(dtype, &dsize);

  while (offset < count) {
    if (offset+length > count-1) length=count-offset;
    else length=MPI_REDUCE_BLOCKSIZE;
    if (MPI_Recv((void*)(sbuf+offset*dsize), length, dtype, source, tag++,
            MPI_COMM_WORLD, MPI_STATUS_IGNORE) != MPI_SUCCESS)
      return MPI_ERR_COUNT;
    offset += length;
  }

  return MPI_SUCCESS;
} /* mc_MPI_Recv */

#endif /* USE_MPI */

/* SECTION: parameters handling ============================================= */

/* Instrument input parameter type handling. */
/*******************************************************************************
* mcparm_double: extract double value from 's' into 'vptr'
*******************************************************************************/
static int
mcparm_double(char *s, void *vptr)
{
  char *p;
  double *v = (double *)vptr;

  if (!s) { *v = 0; return(1); }
  *v = strtod(s, &p);
  if(*s == '\0' || (p != NULL && *p != '\0') || errno == ERANGE)
    return 0;                        /* Failed */
  else
    return 1;                        /* Success */
}

/*******************************************************************************
* mcparminfo_double: display parameter type double
*******************************************************************************/
static char *
mcparminfo_double(char *parmname)
{
  return "double";
}

/*******************************************************************************
* mcparmerror_double: display error message when failed extract double
*******************************************************************************/
static void
mcparmerror_double(char *parm, char *val)
{
  fprintf(stderr, "Error: Invalid value '%s' for floating point parameter %s (mcparmerror_double)\n",
          val, parm);
}

/*******************************************************************************
* mcparmprinter_double: convert double to string
*******************************************************************************/
static void
mcparmprinter_double(char *f, void *vptr)
{
  double *v = (double *)vptr;
  sprintf(f, "%g", *v);
}

/*******************************************************************************
* mcparm_int: extract int value from 's' into 'vptr'
*******************************************************************************/
static int
mcparm_int(char *s, void *vptr)
{
  char *p;
  int *v = (int *)vptr;
  long x;

  if (!s) { *v = 0; return(1); }
  *v = 0;
  x = strtol(s, &p, 10);
  if(x < INT_MIN || x > INT_MAX)
    return 0;                        /* Under/overflow */
  *v = x;
  if(*s == '\0' || (p != NULL && *p != '\0') || errno == ERANGE)
    return 0;                        /* Failed */
  else
    return 1;                        /* Success */
}

/*******************************************************************************
* mcparminfo_int: display parameter type int
*******************************************************************************/
static char *
mcparminfo_int(char *parmname)
{
  return "int";
}

/*******************************************************************************
* mcparmerror_int: display error message when failed extract int
*******************************************************************************/
static void
mcparmerror_int(char *parm, char *val)
{
  fprintf(stderr, "Error: Invalid value '%s' for integer parameter %s (mcparmerror_int)\n",
          val, parm);
}

/*******************************************************************************
* mcparmprinter_int: convert int to string
*******************************************************************************/
static void
mcparmprinter_int(char *f, void *vptr)
{
  int *v = (int *)vptr;
  sprintf(f, "%d", *v);
}

/*******************************************************************************
* mcparm_string: extract char* value from 's' into 'vptr' (copy)
*******************************************************************************/
static int
mcparm_string(char *s, void *vptr)
{
  char **v = (char **)vptr;
  if (!s) { *v = NULL; return(1); }
  *v = (char *)malloc(strlen(s) + 1);
  if(*v == NULL)
  {
    exit(-fprintf(stderr, "Error: Out of memory %li (mcparm_string).\n", (long)strlen(s) + 1));
  }
  strcpy(*v, s);
  return 1;                        /* Success */
}

/*******************************************************************************
* mcparminfo_string: display parameter type string
*******************************************************************************/
static char *
mcparminfo_string(char *parmname)
{
  return "string";
}

/*******************************************************************************
* mcparmerror_string: display error message when failed extract string
*******************************************************************************/
static void
mcparmerror_string(char *parm, char *val)
{
  fprintf(stderr, "Error: Invalid value '%s' for string parameter %s (mcparmerror_string)\n",
          val, parm);
}

/*******************************************************************************
* mcparmprinter_string: convert string to string (including esc chars)
*******************************************************************************/
static void
mcparmprinter_string(char *f, void *vptr)
{
  char **v = (char **)vptr;
  char *p;

  if (!*v) { *f='\0'; return; }
  strcpy(f, "");
  for(p = *v; *p != '\0'; p++)
  {
    switch(*p)
    {
      case '\n':
        strcat(f, "\\n");
        break;
      case '\r':
        strcat(f, "\\r");
        break;
      case '"':
        strcat(f, "\\\"");
        break;
      case '\\':
        strcat(f, "\\\\");
        break;
      default:
        strncat(f, p, 1);
    }
  }
  /* strcat(f, "\""); */
} /* mcparmprinter_string */

/* now we may define the parameter structure, using previous functions */
static struct
  {
    int (*getparm)(char *, void *);
    char * (*parminfo)(char *);
    void (*error)(char *, char *);
    void (*printer)(char *, void *);
} mcinputtypes[] = {
  {
    mcparm_double, mcparminfo_double, mcparmerror_double,
    mcparmprinter_double
  }, {
    mcparm_int, mcparminfo_int, mcparmerror_int,
    mcparmprinter_int
  }, {
    mcparm_string, mcparminfo_string, mcparmerror_string,
    mcparmprinter_string
  }
};

/*******************************************************************************
* mcestimate_error: compute sigma from N,p,p2 in Gaussian large numbers approx
*******************************************************************************/
double mcestimate_error(double N, double p1, double p2)
{
  double pmean, n1;
  if(N <= 1)
    return p1;
  pmean = p1 / N;
  n1 = N - 1;
  /* Note: underflow may cause p2 to become zero; the fabs() below guards
     against this. */
  return sqrt((N/n1)*fabs(p2 - pmean*pmean));
}

double (*mcestimate_error_p)
  (double V2, double psum, double p2sum)=mcestimate_error;

/* ========================================================================== */

/*                               MCCODE_R_IO_C                                */

/* ========================================================================== */

#ifndef MCCODE_R_IO_C
#define MCCODE_R_IO_C "$Revision$"

/* SECTION: file i/o handling ================================================ */

#ifndef HAVE_STRCASESTR
// from msysgit: https://code.google.com/p/msysgit/source/browse/compat/strcasestr.c
char *strcasestr(const char *haystack, const char *needle)
{
  int nlen = strlen(needle);
  int hlen = strlen(haystack) - nlen + 1;
  int i;

  for (i = 0; i < hlen; i++) {
    int j;
    for (j = 0; j < nlen; j++) {
            unsigned char c1 = haystack[i+j];
            unsigned char c2 = needle[j];
            if (toupper(c1) != toupper(c2))
                    goto next;
    }
    return (char *) haystack + i;
  next:
    ;
  }
  return NULL;
}


#endif
#ifndef HAVE_STRCASECMP
int strcasecmp( const char *s1, const char *s2 )
{
  int c1, c2;
  do {
    c1 = tolower( (unsigned char) *s1++ );
    c2 = tolower( (unsigned char) *s2++ );
  } while (c1 == c2 && c1 != 0);
  return c2 > c1 ? -1 : c1 > c2;
}
#endif

/*******************************************************************************
* mcfull_file: allocates a full file name=mcdirname+file. Catenate extension if missing.
*******************************************************************************/
char *mcfull_file(char *name, char *ext)
{
  int   dirlen=0;
  char *mem   =NULL;

  dirlen = mcdirname ? strlen(mcdirname) : 0;
  mem = (char*)malloc(dirlen + strlen(name) + CHAR_BUF_LENGTH);
  if(!mem) {
    exit(-fprintf(stderr, "Error: Out of memory %li (mcfull_file)\n", (long)(dirlen + strlen(name) + 256)));
  }
  strcpy(mem, "");

  /* prepend directory name to path if name does not contain a path */
  if (dirlen > 0 && !strchr(name, MC_PATHSEP_C)) {
    strcat(mem, mcdirname);
    strcat(mem, MC_PATHSEP_S);
  } /* dirlen */

  strcat(mem, name);
  if (!strchr(name, '.') && ext && strlen(ext))
  { /* add extension if not in file name already */
    strcat(mem, ".");
    strcat(mem, ext);
  }
  return(mem);
} /* mcfull_file */

/*******************************************************************************
* mcnew_file: opens a new file within mcdirname if non NULL
*             the file is opened in "a" (append, create if does not exist)
*             the extension 'ext' is added if the file name does not include one.
*             the last argument is set to 0 if file did not exist, else to 1.
*******************************************************************************/
FILE *mcnew_file(char *name, char *ext, int *exists)
{
  char *mem;
  FILE *file=NULL;

  if (!name || strlen(name) == 0 || mcdisable_output_files) return(NULL);
  
  mem  = mcfull_file(name, ext); /* create mcdirname/name.ext */
  
  /* check for existence */
  file = fopen(mem, "r"); /* for reading -> fails if does not exist */
  if (file) {
    fclose(file);
    *exists=1;
  } else
    *exists=0;
  
  /* open the file for writing/appending */
#ifdef USE_NEXUS
  if (mcformat && strcasestr(mcformat, "NeXus")) {
    /* NXhandle nxhandle is defined in the .h with USE_NEXUS */
    NXaccess mode = (*exists ? NXACC_CREATE5 | NXACC_RDWR : NXACC_CREATE5);
      
    if (NXopen(mem, mode, &nxhandle) != NX_OK)
      file = NULL;
    else
      file = (FILE*)&nxhandle; /* to make it non NULL */
  } else
#endif
    file = fopen(mem, "a+"); 
    
  if(!file)
    fprintf(stderr, "Warning: could not open output file '%s' for %s (mcnew_file)\n", 
      mem, *exists ? "append" : "create");
  free(mem);

  return file;
} /* mcnew_file */

/*******************************************************************************
* mcdetector_statistics: compute detector statistics, error bars, [x I I_err N] 1D
* RETURN:            updated detector structure
* Used by: mcdetector_import
*******************************************************************************/
MCDETECTOR mcdetector_statistics(
  MCDETECTOR detector)
{

  if (!detector.p1 || !detector.m || !detector.filename)
    return(detector);
  
  /* compute statistics and update MCDETECTOR structure ===================== */
  double sum_z  = 0, min_z  = 0, max_z  = 0;
  double fmon_x =0,  smon_x = 0, fmon_y =0, smon_y=0, mean_z=0;
  double Nsum=0, P2sum=0;

  double sum_xz = 0, sum_yz = 0, sum_x = 0, sum_y = 0, sum_x2z = 0, sum_y2z = 0;
  int    i,j;
  char   hasnan=0, hasinf=0;
  char   israw = ((char*)strcasestr(detector.format,"raw") != NULL);
  double *this_p1=NULL; /* new 1D McCode array [x I E N]. Freed after writing data */

  /* if McCode/PGPLOT and rank==1 we create a new m*4 data block=[x I E N] */
  if (detector.rank == 1 && strcasestr(detector.format,"McCode")) {
    this_p1 = (double *)calloc(detector.m*detector.n*detector.p*4, sizeof(double));
    if (!this_p1)
      exit(-fprintf(stderr, "Error: Out of memory creating %li 1D " MCCODE_STRING " data set for file '%s' (mcdetector_import)\n",
        detector.m*detector.n*detector.p*4*sizeof(double*), detector.filename));
  }

  max_z = min_z = detector.p1[0];
  
  /* compute sum and moments (not for lists) */
  if (!strcasestr(detector.format,"list") && detector.m)
  for(j = 0; j < detector.n*detector.p; j++)
  {
    for(i = 0; i < detector.m; i++)
    {
      double x,y,z;
      double N, E;
      long   index= !detector.istransposed ? i*detector.n*detector.p + j : i+j*detector.m;
      char   hasnaninf=0;

      if (detector.m) 
        x = detector.xmin + (i + 0.5)/detector.m*(detector.xmax - detector.xmin); 
      else x = 0;
      if (detector.n && detector.p) 
        y = detector.ymin + (j + 0.5)/detector.n/detector.p*(detector.ymax - detector.ymin); 
      else y = 0;
      z = detector.p1[index];
      N = detector.p0 ? detector.p0[index] : 1;
      E = detector.p2 ? detector.p2[index] : 0;
      if (detector.p2 && !israw) 
        detector.p2[index] = (*mcestimate_error_p)(detector.p0[index],detector.p1[index],detector.p2[index]); /* set sigma */
      
      if (detector.rank == 1 && this_p1 && strcasestr(detector.format,"McCode")) {
        /* fill-in 1D McCode array [x I E N] */
        this_p1[index*4]   = x;
        this_p1[index*4+1] = z;
        this_p1[index*4+2] = detector.p2 ? detector.p2[index] : 0;
        this_p1[index*4+3] = N;
      }
      
      if (isnan(z) || isnan(E) || isnan(N)) hasnaninf=hasnan=1;
      if (isinf(z) || isinf(E) || isinf(N)) hasnaninf=hasinf=1;

      /* compute stats integrals */
      if (!hasnaninf) {
        sum_xz += x*z;
        sum_yz += y*z;
        sum_x  += x;
        sum_y  += y;
        sum_z  += z;
        sum_x2z += x*x*z;
        sum_y2z += y*y*z;
        if (z > max_z) max_z = z;
        if (z < min_z) min_z = z;

        Nsum += N;
        P2sum += E;
      }

    }
  } /* for j */

  /* compute 1st and 2nd moments. For lists, sum_z=0 so this is skipped. */
  if (sum_z && detector.n*detector.m*detector.p)
  {
    fmon_x = sum_xz/sum_z;
    fmon_y = sum_yz/sum_z;
    smon_x = sum_x2z/sum_z-fmon_x*fmon_x; smon_x = smon_x > 0 ? sqrt(smon_x) : 0;
    smon_y = sum_y2z/sum_z-fmon_y*fmon_y; smon_y = smon_y > 0 ? sqrt(smon_y) : 0;
    mean_z = sum_z/detector.n/detector.m/detector.p;
  }
  /* store statistics into detector */
  detector.intensity = sum_z;
  detector.error     = Nsum ? (*mcestimate_error_p)(Nsum, sum_z, P2sum) : 0;
  detector.events    = Nsum;
  detector.min       = min_z;
  detector.max       = max_z;
  detector.mean      = mean_z;
  detector.centerX   = fmon_x;
  detector.halfwidthX= smon_x;
  detector.centerY   = fmon_y;
  detector.halfwidthY= smon_y;

  /* if McCode/PGPLOT and rank==1 replace p1 with new m*4 1D McCode and clear others */
  if (detector.rank == 1 && this_p1 && strcasestr(detector.format,"McCode")) {
    
    detector.p1 = this_p1;
    detector.n  = detector.m; detector.m  = 4;
    detector.p0 = detector.p2 = NULL;
    detector.istransposed = 1;
  }

  if (detector.n*detector.m*detector.p > 1)
    snprintf(detector.signal, CHAR_BUF_LENGTH, 
      "Min=%g; Max=%g; Mean=%g;", detector.min, detector.max, detector.mean);
  else
    strcpy(detector.signal, "None");
  snprintf(detector.values, CHAR_BUF_LENGTH,
    "%g %g %g", detector.intensity, detector.error, detector.events);

  switch (detector.rank) {
    case 1:  snprintf(detector.statistics, CHAR_BUF_LENGTH, "X0=%g; dX=%g;",
      detector.centerX, detector.halfwidthX); break;
    case 2:
    case 3:  snprintf(detector.statistics, CHAR_BUF_LENGTH, "X0=%g; dX=%g; Y0=%g; dY=%g;",
      detector.centerX, detector.halfwidthX, detector.centerY, detector.halfwidthY);
      break;
    default: strcpy(detector.statistics, "None");
  }
  
  if (hasnan)
    printf("WARNING: Nan detected in component/file %s %s\n", 
      detector.component, strlen(detector.filename) ? detector.filename : "");
  if (hasinf)
    printf("WARNING: Inf detected in component/file %s %s\n", 
      detector.component, strlen(detector.filename) ? detector.filename : "");
  
  return(detector);
  
} /* mcdetector_statistics */

/*******************************************************************************
* mcdetector_import: build detector structure, merge non-lists from MPI
*                    compute basic stat, write "Detector:" line
* RETURN:            detector structure. Invalid data if detector.p1 == NULL
*                    Invalid detector sets m=0 and filename=""
*                    Simulation data  sets m=0 and filename=mcsiminfo_name
* This function is equivalent to the old 'mcdetector_out', returning a structure
*******************************************************************************/
MCDETECTOR mcdetector_import(
  char *format,
  char *component, char *title,
  long m, long n,  long p,
  char *xlabel, char *ylabel, char *zlabel,
  char *xvar, char *yvar, char *zvar,
  double x1, double x2, double y1, double y2, double z1, double z2,
  char *filename,
  double *p0, double *p1, double *p2,
  Coords position)
{
  time_t t;       /* for detector.date */
  long   date_l;  /* date as a long number */
  char   istransposed=0;
  char   c[CHAR_BUF_LENGTH]; /* temp var for signal label */

  MCDETECTOR detector;

  /* build MCDETECTOR structure ============================================= */
  /* make sure we do not have NULL for char fields */

  /* these also apply to simfile */
  strncpy (detector.filename,  filename ? filename : "",        CHAR_BUF_LENGTH);
  strncpy (detector.format,    format   ? format   : "McCode" , CHAR_BUF_LENGTH);
  /* add extension if missing */
  if (strlen(detector.filename) && !strchr(detector.filename, '.'))
  { /* add extension if not in file name already */
    strcat(detector.filename, ".dat");
  }
  strncpy (detector.component, component ? component : MCCODE_STRING " component", CHAR_BUF_LENGTH);

  snprintf(detector.instrument, CHAR_BUF_LENGTH, "%s (%s)", mcinstrument_name, mcinstrument_source);
  snprintf(detector.user, CHAR_BUF_LENGTH,      "%s on %s",
        getenv("USER") ? getenv("USER") : MCCODE_NAME,
        getenv("HOST") ? getenv("HOST") : "localhost");
  time(&t);         /* get current write time */
  date_l = (long)t; /* same but as a long */
  snprintf(detector.date, CHAR_BUF_LENGTH, "%s", ctime(&t));
  if (strlen(detector.date))   detector.date[strlen(detector.date)-1] = '\0'; /* remove last \n in date */
  detector.date_l = date_l;

  if (!mcget_run_num() || mcget_run_num() >= mcget_ncount())
    snprintf(detector.ncount, CHAR_BUF_LENGTH, "%llu", mcget_ncount()
#ifdef USE_MPI
*mpi_node_count
#endif
  );
  else
    snprintf(detector.ncount, CHAR_BUF_LENGTH, "%g/%g", (double)mcget_run_num(), (double)mcget_ncount());

  detector.p0         = p0;
  detector.p1         = p1;
  detector.p2         = p2;

  /* handle transposition (not for NeXus) */
  if (!strcasestr(detector.format, "NeXus")) {
    if (m<0 || n<0 || p<0)             istransposed = !istransposed;
    if (strcasestr(detector.format, "transpose")) istransposed = !istransposed;
    if (istransposed) { /* do the swap once for all */
      long i=m; m=n; n=i;
    }
  }

  m=labs(m); n=labs(n); p=labs(p); /* make sure dimensions are positive */
  detector.istransposed = istransposed;

  /* determine detector rank (dimensionality) */
  if (!m || !n || !p || !p1) detector.rank = 4; /* invalid: exit with m=0 filename="" */
  else if (m*n*p == 1)       detector.rank = 0; /* 0D */
  else if (n == 1 || m == 1) detector.rank = 1; /* 1D */
  else if (p == 1)           detector.rank = 2; /* 2D */
  else                       detector.rank = 3; /* 3D */

  /* from rank, set type */
  switch (detector.rank) {
    case 0:  strcpy(detector.type,  "array_0d"); m=n=p=1; break;
    case 1:  snprintf(detector.type, CHAR_BUF_LENGTH, "array_1d(%ld)", m*n*p); m *= n*p; n=p=1; break;
    case 2:  snprintf(detector.type, CHAR_BUF_LENGTH, "array_2d(%ld, %ld)", m, n*p); n *= p; p=1; break;
    case 3:  snprintf(detector.type, CHAR_BUF_LENGTH, "array_3d(%ld, %ld, %ld)", m, n, p); break;
    default: m=0; strcpy(detector.type, ""); strcpy(detector.filename, "");/* invalid */
  }

  detector.m    = m;
  detector.n    = n;
  detector.p    = p;

  /* these only apply to detector files ===================================== */

  snprintf(detector.position, CHAR_BUF_LENGTH, "%g %g %g", position.x, position.y, position.z);
  /* may also store actual detector orientation in the future */

  strncpy(detector.title,      title && strlen(title) ? title : component,       CHAR_BUF_LENGTH);
  strncpy(detector.xlabel,     xlabel && strlen(xlabel) ? xlabel : "X", CHAR_BUF_LENGTH); /* axis labels */
  strncpy(detector.ylabel,     ylabel && strlen(ylabel) ? ylabel : "Y", CHAR_BUF_LENGTH);
  strncpy(detector.zlabel,     zlabel && strlen(zlabel) ? zlabel : "Z", CHAR_BUF_LENGTH);
  strncpy(detector.xvar,       xvar && strlen(xvar) ? xvar :       "x", CHAR_BUF_LENGTH); /* axis variables */
  strncpy(detector.yvar,       yvar && strlen(yvar) ? yvar :       detector.xvar, CHAR_BUF_LENGTH);
  strncpy(detector.zvar,       zvar && strlen(zvar) ? zvar :       detector.yvar, CHAR_BUF_LENGTH);

  /* set "variables" as e.g. "I I_err N" */
  strcpy(c, "I ");
  if (strlen(detector.zvar))      strncpy(c, detector.zvar,32);
  else if (strlen(detector.yvar)) strncpy(c, detector.yvar,32);
  else if (strlen(detector.xvar)) strncpy(c, detector.xvar,32);

  if (detector.rank == 1)
    snprintf(detector.variables, CHAR_BUF_LENGTH, "%s %s %s_err N", detector.xvar, c, c);
  else
    snprintf(detector.variables, CHAR_BUF_LENGTH, "%s %s_err N", c, c);

  /* limits */
  detector.xmin = x1;
  detector.xmax = x2;
  detector.ymin = y1;
  detector.ymax = y2;
  detector.zmin = z1;
  detector.zmax = z2;
  if (abs(detector.rank) == 1)
    snprintf(detector.limits, CHAR_BUF_LENGTH, "%g %g", x1, x2);
  else if (detector.rank == 2)
    snprintf(detector.limits, CHAR_BUF_LENGTH, "%g %g %g %g", x1, x2, y1, y2);
  else
    snprintf(detector.limits, CHAR_BUF_LENGTH, "%g %g %g %g %g %g", x1, x2, y1, y2, z1, z2);

  /* if MPI and nodes_nb > 1: reduce data sets when using MPI =============== */
#ifdef USE_MPI
  if (!strcasestr(detector.format,"list") && mpi_node_count > 1 && m) {
    /* we save additive data: reduce everything into mpi_node_root */
    if (p0) mc_MPI_Sum(p0, m*n*p);
    if (p1) mc_MPI_Sum(p1, m*n*p);
    if (p2) mc_MPI_Sum(p2, m*n*p);
    if (!p0) {  /* additive signal must be then divided by the number of nodes */
      int i;
      for (i=0; i<m*n*p; i++) {
        p1[i] /= mpi_node_count;
        if (p2) p2[i] /= mpi_node_count;
      }
    }
  }
#endif /* USE_MPI */

  /* compute statistics, Nsum, intensity, Error bars */
  detector = mcdetector_statistics(detector);

#ifdef USE_MPI
  /* slaves are done */
  if(mpi_node_rank != mpi_node_root) {
    return detector;
  }
#endif

  /* output "Detector:" line ================================================ */
  /* when this is a detector written by a component (not the SAVE from instrument),
     not an event lists */
  if (!m) return(detector);
  if (!strcasestr(detector.format,"list")) {
    if (!strcmp(detector.component, mcinstrument_name)) {
      if (strlen(detector.filename))  /* we name it from its filename, or from its title */
        strncpy(c, detector.filename, CHAR_BUF_LENGTH);
      else
        snprintf(c, CHAR_BUF_LENGTH, "%s", mcinstrument_name);
    } else
      strncpy(c, detector.component, CHAR_BUF_LENGTH);  /* usual detectors written by components */

    printf("Detector: %s_I=%g %s_ERR=%g %s_N=%g",
           c, detector.intensity,
           c, detector.error,
           c, detector.events);
    printf(" \"%s\"\n", strlen(detector.filename) ? detector.filename : detector.component);
  }
  

  return(detector);
} /* mcdetector_import */

/* end MCDETECTOR import section ============================================ */

















/* ========================================================================== */

/*                               ASCII output                                 */
/*     The SIM file is YAML based, the data files have '#' headers            */

/* ========================================================================== */


/*******************************************************************************
* mcinfo_out: output instrument tags/info (only in SIM)
* Used in: mcsiminfo_init (ascii), mcinfo(stdout)
*******************************************************************************/
static void mcinfo_out(char *pre, FILE *f)
{
  char Parameters[CHAR_BUF_LENGTH] = "";
  int  i;

  if (!f || mcdisable_output_files) return;

  /* create parameter string ================================================ */
  for(i = 0; i < mcnumipar; i++)
  {
    char ThisParam[CHAR_BUF_LENGTH];
    if (strlen(mcinputtable[i].name) > CHAR_BUF_LENGTH) break;
    snprintf(ThisParam, CHAR_BUF_LENGTH, " %s(%s)", mcinputtable[i].name,
            (*mcinputtypes[mcinputtable[i].type].parminfo)
                (mcinputtable[i].name));
    strcat(Parameters, ThisParam);
    if (strlen(Parameters) >= CHAR_BUF_LENGTH-64) break;
  }

  /* output data ============================================================ */
  if (f != stdout)
    fprintf(f, "%sFile: %s%c%s\n",    pre, mcdirname, MC_PATHSEP_C, mcsiminfo_name);
  else
    fprintf(f, "%sCreator: %s\n",     pre, MCCODE_STRING);

  fprintf(f, "%sSource: %s\n",   pre, mcinstrument_source);
  fprintf(f, "%sParameters: %s\n",    pre, Parameters);
  
  fprintf(f, "%sTrace_enabled: %s\n", pre, mctraceenabled ? "yes" : "no");
  fprintf(f, "%sDefault_main: %s\n",  pre, mcdefaultmain ?  "yes" : "no");
  fprintf(f, "%sEmbedded_runtime: %s\n", pre, 
#ifdef MC_EMBEDDED_RUNTIME
         "yes"
#else
         "no"
#endif
         );

  fflush(f);
} /* mcinfo_out */

/*******************************************************************************
* mcruninfo_out_backend: output simulation tags/info (both in SIM and data files)
* Used in: mcsiminfo_init (ascii case), mcdetector_out_xD_ascii, mcinfo(stdout)
*******************************************************************************/
static void mcruninfo_out_backend(char *pre, FILE *f, int info)
{
  int i;
  char Parameters[CHAR_BUF_LENGTH];

  if (!f || mcdisable_output_files) return;

  fprintf(f, "%sFormat: %s%s\n",      pre, 
    mcformat && strlen(mcformat) ? mcformat : MCCODE_NAME,
    mcformat && strcasestr(mcformat,"McCode") ? " with text headers" : "");
  fprintf(f, "%sURL: %s\n",         pre, "http://www.mccode.org");
  fprintf(f, "%sCreator: %s\n",     pre, MCCODE_STRING);
  fprintf(f, "%sInstrument: %s\n", pre, mcinstrument_source);
  fprintf(f, "%sNcount: %llu\n",        pre, mcget_ncount());
  fprintf(f, "%sTrace: %s\n",       pre, mcdotrace ? "yes" : "no");
  fprintf(f, "%sGravitation: %s\n", pre, mcgravitation ? "yes" : "no");
  snprintf(Parameters, CHAR_BUF_LENGTH, "%ld", mcseed);
  fprintf(f, "%sSeed: %s\n",        pre, Parameters);
  fprintf(f, "%sDirectory: %s\n",        pre, mcdirname ? mcdirname : ".");
#ifdef USE_MPI
  if (mpi_node_count > 1)
    fprintf(f, "%sNodes: %i\n",        pre, mpi_node_count);
#endif

  /* output parameter string ================================================ */
  for(i = 0; i < mcnumipar; i++) {
      if (!info){
          (*mcinputtypes[mcinputtable[i].type].printer)(Parameters, mcinputtable[i].par);
          fprintf(f, "%sParam: %s=%s\n", pre, mcinputtable[i].name, Parameters);
      }else{
        /*if an info run, some variables might not have values. Flag these by "NULL"*/
	if(mcinputtable[i].val && strlen(mcinputtable[i].val)){
            /* ... those with defautl values*/
            (*mcinputtypes[mcinputtable[i].type].printer)(Parameters, mcinputtable[i].par);
            fprintf(f, "%sParam: %s=%s\n", pre, mcinputtable[i].name, Parameters);
        }else{
            /* ... and those without */
            fprintf(f, "%sParam: %s=NULL\n", pre, mcinputtable[i].name);
	}
      }
  }
} /* mcruninfo_out_backend */

/************************
* wrapper function to mcruninfo_out_backend
*  Regular runs use this whereas the single call from mcinfo is directly to the backend
*************************/
static void mcruninfo_out(char *pre, FILE *f){
    mcruninfo_out_backend(pre,f,0);
}

/*******************************************************************************
* mcsiminfo_out:    wrapper to fprintf(mcsiminfo_file)
*******************************************************************************/
void mcsiminfo_out(char *format, ...)
{
  va_list ap;

  if(mcsiminfo_file && !mcdisable_output_files)
  {
    va_start(ap, format);
    vfprintf(mcsiminfo_file, format, ap);
    va_end(ap);
  }
} /* mcsiminfo_out */


/*******************************************************************************
* mcdatainfo_out: output detector header
*   mcdatainfo_out(prefix, file_handle, detector) writes info to data file
*******************************************************************************/
static void
mcdatainfo_out(char *pre, FILE *f, MCDETECTOR detector)
{
  if (!f || !detector.m || mcdisable_output_files) return;
  
  /* output data ============================================================ */
  fprintf(f, "%sDate: %s (%li)\n",       pre, detector.date, detector.date_l);
  fprintf(f, "%stype: %s\n",       pre, detector.type);
  fprintf(f, "%sSource: %s\n",     pre, detector.instrument);
  fprintf(f, "%scomponent: %s\n",  pre, detector.component);
  fprintf(f, "%sposition: %s\n",   pre, detector.position);

  fprintf(f, "%stitle: %s\n",      pre, detector.title);
  fprintf(f, !mcget_run_num() || mcget_run_num() >= mcget_ncount() ?
             "%sNcount: %s\n" : 
             "%sratio: %s\n",  pre, detector.ncount);

  if (strlen(detector.filename)) {
    fprintf(f, "%sfilename: %s\n", pre, detector.filename);
  }

  fprintf(f, "%sstatistics: %s\n", pre, detector.statistics);
  fprintf(f, "%ssignal: %s\n",     pre, detector.signal);
  fprintf(f, "%svalues: %s\n",     pre, detector.values);

  if (detector.rank >= 1)
  {
    fprintf(f, "%sxvar: %s\n",     pre, detector.xvar);
    fprintf(f, "%syvar: %s\n",     pre, detector.yvar);
    fprintf(f, "%sxlabel: %s\n",   pre, detector.xlabel);
    fprintf(f, "%sylabel: %s\n",   pre, detector.ylabel);
    if (detector.rank > 1) {
      fprintf(f, "%szvar: %s\n",   pre, detector.zvar);
      fprintf(f, "%szlabel: %s\n", pre, detector.zlabel);
    }
  }

  fprintf(f, 
    abs(detector.rank)==1 ?
             "%sxlimits: %s\n" : 
             "%sxylimits: %s\n", pre, detector.limits);
  fprintf(f, "%svariables: %s\n", pre, 
    strcasestr(detector.format, "list") ? detector.ylabel : detector.variables);
    
  fflush(f);

} /* mcdatainfo_out */

/* mcdetector_out_array_ascii: output a single array to a file
 *   m: columns
 *   n: rows
 *   p: array
 *   f: file handle (already opened)
 */
static void mcdetector_out_array_ascii(long m, long n, double *p, FILE *f, char istransposed)
{
  if(f)
  {
    int i,j;
    for(j = 0; j < n; j++)
    {
      for(i = 0; i < m; i++)
      {
          fprintf(f, "%.10g ", p[!istransposed ? i*n + j : j*m+i]);
      }
      fprintf(f,"\n");
    }
  }
} /* mcdetector_out_array_ascii */

/*******************************************************************************
* mcdetector_out_0D_ascii: called by mcdetector_out_0D for ascii output
*******************************************************************************/
MCDETECTOR mcdetector_out_0D_ascii(MCDETECTOR detector)
{
  int exists=0;
  FILE *outfile = NULL;
  
  /* Write data set information to simulation description file. */
  MPI_MASTER(
    mcsiminfo_out("\nbegin data\n"); // detector.component
    mcdatainfo_out("  ", mcsiminfo_file, detector);
    mcsiminfo_out("end data\n");
    /* Don't write if filename is NULL: mcnew_file handles this (return NULL) */
    outfile = mcnew_file(detector.component, "dat", &exists);
    if(outfile)
    {
      /* write data file header and entry in simulation description file */
      mcruninfo_out( "# ", outfile);
      mcdatainfo_out("# ", outfile, detector);
      /* write I I_err N */
      fprintf(outfile, "%g %g %g\n", 
        detector.intensity, detector.error, detector.events);
      fclose(outfile);
    }
  ); /* MPI_MASTER */
  return(detector);
} /* mcdetector_out_0D_ascii */

/*******************************************************************************
* mcdetector_out_1D_ascii: called by mcdetector_out_1D for ascii output
*******************************************************************************/
MCDETECTOR mcdetector_out_1D_ascii(MCDETECTOR detector)
{
  int exists=0;
  FILE *outfile = NULL;

  MPI_MASTER(
    /* Write data set information to simulation description file. */
    mcsiminfo_out("\nbegin data\n"); // detector.filename
    mcdatainfo_out("  ", mcsiminfo_file, detector);
    mcsiminfo_out("end data\n");
    /* Loop over array elements, writing to file. */
    /* Don't write if filename is NULL: mcnew_file handles this (return NULL) */
    outfile = mcnew_file(detector.filename, "dat", &exists);
    if(outfile)
    {
      /* write data file header and entry in simulation description file */
      mcruninfo_out( "# ", outfile);
      mcdatainfo_out("# ", outfile, detector);
      /* output the 1D array columns */
      mcdetector_out_array_ascii(detector.m, detector.n, detector.p1, outfile, detector.istransposed);
      
      fclose(outfile);
    }
  ); /* MPI_MASTER */
  return(detector);
  
}  /* mcdetector_out_1D_ascii */

/*******************************************************************************
* mcdetector_out_2D_ascii: called by mcdetector_out_2D for ascii output
*******************************************************************************/
MCDETECTOR mcdetector_out_2D_ascii(MCDETECTOR detector)
{
  int exists=0;
  FILE *outfile = NULL;
  
  MPI_MASTER(
    /* Loop over array elements, writing to file. */
    /* Don't write if filename is NULL: mcnew_file handles this (return NULL) */
    outfile = mcnew_file(detector.filename, "dat", &exists);
    if(outfile)
    {
      /* write header only if file has just been created (not appending) */
      if (!exists) {
        /* Write data set information to simulation description file. */
        mcsiminfo_out("\nbegin data\n"); // detector.filename
        mcdatainfo_out("  ", mcsiminfo_file, detector);
        mcsiminfo_out("end data\n");
      
        mcruninfo_out( "# ", outfile);
        mcdatainfo_out("# ", outfile,   detector);
      }
      fprintf(outfile, "# Data [%s/%s] %s:\n", detector.component, detector.filename, detector.zvar);
      mcdetector_out_array_ascii(detector.m, detector.n*detector.p, detector.p1, 
        outfile, detector.istransposed);
      if (detector.p2) {
        fprintf(outfile, "# Errors [%s/%s] %s_err:\n", detector.component, detector.filename, detector.zvar);
        mcdetector_out_array_ascii(detector.m, detector.n*detector.p, detector.p2, 
          outfile, detector.istransposed);
      }
      if (detector.p0) {
        fprintf(outfile, "# Events [%s/%s] N:\n", detector.component, detector.filename);
        mcdetector_out_array_ascii(detector.m, detector.n*detector.p, detector.p0, 
          outfile, detector.istransposed);
      }
      fclose(outfile);
      
      if (!exists) {
        if (strcasestr(detector.format, "list"))
          printf("Events:   \"%s\"\n",  
            strlen(detector.filename) ? detector.filename : detector.component);
      }
    } /* if outfile */
  ); /* MPI_MASTER */
#ifdef USE_MPI
  if (strcasestr(detector.format, "list") && mpi_node_count > 1) {
    int node_i=0;
    /* loop along MPI nodes to write sequentially */
    for(node_i=0; node_i<mpi_node_count; node_i++) {
      /* MPI: slaves wait for the master to write its block, then append theirs */
      MPI_Barrier(MPI_COMM_WORLD);
      if (node_i != mpi_node_root && node_i == mpi_node_rank) {
        if(strlen(detector.filename) && !mcdisable_output_files)	/* Don't write if filename is NULL */
          outfile = mcnew_file(detector.filename, "dat", &exists);
        if (!exists)
          fprintf(stderr, "Warning: [MPI node %i] file '%s' does not exist yet, "
                          "MASTER should have opened it before.\n",
            mpi_node_rank, detector.filename);
        if(outfile) {
          mcdetector_out_array_ascii(detector.m, detector.n*detector.p, detector.p1, 
            outfile, detector.istransposed);
          fclose(outfile);
        }
      }
    }
  } /* if strcasestr list */
#endif
  return(detector);
} /* mcdetector_out_2D_ascii */

/*******************************************************************************
* strcpy_valid: makes a valid string for variable names.
*   copy 'original' into 'valid', replacing invalid characters by '_'
*   char arrays must be pre-allocated
*******************************************************************************/
static char *strcpy_valid(char *valid, char *original)
{
  long i;
  int  n=32; /* max length of valid names */

  if (original == NULL || !strlen(original)) return(NULL);

  if (n > strlen(original)) n = strlen(original);
  else original += strlen(original)-n;
  strncpy(valid, original, n);

  for (i=0; i < n; i++)
  {
    if ( (valid[i] > 122)
      || (valid[i] < 32)
      || (strchr("!\"#$%&'()*+,-.:;<=>?@[\\]^`/ \n\r\t", valid[i]) != NULL) )
    {
      if (i) valid[i] = '_'; else valid[i] = 'm';
    }
  }
  valid[i] = '\0';

  return(valid);
} /* strcpy_valid */

/* end ascii output section ================================================= */







#ifdef USE_NEXUS

/* ========================================================================== */

/*                               NeXus output                                 */

/* ========================================================================== */

#define nxprintf(...)    nxstr('d', __VA_ARGS__)
#define nxprintattr(...) nxstr('a', __VA_ARGS__)

/*******************************************************************************
* nxstr: output a tag=value data set (char) in NeXus/current group
*   when 'format' is larger that 1024 chars it is used as value for the 'tag'
*   else the value is assembled with format and following arguments.
*   type='d' -> data set
*        'a' -> attribute for current data set
*******************************************************************************/
static int nxstr(char type, NXhandle *f, char *tag, char *format, ...)
{
  va_list ap;
  char value[CHAR_BUF_LENGTH];
  int  i;
  int  ret=NX_OK;
  
  if (!tag || !format || !strlen(tag) || !strlen(format)) return(NX_OK);
  
  /* assemble the value string */
  if (strlen(format) < CHAR_BUF_LENGTH) {
    va_start(ap, format);
    ret = vsnprintf(value, CHAR_BUF_LENGTH, format, ap);
    va_end(ap);
  
    i = strlen(value);
  } else {
    i = strlen(format);
  }

  if (type == 'd') {
    /* open/put/close data set */
    if (NXmakedata (f, tag, NX_CHAR, 1, &i) != NX_OK) return(NX_ERROR);
    NXopendata (f, tag);
    if (strlen(format) < CHAR_BUF_LENGTH)
      ret = NXputdata  (f, value);
    else
      ret = NXputdata  (f, format);
    NXclosedata(f);
  } else {
    if (strlen(format) < CHAR_BUF_LENGTH)
      ret = NXputattr  (f, tag, value, strlen(value), NX_CHAR);
    else
      ret = NXputattr  (f, tag, format, strlen(format), NX_CHAR);
  }
  
  return(ret);
  
} /* nxstr */

/*******************************************************************************
* mcinfo_readfile: read a full file into a string buffer which is allocated
*   Think to free the buffer after use.
* Used in: mcinfo_out_nexus (nexus)
*******************************************************************************/
char *mcinfo_readfile(char *filename)
{
  FILE *f = fopen(filename, "rb");
  if (!f) return(NULL);
  fseek(f, 0, SEEK_END);
  long fsize = ftell(f);
  rewind(f);
  char *string = malloc(fsize + 1);
  if (string) {
    int n = fread(string, fsize, 1, f);
    fclose(f);

    string[fsize] = 0;
  }
  return(string);
}

/*******************************************************************************
* mcinfo_out: output instrument/simulation groups in NeXus file
* Used in: mcsiminfo_init (nexus)
*******************************************************************************/
static void mcinfo_out_nexus(NXhandle f)
{
  FILE  *fid;     /* for intrument source code/C/IDF */
  char  *buffer=NULL;
  time_t t     =time(NULL); /* for date */
  char   entry0[CHAR_BUF_LENGTH];
  int    count=0;
  char   name[CHAR_BUF_LENGTH];
  char   class[CHAR_BUF_LENGTH];
  
  if (!f || mcdisable_output_files) return;
  
  /* write NeXus NXroot attributes */
  /* automatically added: file_name, HDF5_Version, file_time, NeXus_version */ 
  nxprintattr(f, "creator",   "%s generated with " MCCODE_STRING, mcinstrument_name);
  
  /* count the number of existing NXentry and create the next one */
  NXgetgroupinfo(f, &count, name, class);
  sprintf(entry0, "entry%i", count+1);

  /* create the main NXentry (mandatory in NeXus) */
  if (NXmakegroup(f, entry0, "NXentry") == NX_OK) 
  if (NXopengroup(f, entry0, "NXentry") == NX_OK) {
    
    nxprintf(nxhandle, "program_name", MCCODE_STRING);
    nxprintf(f, "start_time", ctime(&t));
    nxprintf(f, "title", "%s%s%s simulation generated by instrument %s", 
      mcdirname && strlen(mcdirname) ? mcdirname : ".", MC_PATHSEP_S, mcsiminfo_name,
      mcinstrument_name);
    nxprintattr(f, "program_name", MCCODE_STRING);
    nxprintattr(f, "instrument",   mcinstrument_name);
    nxprintattr(f, "simulation",   "%s%s%s",
        mcdirname && strlen(mcdirname) ? mcdirname : ".", MC_PATHSEP_S, mcsiminfo_name);

    /* write NeXus instrument group */
    if (NXmakegroup(f, "instrument", "NXinstrument") == NX_OK)
    if (NXopengroup(f, "instrument", "NXinstrument") == NX_OK) {
      int   i;
      char *string=NULL;

      /* write NeXus parameters(types) data =================================== */
      string = (char*)malloc(CHAR_BUF_LENGTH);
      if (string) {
        strcpy(string, "");
        for(i = 0; i < mcnumipar; i++)
        {
          char ThisParam[CHAR_BUF_LENGTH];
          snprintf(ThisParam, CHAR_BUF_LENGTH, " %s(%s)", mcinputtable[i].name,
                  (*mcinputtypes[mcinputtable[i].type].parminfo)
                      (mcinputtable[i].name));
          if (strlen(string) + strlen(ThisParam) < CHAR_BUF_LENGTH)
            strcat(string, ThisParam);
        }
        nxprintattr(f, "Parameters",    string);
        free(string);
      }
        
      nxprintattr(f, "name",          mcinstrument_name);
      nxprintf   (f, "name",          mcinstrument_name);
      nxprintattr(f, "Source",        mcinstrument_source);
      
      nxprintattr(f, "Trace_enabled", mctraceenabled ? "yes" : "no");
      nxprintattr(f, "Default_main",  mcdefaultmain ?  "yes" : "no");
      nxprintattr(f, "Embedded_runtime",  
  #ifdef MC_EMBEDDED_RUNTIME
           "yes"
  #else
           "no"
  #endif
           );
           
      /* add instrument source code when available */
      buffer = mcinfo_readfile(mcinstrument_source);
      if (buffer && strlen(buffer)) {
        long length=strlen(buffer);
        nxprintf (f, "description", buffer);
        NXopendata(f,"description");
        nxprintattr(f, "file_name", mcinstrument_source);
        nxprintattr(f, "file_size", "%li", length);
        nxprintattr(f, "MCCODE_STRING", MCCODE_STRING);
        NXclosedata(f);
        nxprintf (f,"instrument_source", "%s " MCCODE_NAME " " MCCODE_PARTICLE " Monte Carlo simulation", mcinstrument_name);
        free(buffer);
      } else
        nxprintf (f, "description", "File %s not found (instrument description %s is missing)", 
          mcinstrument_source, mcinstrument_name);
      
      /* add Mantid/IDF.xml when available */
      char *IDFfile=NULL;
      IDFfile = (char*)malloc(CHAR_BUF_LENGTH);
      sprintf(IDFfile,"%s%s",mcinstrument_source,".xml");
      buffer = mcinfo_readfile(IDFfile);
      if (buffer && strlen(buffer)) {
        NXmakegroup (nxhandle, "instrument_xml", "NXnote");
        NXopengroup (nxhandle, "instrument_xml", "NXnote");
        nxprintf(f, "data", buffer);
        nxprintf(f, "description", "IDF.xml file found with instrument %s", mcinstrument_source);
        nxprintf(f, "type", "text/xml");
        NXclosegroup(f); /* instrument_xml */
        free(buffer);
      }
      free(IDFfile);
      NXclosegroup(f); /* instrument */
    } /* NXinstrument */

    /* write NeXus simulation group */
    if (NXmakegroup(f, "simulation", "NXnote") == NX_OK)
    if (NXopengroup(f, "simulation", "NXnote") == NX_OK) {

      nxprintattr(f, "name",   "%s%s%s",
        mcdirname && strlen(mcdirname) ? mcdirname : ".", MC_PATHSEP_S, mcsiminfo_name);
      
      nxprintf   (f, "name",      "%s",     mcsiminfo_name);
      nxprintattr(f, "Format",    mcformat && strlen(mcformat) ? mcformat : MCCODE_NAME);
      nxprintattr(f, "URL",       "http://www.mccode.org");
      nxprintattr(f, "program",   MCCODE_STRING);
      nxprintattr(f, "Instrument",mcinstrument_source);
      nxprintattr(f, "Trace",     mcdotrace ?     "yes" : "no");
      nxprintattr(f, "Gravitation",mcgravitation ? "yes" : "no");
      nxprintattr(f, "Seed",      "%li", mcseed);
      nxprintattr(f, "Directory", mcdirname);
    #ifdef USE_MPI
      if (mpi_node_count > 1)
        nxprintf(f, "Nodes", "%i",        mpi_node_count);
    #endif
    
      /* output parameter string ================================================ */
      if (NXmakegroup(f, "Param", "NXparameters") == NX_OK)
      if (NXopengroup(f, "Param", "NXparameters") == NX_OK) {
        int i;
        char string[CHAR_BUF_LENGTH];
        for(i = 0; i < mcnumipar; i++) {
          if (mcget_run_num() || (mcinputtable[i].val && strlen(mcinputtable[i].val))) {
            if (mcinputtable[i].par == NULL)
              strncpy(string, (mcinputtable[i].val ? mcinputtable[i].val : ""), CHAR_BUF_LENGTH);
            else
              (*mcinputtypes[mcinputtable[i].type].printer)(string, mcinputtable[i].par);

            nxprintf(f,  mcinputtable[i].name, "%s", string);
            nxprintattr(f, mcinputtable[i].name, string);
          }
        }
        NXclosegroup(f); /* Param */
      } /* NXparameters */
      
      NXclosegroup(f); /* simulation */
    } /* NXsimulation */
    
    /* create a group to hold all monitors */
    NXmakegroup(f, "data", "NXdetector");

    /* leave the NXentry opened (closed at exit) */
  } /* NXentry */
} /* mcinfo_out_nexus */

/*******************************************************************************
* mcdatainfo_out_nexus: output detector header
*   mcdatainfo_out_nexus(detector) create group and write info to NeXus data file
*   open data:NXdetector then filename:NXdata and write headers/attributes
*   requires: NXentry to be opened
*******************************************************************************/
static void
mcdatainfo_out_nexus(NXhandle f, MCDETECTOR detector)
{
  char data_name[32];
  if (!f || !detector.m || mcdisable_output_files) return;
  
  strcpy_valid(data_name, 
    detector.filename && strlen(detector.filename) ? 
      detector.filename : detector.component);

  /* the NXdetector group has been created in mcinfo_out_nexus (mcsiminfo_init) */
  if (NXopengroup(f, "data", "NXdetector") == NX_OK) {

    /* create and open the data group */
    /* this may fail when appending to list -> ignore/skip */
    NXMDisableErrorReporting(); /* unactivate NeXus error messages, as creation may fail */
    
    if (NXmakegroup(f, data_name, "NXdata") == NX_OK)
    if (NXopengroup(f, data_name, "NXdata") == NX_OK) {
    
      /* output metadata (as attributes) ======================================== */
      nxprintattr(f, "Date",       detector.date);
      nxprintattr(f, "type",       detector.type);
      nxprintattr(f, "Source",     detector.instrument);
      nxprintattr(f, "component",  detector.component);
      nxprintattr(f, "position",   detector.position);

      nxprintattr(f, "title",      detector.title);
      nxprintattr(f, !mcget_run_num() || mcget_run_num() >= mcget_ncount() ?
                 "Ncount" : 
                 "ratio",  detector.ncount);

      if (strlen(detector.filename)) {
        nxprintattr(f, "filename", detector.filename);
      }

      nxprintattr(f, "statistics", detector.statistics);
      nxprintattr(f, "signal",     detector.signal);
      nxprintattr(f, "values",     detector.values);

      if (detector.rank >= 1)
      {
        nxprintattr(f, "xvar",     detector.xvar);
        nxprintattr(f, "yvar",     detector.yvar);
        nxprintattr(f, "xlabel",   detector.xlabel);
        nxprintattr(f, "ylabel",   detector.ylabel);
        if (detector.rank > 1) {
          nxprintattr(f, "zvar",   detector.zvar);
          nxprintattr(f, "zlabel", detector.zlabel);
        }
      }

      nxprintattr(f, abs(detector.rank)==1 ?
                 "xlimits" : 
                 "xylimits", detector.limits);
      nxprintattr(f, "variables", 
        strcasestr(detector.format, "list") ? detector.ylabel : detector.variables);
      nxprintf(f, "distance", detector.position);
      nxprintf(f, "acquisition_mode",
        strcasestr(detector.format, "list") ? "event" : "summed");
        
      NXclosegroup(f);
    } /* NXdata (filename) */
    NXMEnableErrorReporting();  /* re-enable NeXus error messages */
    NXclosegroup(f);
  } /* NXdetector (data) */
  
} /* mcdatainfo_out_nexus */

/*******************************************************************************
* mcdetector_out_axis_nexus: write detector axis into current NXdata
*   requires: NXdata to be opened
*******************************************************************************/
int mcdetector_out_axis_nexus(NXhandle f, char *label, char *var, int rank, long length, double min, double max)
{
  if (!f || length <= 1 || mcdisable_output_files || max == min) return(NX_OK);
  else {
    double axis[length];
    char valid[32];
    int dim=(int)length;
    int i;
    int nprimary=1;
    /* create an axis from [min:max] */
    for(i = 0; i < length; i++)
      axis[i] = min+(max-min)*(i+0.5)/length;
    /* create the data set */
    strcpy_valid(valid, label);
    NXcompmakedata(f, valid, NX_FLOAT64, 1, &dim, NX_COMP_LZW, &dim);
    /* open it */
    if (NXopendata(f, valid) != NX_OK) {
      fprintf(stderr, "Warning: could not open axis rank %i '%s' (NeXus)\n",
        rank, valid);
      return(NX_ERROR);
    }
    /* put the axis and its attributes */
    NXputdata  (f, axis);
    nxprintattr(f, "long_name",  label);
    nxprintattr(f, "short_name", var);
    NXputattr  (f, "axis",       &rank,     1, NX_INT32);
    nxprintattr(f, "units",      var);
    NXputattr  (f, "primary",    &nprimary, 1, NX_INT32);
    NXclosedata(f);
    
    return(NX_OK);
  }
} /* mcdetector_out_axis_nexus */

/*******************************************************************************
* mcdetector_out_array_nexus: write detector array into current NXdata (1D,2D)
*   requires: NXdata to be opened
*******************************************************************************/
int mcdetector_out_array_nexus(NXhandle f, char *part, double *data, MCDETECTOR detector)
{
  
  int dims[3]={detector.m,detector.n,detector.p};  /* number of elements to write */
  int signal=1;
  int exists=0;
  int current_dims[3]={0,0,0};
  int ret=NX_OK;
  
  if (!f || !data || !detector.m || mcdisable_output_files) return(NX_OK);
  
  /* when this is a list, we set 1st dimension to NX_UNLIMITED for creation */
  if (strcasestr(detector.format, "list")) dims[0] = NX_UNLIMITED;
  
  /* create the data set in NXdata group */
  NXMDisableErrorReporting(); /* unactivate NeXus error messages, as creation may fail */
  /* NXcompmakedata fails with NX_UNLIMITED */
  if (strcasestr(detector.format, "list"))
    ret = NXmakedata(    f, part, NX_FLOAT64, detector.rank, dims);
  else
    ret = NXcompmakedata(f, part, NX_FLOAT64, detector.rank, dims, NX_COMP_LZW, dims);
  if (ret != NX_OK) {
    /* failed: data set already exists */
    int datatype=0;
    int rank=0;
    exists=1;
    /* inquire current size of data set (nb of events stored) */
    NXopendata(f, part);
    NXgetinfo(f, &rank, current_dims, &datatype);
    NXclosedata(f);
  }
  NXMEnableErrorReporting();  /* re-enable NeXus error messages */
  dims[0] = detector.m; /* restore actual dimension from data writing */
  
  /* open the data set */
  if (NXopendata(f, part) == NX_ERROR) {
    fprintf(stderr, "Warning: could not open DataSet %s '%s' (NeXus)\n",
      part, detector.title);
    return(NX_ERROR);
  }
  if (strcasestr(detector.format, "list")) {
    current_dims[1] = current_dims[2] = 0; /* set starting location for writing slab */
    NXputslab(f, data, current_dims, dims);
    if (!exists)
      printf("Events:   \"%s\"\n",  
        strlen(detector.filename) ? detector.filename : detector.component);
  } else {
    NXputdata (f, data);
  }
  
  if (strstr(part,"data") || strstr(part, "events")) {
    NXputattr(f, "signal", &signal, 1, NX_INT32);
    nxprintattr(f, "short_name", detector.filename && strlen(detector.filename) ? 
      detector.filename : detector.component);
  }
  nxprintattr(f, "long_name", "%s '%s'", part, detector.title);
  NXclosedata(f);
  
  return(NX_OK);
} /* mcdetector_out_array_nexus */

/*******************************************************************************
* mcdetector_out_data_nexus: write detector axes+data into current NXdata
*   The data:NXdetector is opened, then filename:NXdata
*   requires: NXentry to be opened
*******************************************************************************/
int mcdetector_out_data_nexus(NXhandle f, MCDETECTOR detector)
{
  char data_name[32];
  
  if (!f || !detector.m || mcdisable_output_files) return(NX_OK);
  
  strcpy_valid(data_name, 
    detector.filename && strlen(detector.filename) ? 
      detector.filename : detector.component);

  /* the NXdetector group has been created in mcinfo_out_nexus (mcsiminfo_init) */
  if (NXopengroup(f, "data", "NXdetector") == NX_OK) {

    /* the NXdata group has been created in mcdatainfo_out_nexus */
    if (NXopengroup(f, data_name, "NXdata") == NX_OK) {
  
      /* write axes, for histogram data sets, not for lists */
      if (!strcasestr(detector.format, "list")) {
        mcdetector_out_axis_nexus(f, detector.xlabel, detector.xvar, 
          1, detector.m, detector.xmin, detector.xmax);
          
        mcdetector_out_axis_nexus(f, detector.ylabel, detector.yvar, 
          2, detector.n, detector.ymin, detector.ymax);
          
        mcdetector_out_axis_nexus(f, detector.zlabel, detector.zvar, 
          3, detector.p, detector.zmin, detector.zmax);

      } /* !list */
      
      /* write the actual data (appended if already exists) */
      if (!strcasestr(detector.format, "list")) {
        mcdetector_out_array_nexus(f, "data", detector.p1, detector);
        mcdetector_out_array_nexus(f, "errors", detector.p2, detector);
        mcdetector_out_array_nexus(f, "ncount", detector.p0, detector);
      } else
        mcdetector_out_array_nexus(  f, "events", detector.p1, detector);
      
      NXclosegroup(f);
    } /* NXdata */
    NXclosegroup(f);
  } /* NXdetector */
  
  return(NX_OK);
} /* mcdetector_out_array_nexus */

#ifdef USE_MPI
/*******************************************************************************
* mcdetector_out_list_slaves: slaves send their list data to master which writes
*   requires: NXentry to be opened
* WARNING: this method has a flaw: it requires all nodes to flush the lists
*   the same number of times. In case one node is just below the buffer size
*   when finishing (e.g. monitor_nd), it may not trigger save but others may. 
*   Then the number of recv/send is not constant along nodes, and simulation stalls.  
*******************************************************************************/
MCDETECTOR mcdetector_out_list_slaves(MCDETECTOR detector)
{
  int     node_i=0;
  MPI_MASTER(
	     printf("\n** MPI master gathering slave node list data ** \n");
  );
  
  if (mpi_node_rank != mpi_node_root) {
    /* MPI slave: slaves send their data to master: 2 MPI_Send calls */
    /* m, n, p must be sent first, since all slaves do not have the same number of events */
    int mnp[3]={detector.m,detector.n,detector.p};

    if (mc_MPI_Send(mnp, 3, MPI_INT, mpi_node_root)!= MPI_SUCCESS)
      fprintf(stderr, "Warning: proc %i to master: MPI_Send mnp list error (mcdetector_out_list_slaves)\n", mpi_node_rank);
    if (!detector.p1
     || mc_MPI_Send(detector.p1, mnp[0]*mnp[1]*mnp[2], MPI_DOUBLE, mpi_node_root) != MPI_SUCCESS)
      fprintf(stderr, "Warning: proc %i to master: MPI_Send p1 list error: mnp=%i (mcdetector_out_list_slaves)\n", mpi_node_rank, abs(mnp[0]*mnp[1]*mnp[2]));
    /* slaves are done: sent mnp and p1 */
    return (detector);
  } /* end slaves */

  /* MPI master: receive data from slaves sequentially: 2 MPI_Recv calls */

  if (mpi_node_rank == mpi_node_root) {
    for(node_i=0; node_i<mpi_node_count; node_i++) {
      double *this_p1=NULL;                               /* buffer to hold the list from slaves */
      int     mnp[3]={0,0,0};  /* size of this buffer */
      if (node_i != mpi_node_root) { /* get data from slaves */
	if (mc_MPI_Recv(mnp, 3, MPI_INT, node_i) != MPI_SUCCESS)
	  fprintf(stderr, "Warning: master from proc %i: "
		  "MPI_Recv mnp list error (mcdetector_write_data)\n", node_i);
	if (mnp[0]*mnp[1]*mnp[2]) {
	  this_p1 = (double *)calloc(mnp[0]*mnp[1]*mnp[2], sizeof(double));
	  if (!this_p1 || mc_MPI_Recv(this_p1, abs(mnp[0]*mnp[1]*mnp[2]), MPI_DOUBLE, node_i)!= MPI_SUCCESS)
	    fprintf(stderr, "Warning: master from proc %i: "
		    "MPI_Recv p1 list error: mnp=%i (mcdetector_write_data)\n", node_i, mnp[0]*mnp[1]*mnp[2]);
	  else {
	    printf(". MPI master writing data for slave node %i\n",node_i);
	    detector.p1 = this_p1;
	    detector.m  = mnp[0]; detector.n  = mnp[1]; detector.p  = mnp[2];
	    
	    mcdetector_out_data_nexus(nxhandle, detector);
	  }
	}
      } /* if not master */
    } /* for */
  MPI_MASTER(
	     printf("\n** Done ** \n");
  );   
  }
}
#endif

MCDETECTOR mcdetector_out_0D_nexus(MCDETECTOR detector)
{
  /* Write data set information to NeXus file. */
  MPI_MASTER(
    mcdatainfo_out_nexus(nxhandle, detector);
  );
  
  return(detector);
} /* mcdetector_out_0D_ascii */

MCDETECTOR mcdetector_out_1D_nexus(MCDETECTOR detector_inc)
{
  MCDETECTOR detector = detector_inc;
  MPI_MASTER(
  mcdatainfo_out_nexus(nxhandle, detector);
  mcdetector_out_data_nexus(nxhandle, detector);
  );
  return(detector);
} /* mcdetector_out_1D_ascii */

MCDETECTOR mcdetector_out_2D_nexus(MCDETECTOR detector_inc)
{
  MCDETECTOR detector = detector_inc;
  MPI_MASTER(
  mcdatainfo_out_nexus(nxhandle, detector);
  mcdetector_out_data_nexus(nxhandle, detector);
  );
  
#ifdef USE_MPI // and USE_NEXUS
  /* NeXus: slave nodes have master write their lists */
  if (strcasestr(detector.format, "list") && mpi_node_count > 1) {
    mcdetector_out_list_slaves(detector);
  }
#endif /* USE_MPI */

  return(detector);
} /* mcdetector_out_2D_nexus */

#endif /* USE_NEXUS*/








/* ========================================================================== */

/*                            Main input functions                            */
/*            DETECTOR_OUT_xD function calls -> ascii or NeXus                */

/* ========================================================================== */

/*******************************************************************************
* mcsiminfo_init:   open SIM and write header
*******************************************************************************/
FILE *mcsiminfo_init(FILE *f)
{
  int exists=0;
  int index;
  
  /* check format */      
  if (!mcformat || !strlen(mcformat) 
   || !strcasecmp(mcformat, "MCSTAS") || !strcasecmp(mcformat, "MCXTRACE") 
   || !strcasecmp(mcformat, "PGPLOT") || !strcasecmp(mcformat, "GNUPLOT") || !strcasecmp(mcformat, "MCCODE")
   || !strcasecmp(mcformat, "MATLAB")) {
    mcformat="McCode";
#ifdef USE_NEXUS
  } else if (strcasestr(mcformat, "NeXus")) {
    /* Do nothing */
#endif
  } else {
    fprintf(stderr,
	    "Warning: You have requested the output format %s which is unsupported by this binary. Resetting to standard %s format.\n",mcformat ,"McCode");
    mcformat="McCode";
  }
  
  /* open the SIM file if not defined yet */
  if (mcsiminfo_file || mcdisable_output_files) 
    return (mcsiminfo_file);
    
#ifdef USE_NEXUS
  /* only master writes NeXus header: calls NXopen(nxhandle) */
  if (mcformat && strcasestr(mcformat, "NeXus")) {
	  MPI_MASTER(
	  mcsiminfo_file = mcnew_file(mcsiminfo_name, "h5", &exists);
    if(!mcsiminfo_file)
      fprintf(stderr,
	      "Warning: could not open simulation description file '%s'\n",
	      mcsiminfo_name);
	  else
	    mcinfo_out_nexus(nxhandle);
	  );
    return(mcsiminfo_file); /* points to nxhandle */
  }
#endif
  
  /* write main description file (only MASTER) */
  MPI_MASTER(

  mcsiminfo_file = mcnew_file(mcsiminfo_name, "sim", &exists);
  if(!mcsiminfo_file)
    fprintf(stderr,
	    "Warning: could not open simulation description file '%s'\n",
	    mcsiminfo_name);
  else
  {
    /* write SIM header */
    time_t t=time(NULL);
    mcsiminfo_out("%s simulation description file for %s.\n", 
      MCCODE_NAME, mcinstrument_name);
    mcsiminfo_out("Date:    %s", ctime(&t)); /* includes \n */
    mcsiminfo_out("Program: %s\n\n", MCCODE_STRING);
    
    mcsiminfo_out("begin instrument: %s\n", mcinstrument_name);
    mcinfo_out(   "  ", mcsiminfo_file);
    mcsiminfo_out("end instrument\n");

    mcsiminfo_out("\nbegin simulation: %s\n", mcdirname);
    mcruninfo_out("  ", mcsiminfo_file);
    mcsiminfo_out("end simulation\n");

  }
  return (mcsiminfo_file);
  
  ); /* MPI_MASTER */
  
} /* mcsiminfo_init */

/*******************************************************************************
*   mcsiminfo_close:  close SIM
*******************************************************************************/
void mcsiminfo_close()
{
  MPI_MASTER(
  if(mcsiminfo_file && !mcdisable_output_files) {
#ifdef USE_NEXUS
    if (mcformat && strcasestr(mcformat, "NeXus")) {
      time_t t=time(NULL);
      nxprintf(nxhandle, "end_time", ctime(&t));
      nxprintf(nxhandle, "duration", "%li", (long)t-mcstartdate);
      NXclosegroup(nxhandle); /* NXentry */
      NXclose(&nxhandle);
    } else
#endif
      fclose(mcsiminfo_file);
    );
    mcsiminfo_file = NULL;
  }
} /* mcsiminfo_close */

/*******************************************************************************
* mcdetector_out_0D: wrapper for 0D (single value).
*   Output single detector/monitor data (p0, p1, p2).
*   Title is t, component name is c.
*******************************************************************************/
MCDETECTOR mcdetector_out_0D(char *t, double p0, double p1, double p2,
                         char *c, Coords posa)
{
  /* import and perform basic detector analysis (and handle MPI reduce) */
  MCDETECTOR detector = mcdetector_import(mcformat,
    c, (t ? t : MCCODE_STRING " data"),
    1, 1, 1,
    "I", "", "",
    "I", "", "",
    0, 0, 0, 0, 0, 0, "",
    &p0, &p1, &p2, posa); /* write Detector: line */

#ifdef USE_NEXUS
  if (strcasestr(detector.format, "NeXus"))
    return(mcdetector_out_0D_nexus(detector));
  else
#endif
    return(mcdetector_out_0D_ascii(detector));
    
} /* mcdetector_out_0D */



/*******************************************************************************
* mcdetector_out_1D: wrapper for 1D.
*   Output 1d detector data (p0, p1, p2) for n bins linearly
*   distributed across the range x1..x2 (x1 is lower limit of first
*   bin, x2 is upper limit of last bin). Title is t, axis labels are xl
*   and yl. File name is f, component name is c.
*******************************************************************************/
MCDETECTOR mcdetector_out_1D(char *t, char *xl, char *yl,
        char *xvar, double x1, double x2,
        long n,
        double *p0, double *p1, double *p2, char *f,
        char *c, Coords posa)
{
  /* import and perform basic detector analysis (and handle MPI_Reduce) */
  MCDETECTOR detector = mcdetector_import(mcformat,
    c, (t ? t : MCCODE_STRING " 1D data"),
    n, 1, 1,
    xl, yl, (n > 1 ? "Signal per bin" : " Signal"),
    xvar, "(I,I_err)", "I",
    x1, x2, 0, 0, 0, 0, f,
    p0, p1, p2, posa); /* write Detector: line */
  if (!detector.p1 || !detector.m) return(detector);

#ifdef USE_NEXUS
  if (strcasestr(detector.format, "NeXus"))
    return(mcdetector_out_1D_nexus(detector));
  else
#endif
    return(mcdetector_out_1D_ascii(detector));
  
} /* mcdetector_out_1D */

/*******************************************************************************
* mcdetector_out_2D: wrapper for 2D.
*   special case for list: master creates file first, then slaves append their blocks without header
*******************************************************************************/
MCDETECTOR mcdetector_out_2D(char *t, char *xl, char *yl,
                  double x1, double x2, double y1, double y2,
                  long m, long n,
                  double *p0, double *p1, double *p2, char *f,
                  char *c, Coords posa)
{
  char xvar[CHAR_BUF_LENGTH];
  char yvar[CHAR_BUF_LENGTH];
  
  /* create short axes labels */
  if (xl && strlen(xl)) { strncpy(xvar, xl, CHAR_BUF_LENGTH); xvar[strcspn(xvar,"\n\r ")]='\0'; }
  else strcpy(xvar, "x");
  if (yl && strlen(yl)) { strncpy(yvar, yl, CHAR_BUF_LENGTH); yvar[strcspn(yvar,"\n\r ")]='\0'; }
  else strcpy(yvar, "y");

  MCDETECTOR detector;

  /* import and perform basic detector analysis (and handle MPI_Reduce) */
  if (labs(m) == 1) {/* n>1 on Y, m==1 on X: 1D, no X axis*/
    detector = mcdetector_import(mcformat,
      c, (t ? t : MCCODE_STRING " 1D data"),
      n, 1, 1,
      yl, "", "Signal per bin",
      yvar, "(I,Ierr)", "I",
      y1, y2, x1, x2, 0, 0, f,
      p0, p1, p2, posa); /* write Detector: line */
  } else if (labs(n)==1) {/* m>1 on X, n==1 on Y: 1D, no Y axis*/
    detector = mcdetector_import(mcformat,
      c, (t ? t : MCCODE_STRING " 1D data"),
      m, 1, 1,
      xl, "", "Signal per bin",
      xvar, "(I,Ierr)", "I",
      x1, x2, y1, y2, 0, 0, f,
      p0, p1, p2, posa); /* write Detector: line */
  }else {
    detector = mcdetector_import(mcformat,
      c, (t ? t : MCCODE_STRING " 2D data"),
      m, n, 1,
      xl, yl, "Signal per bin",
      xvar, yvar, "I",
      x1, x2, y1, y2, 0, 0, f,
      p0, p1, p2, posa); /* write Detector: line */
  }

  if (!detector.p1 || !detector.m) return(detector);

#ifdef USE_NEXUS
  if (strcasestr(detector.format, "NeXus"))
    return(mcdetector_out_2D_nexus(detector));
  else
#endif
    return(mcdetector_out_2D_ascii(detector));
  
} /* mcdetector_out_2D */

/*******************************************************************************
* mcdetector_out_list: wrapper for list output (calls out_2D with mcformat+"list").
*   m=number of events, n=size of each event
*******************************************************************************/
MCDETECTOR mcdetector_out_list(char *t, char *xl, char *yl,
                  long m, long n,
                  double *p1, char *f,
                  char *c, Coords posa)
{
  char       format_new[CHAR_BUF_LENGTH];
  char      *format_org;
  MCDETECTOR detector;
  
  format_org = mcformat;
  strcpy(format_new, mcformat);
  strcat(format_new, " list");
  mcformat = format_new;

  detector = mcdetector_out_2D(t, xl, yl,
                  1,labs(m),1,labs(n),
                  m,n,
                  NULL, p1, NULL, f,
                  c, posa);
  
  mcformat = format_org;
  return(detector);
}

/*******************************************************************************
 * mcuse_dir: set data/sim storage directory and create it,
 * or exit with error if exists
 ******************************************************************************/
static void
mcuse_dir(char *dir)
{
  if (!dir || !strlen(dir)) return;
#ifdef MC_PORTABLE
  fprintf(stderr, "Error: "
          "Directory output cannot be used with portable simulation (mcuse_dir)\n");
  exit(1);
#else  /* !MC_PORTABLE */
  /* handle file://directory URL type */
  if (strncmp(dir, "file://", strlen("file://")))
    mcdirname = dir;
  else
    mcdirname = dir+strlen("file://");
  
  
  
  MPI_MASTER(
    if(mkdir(mcdirname, 0777)) {
#ifndef DANSE
      fprintf(stderr, "Error: unable to create directory '%s' (mcuse_dir)\n", dir);
      fprintf(stderr, "(Maybe the directory already exists?)\n");
#endif
#ifdef USE_MPI
    MPI_Abort(MPI_COMM_WORLD, -1);
#endif
    exit(-1);
    }
  ); /* MPI_MASTER */
  
  /* remove trailing PATHSEP (if any) */
  while (strlen(mcdirname) && mcdirname[strlen(mcdirname) - 1] == MC_PATHSEP_C)
    mcdirname[strlen(mcdirname) - 1]='\0';
#endif /* !MC_PORTABLE */
} /* mcuse_dir */

/*******************************************************************************
* mcinfo: display instrument simulation info to stdout and exit
*******************************************************************************/
static void
mcinfo(void)
{
  fprintf(stdout, "begin instrument: %s\n", mcinstrument_name);
  mcinfo_out("  ", stdout);
  fprintf(stdout, "end instrument\n");
  fprintf(stdout, "begin simulation: %s\n", mcdirname ? mcdirname : ".");
  mcruninfo_out_backend("  ", stdout,1);
  fprintf(stdout, "end simulation\n");
  exit(0); /* includes MPI_Finalize in MPI mode */
} /* mcinfo */

#endif /* ndef MCCODE_R_IO_C */

/* end of the I/O section =================================================== */







/*******************************************************************************
* mcset_ncount: set total number of rays to generate
*******************************************************************************/
void mcset_ncount(unsigned long long int count)
{
  mcncount = count;
}

/* mcget_ncount: get total number of rays to generate */
unsigned long long int mcget_ncount(void)
{
  return mcncount;
}

/* mcget_run_num: get curent number of rays in TRACE */
unsigned long long int mcget_run_num(void)
{
  return mcrun_num;
}

/* mcsetn_arg: get ncount from a string argument */
static void
mcsetn_arg(char *arg)
{
  mcset_ncount((long long int) strtod(arg, NULL));
}

/* mcsetseed: set the random generator seed from a string argument */
static void
mcsetseed(char *arg)
{
  mcseed = atol(arg);
  if(mcseed) {
    srandom(mcseed);
  } else {
    fprintf(stderr, "Error: seed must not be zero (mcsetseed)\n");
    exit(1);
  }
}

/* Following part is only embedded when not redundent with mccode-r.h ========= */

#ifndef MCCODE_H

/* SECTION: MCDISPLAY support. =============================================== */

/*******************************************************************************
* Just output MCDISPLAY keywords to be caught by an external plotter client.
*******************************************************************************/

void mcdis_magnify(char *what){
  // Do nothing here, better use interactive zoom from the tools
}

void mcdis_line(double x1, double y1, double z1,
                double x2, double y2, double z2){
  printf("MCDISPLAY: multiline(2,%g,%g,%g,%g,%g,%g)\n",
         x1,y1,z1,x2,y2,z2);
}

void mcdis_dashed_line(double x1, double y1, double z1,
		       double x2, double y2, double z2, int n){
  int i;
  const double dx = (x2-x1)/(2*n+1);
  const double dy = (y2-y1)/(2*n+1);
  const double dz = (z2-z1)/(2*n+1);

  for(i = 0; i < n+1; i++)
    mcdis_line(x1 + 2*i*dx,     y1 + 2*i*dy,     z1 + 2*i*dz,
	       x1 + (2*i+1)*dx, y1 + (2*i+1)*dy, z1 + (2*i+1)*dz);
}

void mcdis_multiline(int count, ...){
  va_list ap;
  double x,y,z;

  printf("MCDISPLAY: multiline(%d", count);
  va_start(ap, count);
  while(count--)
    {
    x = va_arg(ap, double);
    y = va_arg(ap, double);
    z = va_arg(ap, double);
    printf(",%g,%g,%g", x, y, z);
    }
  va_end(ap);
  printf(")\n");
}

void mcdis_rectangle(char* plane, double x, double y, double z,
		     double width, double height){
  /* draws a rectangle in the plane           */
  /* x is ALWAYS width and y is ALWAYS height */
  if (strcmp("xy", plane)==0) {
    mcdis_multiline(5,
		    x - width/2, y - height/2, z,
		    x + width/2, y - height/2, z,
		    x + width/2, y + height/2, z,
		    x - width/2, y + height/2, z,
		    x - width/2, y - height/2, z);
  } else if (strcmp("xz", plane)==0) {
    mcdis_multiline(5,
		    x - width/2, y, z - height/2,
		    x + width/2, y, z - height/2,
		    x + width/2, y, z + height/2,
		    x - width/2, y, z + height/2,
		    x - width/2, y, z - height/2);
  } else if (strcmp("yz", plane)==0) {
    mcdis_multiline(5,
		    x, y - height/2, z - width/2,
		    x, y - height/2, z + width/2,
		    x, y + height/2, z + width/2,
		    x, y + height/2, z - width/2,
		    x, y - height/2, z - width/2);
  } else {

    fprintf(stderr, "Error: Definition of plane %s unknown\n", plane);
    exit(1);
  }
}

/*  draws a box with center at (x, y, z) and
    width (deltax), height (deltay), length (deltaz) */
void mcdis_box(double x, double y, double z,
	       double width, double height, double length){

  mcdis_rectangle("xy", x, y, z-length/2, width, height);
  mcdis_rectangle("xy", x, y, z+length/2, width, height);
  mcdis_line(x-width/2, y-height/2, z-length/2,
	     x-width/2, y-height/2, z+length/2);
  mcdis_line(x-width/2, y+height/2, z-length/2,
	     x-width/2, y+height/2, z+length/2);
  mcdis_line(x+width/2, y-height/2, z-length/2,
	     x+width/2, y-height/2, z+length/2);
  mcdis_line(x+width/2, y+height/2, z-length/2,
	     x+width/2, y+height/2, z+length/2);
}

void mcdis_circle(char *plane, double x, double y, double z, double r){
  printf("MCDISPLAY: circle('%s',%g,%g,%g,%g)\n", plane, x, y, z, r);
}

/* Draws a circle with center (x,y,z), radius (r), and in the plane
 * with normal (nx,ny,nz)*/
void mcdis_Circle(double x, double y, double z, double r, double nx, double ny, double nz){
    int i;
    if(nx==0 && ny && nz==0){
        for (i=0;i<24; i++){
            mcdis_line(x+r*sin(i*2*M_PI/24),y,z+r*cos(i*2*M_PI/24),
                    x+r*sin((i+1)*2*M_PI/24),y,z+r*cos((i+1)*2*M_PI/24));
        }
    }else{
        double mx,my,mz;
        /*generate perpendicular vector using (nx,ny,nz) and (0,1,0)*/
        vec_prod(mx,my,mz, 0,1,0, nx,ny,nz);
        NORM(mx,my,mz);
        /*draw circle*/
        for (i=0;i<24; i++){
            double ux,uy,uz;
            double wx,wy,wz;
            rotate(ux,uy,uz, mx,my,mz, i*2*M_PI/24, nx,ny,nz);
            rotate(wx,wy,wz, mx,my,mz, (i+1)*2*M_PI/24, nx,ny,nz);
            mcdis_line(x+ux*r,y+uy*r,z+uz*r,
                    x+wx*r,y+wy*r,z+wz*r);
        }
    }
}

/* Draws a cylinder with center at (x,y,z) with extent (r,height).
 * The cylinder axis is along the vector nx,ny,nz.
 * N determines how many vertical lines are drawn.*/
void mcdis_cylinder( double x, double y, double z,
        double r, double height, int N, double nx, double ny, double nz){
    int i;
    /*no lines make little sense - so trigger the default*/
    if(N<=0) N=5;

    NORM(nx,ny,nz);
    double h_2=height/2.0;
    mcdis_Circle(x+nx*h_2,y+ny*h_2,z+nz*h_2,r,nx,ny,nz);
    mcdis_Circle(x-nx*h_2,y-ny*h_2,z-nz*h_2,r,nx,ny,nz);

    double mx,my,mz;
    /*generate perpendicular vector using (nx,ny,nz) and (0,1,0)*/
    if(nx==0 && ny && nz==0){
        mx=my=0;mz=1;
    }else{
        vec_prod(mx,my,mz, 0,1,0, nx,ny,nz);
        NORM(mx,my,mz);
    }
    /*draw circle*/
    for (i=0; i<24; i++){
        double ux,uy,uz;
        rotate(ux,uy,uz, mx,my,mz, i*2*M_PI/24, nx,ny,nz);
        mcdis_line(x+nx*h_2+ux*r, y+ny*h_2+uy*r, z+nz*h_2+uz*r,
                 x-nx*h_2+ux*r, y-ny*h_2+uy*r, z-nz*h_2+uz*r);
    }
}

/* draws a sphere with center at (x,y,z) with extent (r)
 * The sphere is drawn using N longitudes and N latitudes.*/
void mcdis_sphere(double x, double y, double z, double r, int N){
    double nx,ny,nz;
    int i;
    /*no lines make little sense - so trigger the default*/
    if(N<=0) N=5;

    nx=0;ny=0;nz=1;
    mcdis_Circle(x,y,z,r,nx,ny,nz);
    for (i=1;i<N;i++){
        rotate(nx,ny,nz, nx,ny,nz, M_PI/N, 0,1,0);
        mcdis_Circle(x,y,z,r,nx,ny,nz);
    }
    /*lastly draw a great circle perpendicular to all N circles*/
    //mcdis_Circle(x,y,z,radius,1,0,0);

    for (i=1;i<=N;i++){
        double yy=-r+ 2*r*((double)i/(N+1));
        mcdis_Circle(x,y+yy ,z,  sqrt(r*r-yy*yy) ,0,1,0);
    }
}

/* SECTION: coordinates handling ============================================ */

/*******************************************************************************
* Since we use a lot of geometric calculations using Cartesian coordinates,
* we collect some useful routines here. However, it is also permissible to
* work directly on the underlying struct coords whenever that is most
* convenient (that is, the type Coords is not abstract).
*
* Coordinates are also used to store rotation angles around x/y/z axis.
*
* Since coordinates are used much like a basic type (such as double), the
* structure itself is passed and returned, rather than a pointer.
*
* At compile-time, the values of the coordinates may be unknown (for example
* a motor position). Hence coordinates are general expressions and not simple
* numbers. For this we used the type Coords_exp which has three CExp
* fields. For runtime (or calculations possible at compile time), we use
* Coords which contains three double fields.
*******************************************************************************/

/* coords_set: Assign coordinates. */
Coords
coords_set(MCNUM x, MCNUM y, MCNUM z)
{
  Coords a;

  a.x = x;
  a.y = y;
  a.z = z;
  return a;
}

/* coords_get: get coordinates. Required when 'x','y','z' are #defined as ray pars */
Coords
coords_get(Coords a, MCNUM *x, MCNUM *y, MCNUM *z)
{
  *x = a.x;
  *y = a.y;
  *z = a.z;
  return a;
}

/* coords_add: Add two coordinates. */
Coords
coords_add(Coords a, Coords b)
{
  Coords c;

  c.x = a.x + b.x;
  c.y = a.y + b.y;
  c.z = a.z + b.z;
  if (fabs(c.z) < 1e-14) c.z=0.0;
  return c;
}

/* coords_sub: Subtract two coordinates. */
Coords
coords_sub(Coords a, Coords b)
{
  Coords c;

  c.x = a.x - b.x;
  c.y = a.y - b.y;
  c.z = a.z - b.z;
  if (fabs(c.z) < 1e-14) c.z=0.0;
  return c;
}

/* coords_neg: Negate coordinates. */
Coords
coords_neg(Coords a)
{
  Coords b;

  b.x = -a.x;
  b.y = -a.y;
  b.z = -a.z;
  return b;
}

/* coords_scale: Scale a vector. */
Coords coords_scale(Coords b, double scale) {
  Coords a;

  a.x = b.x*scale;
  a.y = b.y*scale;
  a.z = b.z*scale;
  return a;
}

/* coords_sp: Scalar product: a . b */
double coords_sp(Coords a, Coords b) {
  double value;

  value = a.x*b.x + a.y*b.y + a.z*b.z;
  return value;
}

/* coords_xp: Cross product: a = b x c. */
Coords coords_xp(Coords b, Coords c) {
  Coords a;

  a.x = b.y*c.z - c.y*b.z;
  a.y = b.z*c.x - c.z*b.x;
  a.z = b.x*c.y - c.x*b.y;
  return a;
}

/* coords_len: Gives length of coords set. */
double coords_len(Coords a) {
  return sqrt(a.x*a.x + a.y*a.y + a.z*a.z);
}

/* coords_mirror: Mirror a in plane (through the origin) defined by normal n*/
Coords coords_mirror(Coords a, Coords n) {
  double t = scalar_prod(n.x, n.y, n.z, n.x, n.y, n.z);
  Coords b;
  if (t!=1) {
    t = sqrt(t);
    n.x /= t;
    n.y /= t;
    n.z /= t;
  }
  t=scalar_prod(a.x, a.y, a.z, n.x, n.y, n.z);
  b.x = a.x-2*t*n.x;
  b.y = a.y-2*t*n.y;
  b.z = a.z-2*t*n.z;
  return b;
}

/* coords_print: Print out vector values. */
void coords_print(Coords a) {

  fprintf(stdout, "(%f, %f, %f)\n", a.x, a.y, a.z);
  return;
}

mcstatic void coords_norm(Coords* c) {
	double temp = coords_sp(*c,*c);

	// Skip if we will end dividing by zero
	if (temp == 0) return;

	temp = sqrt(temp);

	c->x /= temp;
	c->y /= temp;
	c->z /= temp;
}

/*******************************************************************************
* The Rotation type implements a rotation transformation of a coordinate
* system in the form of a double[3][3] matrix.
*
* Contrary to the Coords type in coords.c, rotations are passed by
* reference. Functions that yield new rotations do so by writing to an
* explicit result parameter; rotations are not returned from functions. The
* reason for this is that arrays cannot by returned from functions (though
* structures can; thus an alternative would have been to wrap the
* double[3][3] array up in a struct). Such are the ways of C programming.
*
* A rotation represents the tranformation of the coordinates of a vector when
* changing between coordinate systems that are rotated with respect to each
* other. For example, suppose that coordinate system Q is rotated 45 degrees
* around the Z axis with respect to coordinate system P. Let T be the
* rotation transformation representing a 45 degree rotation around Z. Then to
* get the coordinates of a vector r in system Q, apply T to the coordinates
* of r in P. If r=(1,0,0) in P, it will be (sqrt(1/2),-sqrt(1/2),0) in
* Q. Thus we should be careful when interpreting the sign of rotation angles:
* they represent the rotation of the coordinate systems, not of the
* coordinates (which has opposite sign).
*******************************************************************************/

/*******************************************************************************
* rot_set_rotation: Get transformation for rotation first phx around x axis,
* then phy around y, then phz around z.
*******************************************************************************/
void
rot_set_rotation(Rotation t, double phx, double phy, double phz)
{
  if ((phx == 0) && (phy == 0) && (phz == 0)) {
    t[0][0] = 1.0;
    t[0][1] = 0.0;
    t[0][2] = 0.0;
    t[1][0] = 0.0;
    t[1][1] = 1.0;
    t[1][2] = 0.0;
    t[2][0] = 0.0;
    t[2][1] = 0.0;
    t[2][2] = 1.0;
  } else {
    double cx = cos(phx);
    double sx = sin(phx);
    double cy = cos(phy);
    double sy = sin(phy);
    double cz = cos(phz);
    double sz = sin(phz);

    t[0][0] = cy*cz;
    t[0][1] = sx*sy*cz + cx*sz;
    t[0][2] = sx*sz - cx*sy*cz;
    t[1][0] = -cy*sz;
    t[1][1] = cx*cz - sx*sy*sz;
    t[1][2] = sx*cz + cx*sy*sz;
    t[2][0] = sy;
    t[2][1] = -sx*cy;
    t[2][2] = cx*cy;
  }
}

/*******************************************************************************
* rot_test_identity: Test if rotation is identity
*******************************************************************************/
int
rot_test_identity(Rotation t)
{
  return (t[0][0] + t[1][1] + t[2][2] == 3);
}

/*******************************************************************************
* rot_mul: Matrix multiplication of transformations (this corresponds to
* combining transformations). After rot_mul(T1, T2, T3), doing T3 is
* equal to doing first T2, then T1.
* Note that T3 must not alias (use the same array as) T1 or T2.
*******************************************************************************/
void
rot_mul(Rotation t1, Rotation t2, Rotation t3)
{
  if (rot_test_identity(t1)) {
    rot_copy(t3, t2);
  } else if (rot_test_identity(t2)) {
    rot_copy(t3, t1);
  } else {
    int i,j;
    for(i = 0; i < 3; i++)
      for(j = 0; j < 3; j++)
	t3[i][j] = t1[i][0]*t2[0][j] + t1[i][1]*t2[1][j] + t1[i][2]*t2[2][j];
  }
}

/*******************************************************************************
* rot_copy: Copy a rotation transformation (arrays cannot be assigned in C).
*******************************************************************************/
void
rot_copy(Rotation dest, Rotation src)
{
  int i,j;
  for(i = 0; i < 3; i++)
    for(j = 0; j < 3; j++)
      dest[i][j] = src[i][j];
}

/*******************************************************************************
* rot_transpose: Matrix transposition, which is inversion for Rotation matrices
*******************************************************************************/
void
rot_transpose(Rotation src, Rotation dst)
{
  dst[0][0] = src[0][0];
  dst[0][1] = src[1][0];
  dst[0][2] = src[2][0];
  dst[1][0] = src[0][1];
  dst[1][1] = src[1][1];
  dst[1][2] = src[2][1];
  dst[2][0] = src[0][2];
  dst[2][1] = src[1][2];
  dst[2][2] = src[2][2];
}

/*******************************************************************************
* rot_apply: returns t*a
*******************************************************************************/
Coords
rot_apply(Rotation t, Coords a)
{
  Coords b;
  if (rot_test_identity(t)) {
    return a;
  } else {
    b.x = t[0][0]*a.x + t[0][1]*a.y + t[0][2]*a.z;
    b.y = t[1][0]*a.x + t[1][1]*a.y + t[1][2]*a.z;
    b.z = t[2][0]*a.x + t[2][1]*a.y + t[2][2]*a.z;
    return b;
  }
}

/**
 * Pretty-printing of rotation matrices.
 */
void rot_print(Rotation rot) {
	printf("[ %4.2f %4.2f %4.2f ]\n",
			rot[0][0], rot[0][1], rot[0][2]);
	printf("[ %4.2f %4.2f %4.2f ]\n",
			rot[1][0], rot[1][1], rot[1][2]);
	printf("[ %4.2f %4.2f %4.2f ]\n\n",
			rot[2][0], rot[2][1], rot[2][2]);
}

/**
 * Vector product: used by vec_prod (mccode-r.h). Use coords_xp for Coords.
 */
mcstatic void vec_prod_func(double *x, double *y, double *z,
		double x1, double y1, double z1,
		double x2, double y2, double z2) {
    *x = (y1)*(z2) - (y2)*(z1);
    *y = (z1)*(x2) - (z2)*(x1);
    *z = (x1)*(y2) - (x2)*(y1);
}

/**
 * Scalar product: use coords_sp for Coords.
 */
mcstatic double scalar_prod(
		double x1, double y1, double z1,
		double x2, double y2, double z2) {
	return ((x1 * x2) + (y1 * y2) + (z1 * z2));
}

/*******************************************************************************
* mccoordschange: applies rotation to (x y z) and (vx vy vz) and Spin (sx,sy,sz)
*******************************************************************************/
void
mccoordschange(Coords a, Rotation t, double *x, double *y, double *z,
               double *vx, double *vy, double *vz, double *sx, double *sy, double *sz)
{
  Coords b, c;

  b.x = *x;
  b.y = *y;
  b.z = *z;
  c = rot_apply(t, b);
  b = coords_add(c, a);
  *x = b.x;
  *y = b.y;
  *z = b.z;

  if ( (vz && vy  && vx) && (*vz != 0.0 || *vx != 0.0 || *vy != 0.0) ) mccoordschange_polarisation(t, vx, vy, vz);

  if ( (sz && sy  && sx) && (*sz != 0.0 || *sx != 0.0 || *sy != 0.0) ) mccoordschange_polarisation(t, sx, sy, sz);

}

/*******************************************************************************
* mccoordschange_polarisation: applies rotation to vector (sx sy sz)
*******************************************************************************/
void
mccoordschange_polarisation(Rotation t, double *sx, double *sy, double *sz)
{
  Coords b, c;

  b.x = *sx;
  b.y = *sy;
  b.z = *sz;
  c = rot_apply(t, b);
  *sx = c.x;
  *sy = c.y;
  *sz = c.z;
}

/* SECTION: vector math  ==================================================== */

/* normal_vec_func: Compute normal vector to (x,y,z). */
mcstatic void normal_vec_func(double *nx, double *ny, double *nz,
                double x, double y, double z)
{
  double ax = fabs(x);
  double ay = fabs(y);
  double az = fabs(z);
  double l;
  if(x == 0 && y == 0 && z == 0)
  {
    *nx = 0;
    *ny = 0;
    *nz = 0;
    return;
  }
  if(ax < ay)
  {
    if(ax < az)
    {                           /* Use X axis */
      l = sqrt(z*z + y*y);
      *nx = 0;
      *ny = z/l;
      *nz = -y/l;
      return;
    }
  }
  else
  {
    if(ay < az)
    {                           /* Use Y axis */
      l = sqrt(z*z + x*x);
      *nx = z/l;
      *ny = 0;
      *nz = -x/l;
      return;
    }
  }
  /* Use Z axis */
  l = sqrt(y*y + x*x);
  *nx = y/l;
  *ny = -x/l;
  *nz = 0;
} /* normal_vec */

/*******************************************************************************
 * solve_2nd_order: second order equation solve: A*t^2 + B*t + C = 0
 * solve_2nd_order(&t1, NULL, A,B,C)
 *   returns 0 if no solution was found, or set 't1' to the smallest positive
 *   solution.
 * solve_2nd_order(&t1, &t2, A,B,C)
 *   same as with &t2=NULL, but also returns the second solution.
 * EXAMPLE usage for intersection of a trajectory with a plane in gravitation
 * field (gx,gy,gz):
 * The neutron starts at point r=(x,y,z) with velocityv=(vx vy vz). The plane
 * has a normal vector n=(nx,ny,nz) and contains the point W=(wx,wy,wz).
 * The problem consists in solving the 2nd order equation:
 *      1/2.n.g.t^2 + n.v.t + n.(r-W) = 0
 * so that A = 0.5 n.g; B = n.v; C = n.(r-W);
 * Without acceleration, t=-n.(r-W)/n.v
 ******************************************************************************/
int solve_2nd_order(double *t1, double *t2,
                  double A,  double B,  double C)
{
  int ret=0;

  if (!t1) return 0;
  *t1 = 0;
  if (t2) *t2=0;

  if (fabs(A) < 1E-10) /* approximate to linear equation: A ~ 0 */
  {
    if (B) {  *t1 = -C/B; ret=1; if (t2) *t2=*t1; }
    /* else no intersection: A=B=0 ret=0 */
  }
  else
  {
    double D;
    D = B*B - 4*A*C;
    if (D >= 0) /* Delta > 0: two solutions */
    {
      double sD, dt1, dt2;
      sD = sqrt(D);
      dt1 = (-B + sD)/2/A;
      dt2 = (-B - sD)/2/A;
      /* we identify very small values with zero */
      if (fabs(dt1) < 1e-10) dt1=0.0;
      if (fabs(dt2) < 1e-10) dt2=0.0;

      /* now we choose the smallest positive solution */
      if      (dt1<=0.0 && dt2>0.0) ret=2; /* dt2 positive */
      else if (dt2<=0.0 && dt1>0.0) ret=1; /* dt1 positive */
      else if (dt1> 0.0 && dt2>0.0)
      {  if (dt1 < dt2) ret=1; else ret=2; } /* all positive: min(dt1,dt2) */
      /* else two solutions are negative. ret=-1 */
      if (ret==1) { *t1 = dt1;  if (t2) *t2=dt2; }
      else        { *t1 = dt2;  if (t2) *t2=dt1; }
      ret=2;  /* found 2 solutions and t1 is the positive one */
    } /* else Delta <0: no intersection. ret=0 */
  }
  return(ret);
} /* solve_2nd_order */

/*******************************************************************************
 * randvec_target_circle: Choose random direction towards target at (x,y,z)
 * with given radius.
 * If radius is zero, choose random direction in full 4PI, no target.
 ******************************************************************************/
void
randvec_target_circle(double *xo, double *yo, double *zo, double *solid_angle,
               double xi, double yi, double zi, double radius)
{
  double l2, phi, theta, nx, ny, nz, xt, yt, zt, xu, yu, zu;

  if(radius == 0.0)
  {
    /* No target, choose uniformly a direction in full 4PI solid angle. */
    theta = acos (1 - rand0max(2));
    phi = rand0max(2 * PI);
    if(solid_angle)
      *solid_angle = 4*PI;
    nx = 1;
    ny = 0;
    nz = 0;
    yi = sqrt(xi*xi+yi*yi+zi*zi);
    zi = 0;
    xi = 0;
  }
  else
  {
    double costheta0;
    l2 = xi*xi + yi*yi + zi*zi; /* sqr Distance to target. */
    costheta0 = sqrt(l2/(radius*radius+l2));
    if (radius < 0) costheta0 *= -1;
    if(solid_angle)
    {
      /* Compute solid angle of target as seen from origin. */
        *solid_angle = 2*PI*(1 - costheta0);
    }

    /* Now choose point uniformly on circle surface within angle theta0 */
    theta = acos (1 - rand0max(1 - costheta0)); /* radius on circle */
    phi = rand0max(2 * PI); /* rotation on circle at given radius */
    /* Now, to obtain the desired vector rotate (xi,yi,zi) angle theta around a
       perpendicular axis u=i x n and then angle phi around i. */
    if(xi == 0 && zi == 0)
    {
      nx = 1;
      ny = 0;
      nz = 0;
    }
    else
    {
      nx = -zi;
      nz = xi;
      ny = 0;
    }
  }

  /* [xyz]u = [xyz]i x n[xyz] (usually vertical) */
  vec_prod(xu,  yu,  zu, xi, yi, zi,        nx, ny, nz);
  /* [xyz]t = [xyz]i rotated theta around [xyz]u */
  rotate  (xt,  yt,  zt, xi, yi, zi, theta, xu, yu, zu);
  /* [xyz]o = [xyz]t rotated phi around n[xyz] */
  rotate (*xo, *yo, *zo, xt, yt, zt, phi, xi, yi, zi);
} /* randvec_target_circle */

/*******************************************************************************
 * randvec_target_rect_angular: Choose random direction towards target at
 * (xi,yi,zi) with given ANGULAR dimension height x width. height=phi_x=[0,PI],
 * width=phi_y=[0,2*PI] (radians)
 * If height or width is zero, choose random direction in full 4PI, no target.
 *******************************************************************************/
void
randvec_target_rect_angular(double *xo, double *yo, double *zo, double *solid_angle,
               double xi, double yi, double zi, double width, double height, Rotation A)
{
  double theta, phi, nx, ny, nz, xt, yt, zt, xu, yu, zu;
  Coords tmp;
  Rotation Ainverse;

  rot_transpose(A, Ainverse);

  if(height == 0.0 || width == 0.0)
  {
    randvec_target_circle(xo, yo, zo, solid_angle,
               xi, yi, zi, 0);
    return;
  }
  else
  {
    if(solid_angle)
    {
      /* Compute solid angle of target as seen from origin. */
      *solid_angle = 2*fabs(width*sin(height/2));
    }

    /* Go to global coordinate system */

    tmp = coords_set(xi, yi, zi);
    tmp = rot_apply(Ainverse, tmp);
    coords_get(tmp, &xi, &yi, &zi);

    /* Now choose point uniformly on the unit sphere segment with angle theta/phi */
    phi   = width*randpm1()/2.0;
    theta = asin(randpm1()*sin(height/2.0));
    /* Now, to obtain the desired vector rotate (xi,yi,zi) angle theta around
       n, and then phi around u. */
    if(xi == 0 && zi == 0)
    {
      nx = 1;
      ny = 0;
      nz = 0;
    }
    else
    {
      nx = -zi;
      nz = xi;
      ny = 0;
    }
  }

  /* [xyz]u = [xyz]i x n[xyz] (usually vertical) */
  vec_prod(xu,  yu,  zu, xi, yi, zi,        nx, ny, nz);
  /* [xyz]t = [xyz]i rotated theta around [xyz]u */
  rotate  (xt,  yt,  zt, xi, yi, zi, theta, nx, ny, nz);
  /* [xyz]o = [xyz]t rotated phi around n[xyz] */
  rotate (*xo, *yo, *zo, xt, yt, zt, phi, xu,  yu,  zu);

  /* Go back to local coordinate system */
  tmp = coords_set(*xo, *yo, *zo);
  tmp = rot_apply(A, tmp);
  coords_get(tmp, &*xo, &*yo, &*zo);

} /* randvec_target_rect_angular */

/*******************************************************************************
 * randvec_target_rect_real: Choose random direction towards target at (xi,yi,zi)
 * with given dimension height x width (in meters !).
 *
 * Local emission coordinate is taken into account and corrected for 'order' times.
 * (See remarks posted to mcstas-users by George Apostolopoulus <gapost@ipta.demokritos.gr>)
 *
 * If height or width is zero, choose random direction in full 4PI, no target.
 *
 * Traditionally, this routine had the name randvec_target_rect - this is now a
 * a define (see mcstas-r.h) pointing here. If you use the old rouine, you are NOT
 * taking the local emmission coordinate into account.
*******************************************************************************/

void
randvec_target_rect_real(double *xo, double *yo, double *zo, double *solid_angle,
               double xi, double yi, double zi,
               double width, double height, Rotation A,
               double lx, double ly, double lz, int order)
{
  double dx, dy, dist, dist_p, nx, ny, nz, mx, my, mz, n_norm, m_norm;
  double cos_theta;
  Coords tmp;
  Rotation Ainverse;

  rot_transpose(A, Ainverse);

  if(height == 0.0 || width == 0.0)
  {
    randvec_target_circle(xo, yo, zo, solid_angle,
               xi, yi, zi, 0);
    return;
  }
  else
  {

    /* Now choose point uniformly on rectangle within width x height */
    dx = width*randpm1()/2.0;
    dy = height*randpm1()/2.0;

    /* Determine distance to target plane*/
    dist = sqrt(xi*xi + yi*yi + zi*zi);
    /* Go to global coordinate system */

    tmp = coords_set(xi, yi, zi);
    tmp = rot_apply(Ainverse, tmp);
    coords_get(tmp, &xi, &yi, &zi);

    /* Determine vector normal to trajectory axis (z) and gravity [0 1 0] */
    vec_prod(nx, ny, nz, xi, yi, zi, 0, 1, 0);

    /* This now defines the x-axis, normalize: */
    n_norm=sqrt(nx*nx + ny*ny + nz*nz);
    nx = nx/n_norm;
    ny = ny/n_norm;
    nz = nz/n_norm;

    /* Now, determine our y-axis (vertical in many cases...) */
    vec_prod(mx, my, mz, xi, yi, zi, nx, ny, nz);
    m_norm=sqrt(mx*mx + my*my + mz*mz);
    mx = mx/m_norm;
    my = my/m_norm;
    mz = mz/m_norm;

    /* Our output, random vector can now be defined by linear combination: */

    *xo = xi + dx * nx + dy * mx;
    *yo = yi + dx * ny + dy * my;
    *zo = zi + dx * nz + dy * mz;

    /* Go back to local coordinate system */
    tmp = coords_set(*xo, *yo, *zo);
    tmp = rot_apply(A, tmp);
    coords_get(tmp, &*xo, &*yo, &*zo);

    /* Go back to local coordinate system */
    tmp = coords_set(xi, yi, zi);
    tmp = rot_apply(A, tmp);
    coords_get(tmp, &xi, &yi, &zi);

    if (solid_angle) {
      /* Calculate vector from local point to remote random point */
      lx = *xo - lx;
      ly = *yo - ly;
      lz = *zo - lz;
      dist_p = sqrt(lx*lx + ly*ly + lz*lz);

      /* Adjust the 'solid angle' */
      /* 1/r^2 to the chosen point times cos(\theta) between the normal */
      /* vector of the target rectangle and direction vector of the chosen point. */
      cos_theta = (xi * lx + yi * ly + zi * lz) / (dist * dist_p);
      *solid_angle = width * height / (dist_p * dist_p);
      int counter;
      for (counter = 0; counter < order; counter++) {
	*solid_angle = *solid_angle * cos_theta;
      }
    }
  }
} /* randvec_target_rect_real */

/* SECTION: random numbers ================================================== */

/*
 * Copyright (c) 1983 Regents of the University of California.
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms are permitted
 * provided that the above copyright notice and this paragraph are
 * duplicated in all such forms and that any documentation,
 * advertising materials, and other materials related to such
 * distribution and use acknowledge that the software was developed
 * by the University of California, Berkeley.  The name of the
 * University may not be used to endorse or promote products derived
 * from this software without specific prior written permission.
 * THIS SOFTWARE IS PROVIDED ``AS IS'' AND WITHOUT ANY EXPRESS OR
 * IMPLIED WARRANTIES, INCLUDING, WITHOUT LIMITATION, THE IMPLIED
 * WARRANTIES OF MERCHANTIBILITY AND FITNESS FOR A PARTICULAR PURPOSE.
 */

/*
 * This is derived from the Berkeley source:
 *        @(#)random.c        5.5 (Berkeley) 7/6/88
 * It was reworked for the GNU C Library by Roland McGrath.
 * Rewritten to use reentrant functions by Ulrich Drepper, 1995.
 */

/*******************************************************************************
* Modified for McStas from glibc 2.0.7pre1 stdlib/random.c and
* stdlib/random_r.c.
*
* This way random() is more than four times faster compared to calling
* standard glibc random() on ix86 Linux, probably due to multithread support,
* ELF shared library overhead, etc. It also makes McStas generated
* simulations more portable (more likely to behave identically across
* platforms, important for parrallel computations).
*******************************************************************************/


#define        TYPE_3                3
#define        BREAK_3                128
#define        DEG_3                31
#define        SEP_3                3

static mc_int32_t randtbl[DEG_3 + 1] =
  {
    TYPE_3,

    -1726662223, 379960547, 1735697613, 1040273694, 1313901226,
    1627687941, -179304937, -2073333483, 1780058412, -1989503057,
    -615974602, 344556628, 939512070, -1249116260, 1507946756,
    -812545463, 154635395, 1388815473, -1926676823, 525320961,
    -1009028674, 968117788, -123449607, 1284210865, 435012392,
    -2017506339, -911064859, -370259173, 1132637927, 1398500161,
    -205601318,
  };

static mc_int32_t *fptr = &randtbl[SEP_3 + 1];
static mc_int32_t *rptr = &randtbl[1];
static mc_int32_t *state = &randtbl[1];
#define rand_deg DEG_3
#define rand_sep SEP_3
static mc_int32_t *end_ptr = &randtbl[sizeof (randtbl) / sizeof (randtbl[0])];

mc_int32_t
mc_random (void)
{
  mc_int32_t result;

  *fptr += *rptr;
  /* Chucking least random bit.  */
  result = (*fptr >> 1) & 0x7fffffff;
  ++fptr;
  if (fptr >= end_ptr)
  {
    fptr = state;
    ++rptr;
  }
  else
  {
    ++rptr;
    if (rptr >= end_ptr)
      rptr = state;
  }
  return result;
}

void
mc_srandom (unsigned int x)
{
  /* We must make sure the seed is not 0.  Take arbitrarily 1 in this case.  */
  state[0] = x ? x : 1;
  {
    long int i;
    for (i = 1; i < rand_deg; ++i)
    {
      /* This does:
         state[i] = (16807 * state[i - 1]) % 2147483647;
         but avoids overflowing 31 bits.  */
      long int hi = state[i - 1] / 127773;
      long int lo = state[i - 1] % 127773;
      long int test = 16807 * lo - 2836 * hi;
      state[i] = test + (test < 0 ? 2147483647 : 0);
    }
    fptr = &state[rand_sep];
    rptr = &state[0];
    for (i = 0; i < 10 * rand_deg; ++i)
      random ();
  }
}

/* "Mersenne Twister", by Makoto Matsumoto and Takuji Nishimura. */
/* See http://www.math.keio.ac.jp/~matumoto/emt.html for original source. */


/*
   A C-program for MT19937, with initialization improved 2002/1/26.
   Coded by Takuji Nishimura and Makoto Matsumoto.

   Before using, initialize the state by using mt_srandom(seed)
   or init_by_array(init_key, key_length).

   Copyright (C) 1997 - 2002, Makoto Matsumoto and Takuji Nishimura,
   All rights reserved.

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions
   are met:

     1. Redistributions of source code must retain the above copyright
        notice, this list of conditions and the following disclaimer.

     2. Redistributions in binary form must reproduce the above copyright
        notice, this list of conditions and the following disclaimer in the
        documentation and/or other materials provided with the distribution.

     3. The names of its contributors may not be used to endorse or promote
        products derived from this software without specific prior written
        permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
   A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
   CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
   EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
   PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
   PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
   LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


   Any feedback is very welcome.
   http://www.math.keio.ac.jp/matumoto/emt.html
   email: matumoto@math.keio.ac.jp
*/

#include <stdio.h>

/* Period parameters */
#define N 624
#define M 397
#define MATRIX_A 0x9908b0dfUL   /* constant vector a */
#define UPPER_MASK 0x80000000UL /* most significant w-r bits */
#define LOWER_MASK 0x7fffffffUL /* least significant r bits */

static unsigned long mt[N]; /* the array for the state vector  */
static int mti=N+1; /* mti==N+1 means mt[N] is not initialized */

/* initializes mt[N] with a seed */
void mt_srandom(unsigned long s)
{
    mt[0]= s & 0xffffffffUL;
    for (mti=1; mti<N; mti++) {
        mt[mti] =
            (1812433253UL * (mt[mti-1] ^ (mt[mti-1] >> 30)) + mti);
        /* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
        /* In the previous versions, MSBs of the seed affect   */
        /* only MSBs of the array mt[].                        */
        /* 2002/01/09 modified by Makoto Matsumoto             */
        mt[mti] &= 0xffffffffUL;
        /* for >32 bit machines */
    }
}

/* initialize by an array with array-length */
/* init_key is the array for initializing keys */
/* key_length is its length */
void init_by_array(unsigned long init_key[], unsigned long key_length)
{
    int i, j, k;
    mt_srandom(19650218UL);
    i=1; j=0;
    k = (N>key_length ? N : key_length);
    for (; k; k--) {
        mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1664525UL))
          + init_key[j] + j; /* non linear */
        mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
        i++; j++;
        if (i>=N) { mt[0] = mt[N-1]; i=1; }
        if (j>=key_length) j=0;
    }
    for (k=N-1; k; k--) {
        mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1566083941UL))
          - i; /* non linear */
        mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
        i++;
        if (i>=N) { mt[0] = mt[N-1]; i=1; }
    }

    mt[0] = 0x80000000UL; /* MSB is 1; assuring non-zero initial array */
}

/* generates a random number on [0,0xffffffff]-interval */
unsigned long mt_random(void)
{
    unsigned long y;
    static unsigned long mag01[2]={0x0UL, MATRIX_A};
    /* mag01[x] = x * MATRIX_A  for x=0,1 */

    if (mti >= N) { /* generate N words at one time */
        int kk;

        if (mti == N+1)   /* if mt_srandom() has not been called, */
            mt_srandom(5489UL); /* a default initial seed is used */

        for (kk=0;kk<N-M;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+M] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        for (;kk<N-1;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+(M-N)] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        y = (mt[N-1]&UPPER_MASK)|(mt[0]&LOWER_MASK);
        mt[N-1] = mt[M-1] ^ (y >> 1) ^ mag01[y & 0x1UL];

        mti = 0;
    }

    y = mt[mti++];

    /* Tempering */
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);

    return y;
}

#undef N
#undef M
#undef MATRIX_A
#undef UPPER_MASK
#undef LOWER_MASK

/* End of "Mersenne Twister". */

/* End of McCode random number routine. */

/* randnorm: generate a random number from normal law */
double
randnorm(void)
{
  static double v1, v2, s;
  static int phase = 0;
  double X, u1, u2;

  if(phase == 0)
  {
    do
    {
      u1 = rand01();
      u2 = rand01();
      v1 = 2*u1 - 1;
      v2 = 2*u2 - 1;
      s = v1*v1 + v2*v2;
    } while(s >= 1 || s == 0);

    X = v1*sqrt(-2*log(s)/s);
  }
  else
  {
    X = v2*sqrt(-2*log(s)/s);
  }

  phase = 1 - phase;
  return X;
}

/**
 * Generate a random number from -1 to 1 with triangle distribution
 */
double randtriangle(void) {
	double randnum = rand01();
	if (randnum>0.5) return(1-sqrt(2*(randnum-0.5)));
	else return(sqrt(2*randnum)-1);
}

/**
 * Random number between 0.0 and 1.0 (including?)
 */
double rand01() {
	double randnum;
	randnum = (double) random();
	randnum /= (double) MC_RAND_MAX + 1;
	return randnum;
}

/**
 * Return a random number between 1 and -1
 */
double randpm1() {
	double randnum;
	randnum = (double) random();
	randnum /= ((double) MC_RAND_MAX + 1) / 2;
	randnum -= 1;
	return randnum;
}

/**
 * Return a random number between 0 and max.
 */
double rand0max(double max) {
	double randnum;
	randnum = (double) random();
	randnum /= ((double) MC_RAND_MAX + 1) / max;
	return randnum;
}

/**
 * Return a random number between min and max.
 */
double randminmax(double min, double max) {
	return rand0max(max - min) + max;
}

/* SECTION: main and signal handlers ======================================== */

/*******************************************************************************
* mchelp: displays instrument executable help with possible options
*******************************************************************************/
static void
mchelp(char *pgmname)
{
  int i;

  fprintf(stderr, "%s (%s) instrument simulation, generated with " MCCODE_STRING " (" MCCODE_DATE ")\n", mcinstrument_name, mcinstrument_source);
  fprintf(stderr, "Usage: %s [options] [parm=value ...]\n", pgmname);
  fprintf(stderr,
"Options are:\n"
"  -s SEED   --seed=SEED      Set random seed (must be != 0)\n"
"  -n COUNT  --ncount=COUNT   Set number of " MCCODE_PARTICLE "s to simulate.\n"
"  -d DIR    --dir=DIR        Put all data files in directory DIR.\n"
"  -t        --trace          Enable trace of " MCCODE_PARTICLE "s through instrument.\n"
"  -g        --gravitation    Enable gravitation for all trajectories.\n"
"  --no-output-files          Do not write any data files.\n"
"  -h        --help           Show this help message.\n"
"  -i        --info           Detailed instrument information.\n"
"  --format=FORMAT            Output data files using FORMAT="
   FLAVOR_UPPER
#ifdef USE_NEXUS
   " NEXUS"
#endif
"\n\n"
);
#ifdef USE_MPI
  fprintf(stderr,
  "This instrument has been compiled with MPI support.\n  Use 'mpirun %s [options] [parm=value ...]'.\n", pgmname);
#endif
  if(mcnumipar > 0)
  {
    fprintf(stderr, "Instrument parameters are:\n");
    for(i = 0; i < mcnumipar; i++)
      if (mcinputtable[i].val && strlen(mcinputtable[i].val))
        fprintf(stderr, "  %-16s(%s) [default='%s']\n", mcinputtable[i].name,
        (*mcinputtypes[mcinputtable[i].type].parminfo)(mcinputtable[i].name),
        mcinputtable[i].val);
      else
        fprintf(stderr, "  %-16s(%s)\n", mcinputtable[i].name,
        (*mcinputtypes[mcinputtable[i].type].parminfo)(mcinputtable[i].name));
  }

#ifndef NOSIGNALS
  fprintf(stderr, "Known signals are: "
#ifdef SIGUSR1
  "USR1 (status) "
#endif
#ifdef SIGUSR2
  "USR2 (save) "
#endif
#ifdef SIGBREAK
  "BREAK (save) "
#endif
#ifdef SIGTERM
  "TERM (save and exit)"
#endif
  "\n");
#endif /* !NOSIGNALS */
} /* mchelp */


/* mcshowhelp: show help and exit with 0 */
static void
mcshowhelp(char *pgmname)
{
  mchelp(pgmname);
  exit(0);
}

/* mcusage: display usage when error in input arguments and exit with 1 */
static void
mcusage(char *pgmname)
{
  fprintf(stderr, "Error: incorrect command line arguments\n");
  mchelp(pgmname);
  exit(1);
}

/* mcenabletrace: enable trace/mcdisplay or error if requires recompile */
static void
mcenabletrace(void)
{
 if(mctraceenabled)
  mcdotrace = 1;
 else
 {
   fprintf(stderr,
           "Error: trace not enabled (mcenabletrace)\n"
           "Please re-run the " MCCODE_NAME " compiler "
                   "with the --trace option, or rerun the\n"
           "C compiler with the MC_TRACE_ENABLED macro defined.\n");
   exit(1);
 }
}

/*******************************************************************************
* mcreadparams: request parameters from the prompt (or use default)
*******************************************************************************/
void
mcreadparams(void)
{
  int i,j,status;
  static char buf[CHAR_BUF_LENGTH];
  char *p;
  int len;

  MPI_MASTER(printf("Instrument parameters for %s (%s)\n",
                    mcinstrument_name, mcinstrument_source));

  for(i = 0; mcinputtable[i].name != 0; i++)
  {
    do
    {
      MPI_MASTER(
                 if (mcinputtable[i].val && strlen(mcinputtable[i].val))
                   printf("Set value of instrument parameter %s (%s) [default='%s']:\n",
                          mcinputtable[i].name,
                          (*mcinputtypes[mcinputtable[i].type].parminfo)
                          (mcinputtable[i].name), mcinputtable[i].val);
                 else
                   printf("Set value of instrument parameter %s (%s):\n",
                          mcinputtable[i].name,
                          (*mcinputtypes[mcinputtable[i].type].parminfo)
                          (mcinputtable[i].name));
                 fflush(stdout);
                 );
#ifdef USE_MPI
      if(mpi_node_rank == mpi_node_root)
        {
          p = fgets(buf, CHAR_BUF_LENGTH, stdin);
          if(p == NULL)
            {
              fprintf(stderr, "Error: empty input for paramater %s (mcreadparams)\n", mcinputtable[i].name);
              exit(1);
            }
        }
      else
        p = buf;
      MPI_Bcast(buf, CHAR_BUF_LENGTH, MPI_CHAR, mpi_node_root, MPI_COMM_WORLD);
#else /* !USE_MPI */
      p = fgets(buf, CHAR_BUF_LENGTH, stdin);
      if(p == NULL)
        {
          fprintf(stderr, "Error: empty input for paramater %s (mcreadparams)\n", mcinputtable[i].name);
          exit(1);
        }
#endif /* USE_MPI */
      len = strlen(buf);
      if (!len || (len == 1 && (buf[0] == '\n' || buf[0] == '\r')))
      {
        if (mcinputtable[i].val && strlen(mcinputtable[i].val)) {
          strncpy(buf, mcinputtable[i].val, CHAR_BUF_LENGTH);  /* use default value */
          len = strlen(buf);
        }
      }
      for(j = 0; j < 2; j++)
      {
        if(len > 0 && (buf[len - 1] == '\n' || buf[len - 1] == '\r'))
        {
          len--;
          buf[len] = '\0';
        }
      }

      status = (*mcinputtypes[mcinputtable[i].type].getparm)
                   (buf, mcinputtable[i].par);
      if(!status)
      {
        (*mcinputtypes[mcinputtable[i].type].error)(mcinputtable[i].name, buf);
        if (!mcinputtable[i].val || strlen(mcinputtable[i].val)) {
          fprintf(stderr, "       Change %s default value in instrument definition.\n", mcinputtable[i].name);
          exit(1);
        }
      }
    } while(!status);
  }
} /* mcreadparams */

/*******************************************************************************
* mcparseoptions: parse command line arguments (options, parameters)
*******************************************************************************/
void
mcparseoptions(int argc, char *argv[])
{
  int i, j;
  char *p;
  int paramset = 0, *paramsetarray;
  char *usedir=NULL;

  /* Add one to mcnumipar to avoid allocating zero size memory block. */
  paramsetarray = (int*)malloc((mcnumipar + 1)*sizeof(*paramsetarray));
  if(paramsetarray == NULL)
  {
    fprintf(stderr, "Error: insufficient memory (mcparseoptions)\n");
    exit(1);
  }
  for(j = 0; j < mcnumipar; j++)
    {
      paramsetarray[j] = 0;
      if (mcinputtable[j].val != NULL && strlen(mcinputtable[j].val))
      {
        int  status;
        char buf[CHAR_BUF_LENGTH];
        strncpy(buf, mcinputtable[j].val, CHAR_BUF_LENGTH);
        status = (*mcinputtypes[mcinputtable[j].type].getparm)
                   (buf, mcinputtable[j].par);
        if(!status) fprintf(stderr, "Invalid '%s' default value %s in instrument definition (mcparseoptions)\n", mcinputtable[j].name, buf);
        else paramsetarray[j] = 1;
      } else {
        (*mcinputtypes[mcinputtable[j].type].getparm)
          (NULL, mcinputtable[j].par);
        paramsetarray[j] = 0;
      }
    }
  for(i = 1; i < argc; i++)
  {
    if(!strcmp("-s", argv[i]) && (i + 1) < argc)
      mcsetseed(argv[++i]);
    else if(!strncmp("-s", argv[i], 2))
      mcsetseed(&argv[i][2]);
    else if(!strcmp("--seed", argv[i]) && (i + 1) < argc)
      mcsetseed(argv[++i]);
    else if(!strncmp("--seed=", argv[i], 7))
      mcsetseed(&argv[i][7]);
    else if(!strcmp("-n", argv[i]) && (i + 1) < argc)
      mcsetn_arg(argv[++i]);
    else if(!strncmp("-n", argv[i], 2))
      mcsetn_arg(&argv[i][2]);
    else if(!strcmp("--ncount", argv[i]) && (i + 1) < argc)
      mcsetn_arg(argv[++i]);
    else if(!strncmp("--ncount=", argv[i], 9))
      mcsetn_arg(&argv[i][9]);
    else if(!strcmp("-d", argv[i]) && (i + 1) < argc)
      usedir=argv[++i];  /* will create directory after parsing all arguments (end of this function) */
    else if(!strncmp("-d", argv[i], 2))
      usedir=&argv[i][2];
    else if(!strcmp("--dir", argv[i]) && (i + 1) < argc)
      usedir=argv[++i];
    else if(!strncmp("--dir=", argv[i], 6))
      usedir=&argv[i][6];
    else if(!strcmp("-h", argv[i]))
      mcshowhelp(argv[0]);
    else if(!strcmp("--help", argv[i]))
      mcshowhelp(argv[0]);
    else if(!strcmp("-i", argv[i])) {
      mcformat=FLAVOR_UPPER;
      mcinfo();
    }
    else if(!strcmp("--info", argv[i]))
      mcinfo();
    else if(!strcmp("-t", argv[i]))
      mcenabletrace();
    else if(!strcmp("--trace", argv[i]))
      mcenabletrace();
    else if(!strcmp("--gravitation", argv[i]))
      mcgravitation = 1;
    else if(!strcmp("-g", argv[i]))
      mcgravitation = 1;
    else if(!strncmp("--format=", argv[i], 9)) {
      mcformat=&argv[i][9];
    }
    else if(!strcmp("--format", argv[i]) && (i + 1) < argc) {
      mcformat=argv[++i];
    }
    else if(!strcmp("--no-output-files", argv[i]))
      mcdisable_output_files = 1;
    else if(argv[i][0] != '-' && (p = strchr(argv[i], '=')) != NULL)
    {
      *p++ = '\0';

      for(j = 0; j < mcnumipar; j++)
        if(!strcmp(mcinputtable[j].name, argv[i]))
        {
          int status;
          status = (*mcinputtypes[mcinputtable[j].type].getparm)(p,
                        mcinputtable[j].par);
          if(!status || !strlen(p))
          {
            (*mcinputtypes[mcinputtable[j].type].error)
              (mcinputtable[j].name, p);
            exit(1);
          }
          paramsetarray[j] = 1;
          paramset = 1;
          break;
        }
      if(j == mcnumipar)
      {                                /* Unrecognized parameter name */
        fprintf(stderr, "Error: unrecognized parameter %s (mcparseoptions)\n", argv[i]);
        exit(1);
      }
    }
    else if(argv[i][0] == '-') {
      fprintf(stderr, "Error: unrecognized option argument %s (mcparseoptions). Ignored.\n", argv[i++]);
    }
    else {
      fprintf(stderr, "Error: unrecognized argument %s (mcparseoptions). Aborting.\n", argv[i]);
      mcusage(argv[0]);
    }
  }
  if(!paramset)
    mcreadparams();                /* Prompt for parameters if not specified. */
  else
  {
    for(j = 0; j < mcnumipar; j++)
      if(!paramsetarray[j])
      {
        fprintf(stderr, "Error: Instrument parameter %s left unset (mcparseoptions)\n",
                mcinputtable[j].name);
        exit(1);
      }
  }
  free(paramsetarray);
#ifdef USE_MPI
  if (mcdotrace) mpi_node_count=1; /* disable threading when in trace mode */
#endif
  if (usedir && strlen(usedir) && !mcdisable_output_files) mcuse_dir(usedir);
} /* mcparseoptions */

#ifndef NOSIGNALS
mcstatic char  mcsig_message[256];


/*******************************************************************************
* sighandler: signal handler that makes simulation stop, and save results
*******************************************************************************/
void sighandler(int sig)
{
  /* MOD: E. Farhi, Sep 20th 2001: give more info */
  time_t t1, t0;
#define SIG_SAVE 0
#define SIG_TERM 1
#define SIG_STAT 2
#define SIG_ABRT 3

  printf("\n# " MCCODE_STRING ": [pid %i] Signal %i detected", getpid(), sig);
#ifdef USE_MPI
  printf(" [proc %i]", mpi_node_rank);
#endif
#if defined(SIGUSR1) && defined(SIGUSR2) && defined(SIGKILL)
  if (!strcmp(mcsig_message, "sighandler") && (sig != SIGUSR1) && (sig != SIGUSR2))
  {
    printf("\n# Fatal : unrecoverable loop ! Suicide (naughty boy).\n");
    kill(0, SIGKILL); /* kill myself if error occurs within sighandler: loops */
  }
#endif
  switch (sig) {
#ifdef SIGINT
    case SIGINT : printf(" SIGINT (interrupt from terminal, Ctrl-C)"); sig = SIG_TERM; break;
#endif
#ifdef SIGILL
    case SIGILL  : printf(" SIGILL (Illegal instruction)"); sig = SIG_ABRT; break;
#endif
#ifdef SIGFPE
    case SIGFPE  : printf(" SIGFPE (Math Error)"); sig = SIG_ABRT; break;
#endif
#ifdef SIGSEGV
    case SIGSEGV : printf(" SIGSEGV (Mem Error)"); sig = SIG_ABRT; break;
#endif
#ifdef SIGTERM
    case SIGTERM : printf(" SIGTERM (Termination)"); sig = SIG_TERM; break;
#endif
#ifdef SIGABRT
    case SIGABRT : printf(" SIGABRT (Abort)"); sig = SIG_ABRT; break;
#endif
#ifdef SIGQUIT
    case SIGQUIT : printf(" SIGQUIT (Quit from terminal)"); sig = SIG_TERM; break;
#endif
#ifdef SIGTRAP
    case SIGTRAP : printf(" SIGTRAP (Trace trap)"); sig = SIG_ABRT; break;
#endif
#ifdef SIGPIPE
    case SIGPIPE : printf(" SIGPIPE (Broken pipe)"); sig = SIG_ABRT; break;
#endif
#ifdef SIGUSR1
    case SIGUSR1 : printf(" SIGUSR1 (Display info)"); sig = SIG_STAT; break;
#endif
#ifdef SIGUSR2
    case SIGUSR2 : printf(" SIGUSR2 (Save simulation)"); sig = SIG_SAVE; break;
#endif
#ifdef SIGHUP
    case SIGHUP  : printf(" SIGHUP (Hangup/update)"); sig = SIG_SAVE; break;
#endif
#ifdef SIGBUS
    case SIGBUS  : printf(" SIGBUS (Bus error)"); sig = SIG_ABRT; break;
#endif
#ifdef SIGURG
    case SIGURG  : printf(" SIGURG (Urgent socket condition)"); sig = SIG_ABRT; break;
#endif
#ifdef SIGBREAK
    case SIGBREAK: printf(" SIGBREAK (Break signal, Ctrl-Break)"); sig = SIG_SAVE; break;
#endif
    default : printf(" (look at signal list for signification)"); sig = SIG_ABRT; break;
  }
  printf("\n");
  printf("# Simulation: %s (%s) \n", mcinstrument_name, mcinstrument_source);
  printf("# Breakpoint: %s ", mcsig_message);
  if (strstr(mcsig_message, "Save") && (sig == SIG_SAVE))
    sig = SIG_STAT;
  SIG_MESSAGE("sighandler");
  if (mcget_ncount() == 0)
    printf("(0 %%)\n" );
  else
  {
    printf("%.2f %% (%10.1f/%10.1f)\n", 100.0*mcget_run_num()/mcget_ncount(), 1.0*mcget_run_num(), 1.0*mcget_ncount());
  }
  t0 = (time_t)mcstartdate;
  t1 = time(NULL);
  printf("# Date:      %s", ctime(&t1));
  printf("# Started:   %s", ctime(&t0));

  if (sig == SIG_STAT)
  {
    printf("# " MCCODE_STRING ": Resuming simulation (continue)\n");
    fflush(stdout);
    return;
  }
  else
  if (sig == SIG_SAVE)
  {
    printf("# " MCCODE_STRING ": Saving data and resume simulation (continue)\n");
    mcsave(NULL);
    fflush(stdout);
    return;
  }
  else
  if (sig == SIG_TERM)
  {
    printf("# " MCCODE_STRING ": Finishing simulation (save results and exit)\n");
    mcfinally();
    exit(0);
  }
  else
  {
    fflush(stdout);
    perror("# Last I/O Error");
    printf("# " MCCODE_STRING ": Simulation stop (abort).\n");
// This portion of the signal handling only works on UNIX
#if defined(__unix__) || defined(__APPLE__)
    signal(sig, SIG_DFL); /* force to use default sighandler now */
    kill(getpid(), sig);  /* and trigger it with the current signal */
#endif
    exit(-1);
  }
#undef SIG_SAVE
#undef SIG_TERM
#undef SIG_STAT
#undef SIG_ABRT

} /* sighandler */
#endif /* !NOSIGNALS */

/*******************************************************************************
* mccode_main: McCode main() function.
*******************************************************************************/
int mccode_main(int argc, char *argv[])
{
/*  double run_num = 0; */
  time_t  t;
#ifdef USE_MPI
  char mpi_node_name[MPI_MAX_PROCESSOR_NAME];
  int  mpi_node_name_len;
#endif /* USE_MPI */

#ifdef MAC
  argc = ccommand(&argv);
#endif

#ifdef USE_MPI
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_node_count); /* get number of nodes */
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_node_rank);
  MPI_Comm_set_name(MPI_COMM_WORLD, mcinstrument_name);
  MPI_Get_processor_name(mpi_node_name, &mpi_node_name_len);
#endif /* USE_MPI */

t = time(NULL);
mcseed = (long)t+(long)getpid();

#ifdef USE_MPI
/* *** print number of nodes *********************************************** */
  if (mpi_node_count > 1) {
    MPI_MASTER(
    printf("Simulation '%s' (%s): running on %i nodes (master is '%s', MPI version %i.%i).\n",
      mcinstrument_name, mcinstrument_source, mpi_node_count, mpi_node_name, MPI_VERSION, MPI_SUBVERSION);
    );
  }
#endif /* USE_MPI */
  
  mcstartdate = (long)t;  /* set start date before parsing options and creating sim file */

/* *** parse options ******************************************************* */
  SIG_MESSAGE("main (Start)");
  mcformat=getenv(FLAVOR_UPPER "_FORMAT") ?
           getenv(FLAVOR_UPPER "_FORMAT") : FLAVOR_UPPER;
  mcinstrument_exe = argv[0]; /* store the executable path */
  /* read simulation parameters and options */
  mcparseoptions(argc, argv); /* sets output dir and format */
  
#ifdef USE_MPI
  if (mpi_node_count > 1) {
    /* share the same seed, then adapt random seed for each node */
    MPI_Bcast(&mcseed, 1, MPI_LONG, 0, MPI_COMM_WORLD); /* root sends its seed to slaves */
    mcseed += mpi_node_rank; /* make sure we use different seeds per node */
  }
#endif
  srandom(mcseed);

/* *** install sig handler, but only once !! after parameters parsing ******* */
#ifndef NOSIGNALS
#ifdef SIGQUIT
  if (signal( SIGQUIT ,sighandler) == SIG_IGN)
    signal( SIGQUIT,SIG_IGN);   /* quit (ASCII FS) */
#endif
#ifdef SIGABRT
  if (signal( SIGABRT ,sighandler) == SIG_IGN)
    signal( SIGABRT,SIG_IGN);   /* used by abort, replace SIGIOT in the future */
#endif
#ifdef SIGTERM
  if (signal( SIGTERM ,sighandler) == SIG_IGN)
    signal( SIGTERM,SIG_IGN);   /* software termination signal from kill */
#endif
#ifdef SIGUSR1
  if (signal( SIGUSR1 ,sighandler) == SIG_IGN)
    signal( SIGUSR1,SIG_IGN);   /* display simulation status */
#endif
#ifdef SIGUSR2
  if (signal( SIGUSR2 ,sighandler) == SIG_IGN)
    signal( SIGUSR2,SIG_IGN);
#endif
#ifdef SIGHUP
  if (signal( SIGHUP ,sighandler) == SIG_IGN)
    signal( SIGHUP,SIG_IGN);
#endif
#ifdef SIGILL
  if (signal( SIGILL ,sighandler) == SIG_IGN)
    signal( SIGILL,SIG_IGN);    /* illegal instruction (not reset when caught) */
#endif
#ifdef SIGFPE
  if (signal( SIGFPE ,sighandler) == SIG_IGN)
    signal( SIGSEGV,SIG_IGN);    /* floating point exception */
#endif
#ifdef SIGBUS
  if (signal( SIGBUS ,sighandler) == SIG_IGN)
    signal( SIGSEGV,SIG_IGN);    /* bus error */
#endif
#ifdef SIGSEGV
  if (signal( SIGSEGV ,sighandler) == SIG_IGN)
    signal( SIGSEGV,SIG_IGN);   /* segmentation violation */
#endif
#endif /* !NOSIGNALS */
  mcsiminfo_init(NULL); /* open SIM */
  SIG_MESSAGE("main (Init)");
  mcinit();
#ifndef NOSIGNALS
#ifdef SIGINT
  if (signal( SIGINT ,sighandler) == SIG_IGN)
    signal( SIGINT,SIG_IGN);    /* interrupt (rubout) only after INIT */
#endif
#endif /* !NOSIGNALS */

/* ================ main particle generation/propagation loop ================ */
#if defined (USE_MPI)
  /* sliced Ncount on each MPI node */
  mcncount = mpi_node_count > 1 ?
    floor(mcncount / mpi_node_count) :
    mcncount; /* number of rays per node */
#endif

/* main particle event loop */
while(mcrun_num < mcncount || mcrun_num < mcget_ncount())
  {
#ifndef NEUTRONICS
    mcgenstate();
#endif
    /* old init: mcsetstate(0, 0, 0, 0, 0, 1, 0, sx=0, sy=1, sz=0, 1); */
    mcraytrace();
    mcrun_num++;
  }

#ifdef USE_MPI
 /* merge run_num from MPI nodes */
  if (mpi_node_count > 1) {
  double mcrun_num_double = (double)mcrun_num;
  mc_MPI_Sum(&mcrun_num_double, 1);
  mcrun_num = (unsigned long long)mcrun_num_double;
  }
#endif

/* save/finally executed by master node/thread */
  mcfinally();

#ifdef USE_MPI
  MPI_Finalize();
#endif /* USE_MPI */

  return 0;
} /* mccode_main */

#ifdef NEUTRONICS
/*Main neutronics function steers the McStas calls, initializes parameters etc */
/* Only called in case NEUTRONICS = TRUE */
void neutronics_main_(float *inx, float *iny, float *inz, float *invx, float *invy, float *invz, float *intime, float *insx, float *insy, float *insz, float *inw, float *outx, float *outy, float *outz, float *outvx, float *outvy, float *outvz, float *outtime, float *outsx, float *outsy, float *outsz, float *outwgt)
{

  extern double mcnx, mcny, mcnz, mcnvx, mcnvy, mcnvz;
  extern double mcnt, mcnsx, mcnsy, mcnsz, mcnp;

  /* External code governs iteration - McStas is iterated once per call to neutronics_main. I.e. below counter must be initiancated for each call to neutronics_main*/
  mcrun_num=0;

  time_t t;
  t = (time_t)mcstartdate;
  mcstartdate = t;  /* set start date before parsing options and creating sim file */
  mcinit();

  /* *** parse options *** */
  SIG_MESSAGE("main (Start)");
  mcformat=getenv(FLAVOR_UPPER "_FORMAT") ?
           getenv(FLAVOR_UPPER "_FORMAT") : FLAVOR_UPPER;

  /* Set neutron state based on input from neutronics code */
  mcsetstate(*inx,*iny,*inz,*invx,*invy,*invz,*intime,*insx,*insy,*insz,*inw);

  /* main neutron event loop - runs only one iteration */

  //mcstas_raytrace(&mcncount); /* prior to McStas 1.12 */

  mcallowbackprop = 1; //avoid absorbtion from negative dt
  int argc=1;
  char *argv[0];
  int dummy = mccode_main(argc, argv);

  *outx =  mcnx;
  *outy =  mcny;
  *outz =  mcnz;
  *outvx =  mcnvx;
  *outvy =  mcnvy;
  *outvz =  mcnvz;
  *outtime =  mcnt;
  *outsx =  mcnsx;
  *outsy =  mcnsy;
  *outsz =  mcnsz;
  *outwgt =  mcnp;

  return;
} /* neutronics_main */

#endif /*NEUTRONICS*/

#endif /* !MCCODE_H */
/* End of file "mccode-r.c". */
/* End of file "mccode-r.c". */

#line 4977 "./ess_extraction.c"

#line 1 "mcstas-r.c"
/*******************************************************************************
*
* McStas, neutron ray-tracing package
*         Copyright (C) 1997-2009, All rights reserved
*         Risoe National Laboratory, Roskilde, Denmark
*         Institut Laue Langevin, Grenoble, France
*
* Runtime: share/mcstas-r.c
*
* %Identification
* Written by: KN
* Date:    Aug 29, 1997
* Release: McStas X.Y
* Version: $Revision$
*
* Runtime system for McStas.
* Embedded within instrument in runtime mode.
*
* Usage: Automatically embbeded in the c code whenever required.
*
* $Id$
*
*******************************************************************************/

#ifndef MCSTAS_R_H
#include "mcstas-r.h"
#endif
#ifdef DANSE
#include "mcstas-globals.h"
#endif

/*******************************************************************************
* The I/O format definitions and functions
*******************************************************************************/

/*the magnet stack*/
#ifdef MC_POL_COMPAT
void (*mcMagnetPrecession) (double, double, double, double, double, double,
    double, double*, double*, double*, double, Coords, Rotation)=NULL;
Coords   mcMagnetPos;
Rotation mcMagnetRot;
double*  mcMagnetData                = NULL;
/* mcMagneticField(x, y, z, t, Bx, By, Bz) */
int (*mcMagneticField) (double, double, double, double,
    double*, double*, double*, void *) = NULL;
#endif

#ifndef MCSTAS_H

/*******************************************************************************
* mcstore_neutron: stores neutron coodinates into global array (per component)
*******************************************************************************/
void
mcstore_neutron(MCNUM *s, int index, double x, double y, double z,
               double vx, double vy, double vz, double t,
               double sx, double sy, double sz, double p)
{
    double *dptr = &s[11*index];
    *dptr++  = x;
    *dptr++  = y ;
    *dptr++  = z ;
    *dptr++  = vx;
    *dptr++  = vy;
    *dptr++  = vz;
    *dptr++  = t ;
    *dptr++  = sx;
    *dptr++  = sy;
    *dptr++  = sz;
    *dptr    = p ;
} /* mcstore_neutron */

/*******************************************************************************
* mcrestore_neutron: restores neutron coodinates from global array
*******************************************************************************/
void
mcrestore_neutron(MCNUM *s, int index, double *x, double *y, double *z,
               double *vx, double *vy, double *vz, double *t,
               double *sx, double *sy, double *sz, double *p)
{
    double *dptr = &s[11*index];
    *x  =  *dptr++;
    *y  =  *dptr++;
    *z  =  *dptr++;
    *vx =  *dptr++;
    *vy =  *dptr++;
    *vz =  *dptr++;
    *t  =  *dptr++;
    *sx =  *dptr++;
    *sy =  *dptr++;
    *sz =  *dptr++;
    *p  =  *dptr;
} /* mcrestore_neutron */

/*******************************************************************************
* mcsetstate: transfer parameters into global McStas variables 
*******************************************************************************/
void
mcsetstate(double x, double y, double z, double vx, double vy, double vz,
           double t, double sx, double sy, double sz, double p)
{
  extern double mcnx, mcny, mcnz, mcnvx, mcnvy, mcnvz;
  extern double mcnt, mcnsx, mcnsy, mcnsz, mcnp;

  mcnx = x;
  mcny = y;
  mcnz = z;
  mcnvx = vx;
  mcnvy = vy;
  mcnvz = vz;
  mcnt = t;
  mcnsx = sx;
  mcnsy = sy;
  mcnsz = sz;
  mcnp = p;
} /* mcsetstate */

/*******************************************************************************
* mcgenstate: set default neutron parameters 
*******************************************************************************/
void
mcgenstate(void)
{
  mcsetstate(0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1);
  /* old initialisation: mcsetstate(0, 0, 0, 0, 0, 1, 0, sx=0, sy=1, sz=0, 1); */
}

/* intersection routines ==================================================== */

/*******************************************************************************
* inside_rectangle: Check if (x,y) is inside rectangle (xwidth, yheight) 
* return 0 if outside and 1 if inside 
*******************************************************************************/
int inside_rectangle(double x, double y, double xwidth, double yheight)
{
  if (x>-xwidth/2 && x<xwidth/2 && y>-yheight/2 && y<yheight/2)
    return 1;
  else
    return 0;
}

/*******************************************************************************
 * box_intersect: compute time intersection with a box
 * returns 0 when no intersection is found
 *      or 1 in case of intersection with resulting times dt_in and dt_out
 * This function written by Stine Nyborg, 1999. 
 *******************************************************************************/
int box_intersect(double *dt_in, double *dt_out,
                  double x, double y, double z,
                  double vx, double vy, double vz,
                  double dx, double dy, double dz)
{
  double x_in, y_in, z_in, tt, t[6], a, b;
  int i, count, s;

      /* Calculate intersection time for each of the six box surface planes
       *  If the box surface plane is not hit, the result is zero.*/

  if(vx != 0)
   {
    tt = -(dx/2 + x)/vx;
    y_in = y + tt*vy;
    z_in = z + tt*vz;
    if( y_in > -dy/2 && y_in < dy/2 && z_in > -dz/2 && z_in < dz/2)
      t[0] = tt;
    else
      t[0] = 0;

    tt = (dx/2 - x)/vx;
    y_in = y + tt*vy;
    z_in = z + tt*vz;
    if( y_in > -dy/2 && y_in < dy/2 && z_in > -dz/2 && z_in < dz/2)
      t[1] = tt;
    else
      t[1] = 0;
   }
  else
    t[0] = t[1] = 0;

  if(vy != 0)
   {
    tt = -(dy/2 + y)/vy;
    x_in = x + tt*vx;
    z_in = z + tt*vz;
    if( x_in > -dx/2 && x_in < dx/2 && z_in > -dz/2 && z_in < dz/2)
      t[2] = tt;
    else
      t[2] = 0;

    tt = (dy/2 - y)/vy;
    x_in = x + tt*vx;
    z_in = z + tt*vz;
    if( x_in > -dx/2 && x_in < dx/2 && z_in > -dz/2 && z_in < dz/2)
      t[3] = tt;
    else
      t[3] = 0;
   }
  else
    t[2] = t[3] = 0;

  if(vz != 0)
   {
    tt = -(dz/2 + z)/vz;
    x_in = x + tt*vx;
    y_in = y + tt*vy;
    if( x_in > -dx/2 && x_in < dx/2 && y_in > -dy/2 && y_in < dy/2)
      t[4] = tt;
    else
      t[4] = 0;

    tt = (dz/2 - z)/vz;
    x_in = x + tt*vx;
    y_in = y + tt*vy;
    if( x_in > -dx/2 && x_in < dx/2 && y_in > -dy/2 && y_in < dy/2)
      t[5] = tt;
    else
      t[5] = 0;
   }
  else
    t[4] = t[5] = 0;

  /* The intersection is evaluated and *dt_in and *dt_out are assigned */

  a = b = s = 0;
  count = 0;

  for( i = 0; i < 6; i = i + 1 )
    if( t[i] == 0 )
      s = s+1;
    else if( count == 0 )
    {
      a = t[i];
      count = 1;
    }
    else
    {
      b = t[i];
      count = 2;
    }

  if ( a == 0 && b == 0 )
    return 0;
  else if( a < b )
  {
    *dt_in = a;
    *dt_out = b;
    return 1;
  }
  else
  {
    *dt_in = b;
    *dt_out = a;
    return 1;
  }

} /* box_intersect */

/*******************************************************************************
 * cylinder_intersect: compute intersection with a cylinder
 * returns 0 when no intersection is found
 *      or 2/4/8/16 bits depending on intersection,
 *     and resulting times t0 and t1
 * Written by: EM,NB,ABA 4.2.98 
  *******************************************************************************/
int
cylinder_intersect(double *t0, double *t1, double x, double y, double z,
                   double vx, double vy, double vz, double r, double h)
{
  double D, t_in, t_out, y_in, y_out;
  int ret=1;

  D = (2*vx*x + 2*vz*z)*(2*vx*x + 2*vz*z)
    - 4*(vx*vx + vz*vz)*(x*x + z*z - r*r);

  if (D>=0)
  {
    if (vz*vz + vx*vx) {
      t_in  = (-(2*vz*z + 2*vx*x) - sqrt(D))/(2*(vz*vz + vx*vx));
      t_out = (-(2*vz*z + 2*vx*x) + sqrt(D))/(2*(vz*vz + vx*vx));
    } else if (vy) { /* trajectory parallel to cylinder axis */
      t_in = (-h/2-y)/vy;
      t_out = (h/2-y)/vy;
      if (t_in>t_out){
        double tmp=t_in;
        t_in=t_out;t_out=tmp;
      }
    } else return 0;
    y_in = vy*t_in + y;
    y_out =vy*t_out + y;

    if ( (y_in > h/2 && y_out > h/2) || (y_in < -h/2 && y_out < -h/2) )
      return 0;
    else
    {
      if (y_in > h/2)
        { t_in = ((h/2)-y)/vy; ret += 2; }
      else if (y_in < -h/2)
        { t_in = ((-h/2)-y)/vy; ret += 4; }
      if (y_out > h/2)
        { t_out = ((h/2)-y)/vy; ret += 8; }
      else if (y_out < -h/2)
        { t_out = ((-h/2)-y)/vy; ret += 16; }
    }
    *t0 = t_in;
    *t1 = t_out;
    return ret;
  }
  else
  {
    *t0 = *t1 = 0;
    return 0;
  }
} /* cylinder_intersect */


/*******************************************************************************
 * sphere_intersect: Calculate intersection between a line and a sphere.
 * returns 0 when no intersection is found
 *      or 1 in case of intersection with resulting times t0 and t1 
 *******************************************************************************/
int
sphere_intersect(double *t0, double *t1, double x, double y, double z,
                 double vx, double vy, double vz, double r)
{
  double A, B, C, D, v;

  v = sqrt(vx*vx + vy*vy + vz*vz);
  A = v*v;
  B = 2*(x*vx + y*vy + z*vz);
  C = x*x + y*y + z*z - r*r;
  D = B*B - 4*A*C;
  if(D < 0)
    return 0;
  D = sqrt(D);
  *t0 = (-B - D) / (2*A);
  *t1 = (-B + D) / (2*A);
  return 1;
} /* sphere_intersect */

/*******************************************************************************
 * plane_intersect: Calculate intersection between a plane and a line.
 * returns 0 when no intersection is found (i.e. line is parallel to the plane)
 * returns 1 or -1 when intersection time is positive and negative respectively
 *******************************************************************************/
int
plane_intersect(double *t, double x, double y, double z,
                 double vx, double vy, double vz, double nx, double ny, double nz, double wx, double wy, double wz)
{
  double s;
  if (fabs(s=scalar_prod(nx,ny,nz,vx,vy,vz))<FLT_EPSILON) return 0;
  *t = - scalar_prod(nx,ny,nz,x-wx,y-wy,z-wz)/s;
  if (*t<0) return -1;
  else return 1;
} /* plane_intersect */

#endif /* !MCSTAS_H */
/* End of file "mcstas-r.c". */

#line 5337 "./ess_extraction.c"
#ifdef MC_TRACE_ENABLED
int mctraceenabled = 1;
#else
int mctraceenabled = 0;
#endif
#define MCSTAS "/usr/share/mcstas/2.6.1/tools/Python/mcrun/../mccodelib/../../../"
int mcdefaultmain = 1;
char mcinstrument_name[] = "template_simple";
char mcinstrument_source[] = "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr";
char *mcinstrument_exe=NULL; /* will be set to argv[0] in main */
int main(int argc, char *argv[]){return mccode_main(argc, argv);}
void mcinit(void);
void mcraytrace(void);
void mcsave(FILE *);
void mcfinally(void);
void mcdisplay(void);

/* Shared user declarations for all components 'PSD_monitor'. */
#line 57 "/usr/share/mcstas/2.6.1/tools/Python/mcrun/../mccodelib/../../../monitors/PSD_monitor.comp"

#ifndef ARRAYS_H
#define ARRAYS_H
typedef double* DArray1d;
DArray1d create_darr1d(int n);
void destroy_darr1d(DArray1d a);

typedef double** DArray2d;
DArray2d create_darr2d(int nx, int ny);
void destroy_darr2d(DArray2d a);

typedef double*** DArray3d;
DArray3d create_darr3d(int nx, int ny, int nz);
void destroy_darr3d(DArray3d a);
#endif
#ifndef ARRAYS_C
#define ARRAYS_C
#include <stdlib.h>

DArray1d create_darr1d(int n){
  DArray1d arr2d;
  arr2d = calloc(n, sizeof(double));
  return arr2d;
}
void destroy_darr1d(DArray1d a){
  free(a);
}

DArray2d create_darr2d(int nx, int ny){
  DArray2d arr2d;
  arr2d = calloc(nx, sizeof(double *));

  double *p1;
  p1 = calloc(nx*ny, sizeof(double));

  int i;
  for (i=0; i<nx; i++){
    arr2d[i] = &(p1[i*ny]);
  }
  return arr2d;
}
void destroy_darr2d(DArray2d a){
  free(a[0]);
  free(a);
}

DArray3d create_darr3d(int nx, int ny, int nz){
  DArray3d arr3d;
  int i, j;

  // 1d
  arr3d = calloc(nx, sizeof(double **));

  // d2
  double **p1;
  p1 = calloc(nx*ny, sizeof(double *));

  for (i=0; i<nx; i++){
    arr3d[i] = &(p1[i*ny]);
  }

  // 3d
  double *p2;
  p2 = calloc(nx*ny*nz, sizeof(double));
  for (i=0; i<nx; i++){
    for (j=0; j<ny; j++){
      arr3d[i][j] = &(p2[(i*ny+j)*nz]);
    }
  }
  return arr3d;
}
void destroy_darr3d(DArray3d a){
  free(a[0][0]);
  free(a[0]);
  free(a);
}
#endif

#line 5435 "./ess_extraction.c"

/* Shared user declarations for all components 'FlatEllipse_finite_mirror'. */
#line 53 "FlatEllipse_finite_mirror.comp"
/**
\mainpage
Simple Meta-Conic Neutron Raytracer is a framework for raytracing geometries of the form: @f$ r^2=k_1 + k_2 z + k_3 z^2 @f$.

<h3>General Notes</h3>
To use the software you must make a Scene element using the function makeScene(). You must then add items to this scene element using the various add function (addDisk(), addParaboloid(), etc...). Next you must call the function traceSingleNeutron() for every neutron you would like to trace through the geometry. The maximum number of each geometry you can place in a scene is defined by the MAX_CONICSURF, MAX_DISK and MAX_DETECTOR definitions in the conic.h file.

<h3>TODO</h3>

@todo
       Name variable for each component <br/>
       Normalize Detector Events by weight of neutron <br/>

<h3>Known Bugs</h3>
@bug  HPPH works incorrectly for M != 1 <br/>
      Neutrons t=1e-11 away from a surface will pass through <br/>

<h3>Note on Pointers</h3>
This framework uses pointers extensivly. It is important to be familiar with how to use pointers.
<br/>Here is a quick overview.

@code
//Making an Element
ConicSurf c = makeParaboloid(...);
int k = 10;

//Making a Pointer from an Element
ConicSurf * c_pointer = &c;
int * k_pointer = &k;

//Making an Element from a Pointer
ConicSurf c2 = *c_pointer;
int ten = *k_pointer;

//Getting Item in Element
double k1 = c.k1;

//Getting Item in Element from a Pointer
double k1 = c_pointer->k1;

//Functions that have pointers as parameters can modify the element being pointed to.
//Functions that have elements as parameters can not modify the value of the element.
@endcode

<h3>Stand Alone Example Code</h3>
This framework can be used entirely by itself as the following example demonstrates.
@code
/////////////////////////////////////////////////////////////////////////////////
// Giacomo Resta <gresta@mit.edu>
//
// Basic standalone exmaple of a raytracer. It is advised to modify
// the getRandom() function to use a better random number generator (i.e. glib).
// The systems random generator was used to preserve clarity and conciseness.
/////////////////////////////////////////////////////////////////////////////////

#include <math.h>
#include <stdio.h>

#include "conic.h"
#include "w1_general.h"

#define NUM_NEUTRON 1000000
*/
//Function to get random number
double getRandom() {
    return (double)lrand48()/RAND_MAX;
}
/*
//Function to get new particle from source
Particle generateParticleFromSource(double radius, double div, double Lambda0,
        double num_neutrons) {

    double chi = 2*M_PI*(getRandom());
    double r = sqrt(getRandom()) * radius;

    double v = 3956.036 / Lambda0;

    double theta = sqrt(getRandom())*div;
    double phi = getRandom()*2*M_PI;
    double tan_r = tan(theta);

    double vz = v/sqrt(1+tan_r*tan_r);
    double vr = tan_r * vz;

    return makeParticle(r*cos(chi),r*sin(chi),0,cos(phi)*vr,sin(phi)*vr,
                vz,0,0.0,0.0,0.0,1.0/num_neutrons);
}

//Function to add items to scene
void addItems(Scene* s, double instr_len,double r,double f,double M,
        double max_mir_len,double m, double mirr_thick) {

    //Change code here to make different Scenes
    PP p = addPPShell(0.0, instr_len, r, f, M, max_mir_len, m, m, mirr_thick, s);
    addDisk(p.p0->zs, 0.0, rConic(p.p0->ze, *p.p0)*p.p0->zs/p.p0->ze, s);
    addDisk(p.p1->zs, rConic(p.p1->zs, *p.p1), 10000,s);
    addDetector(10.0, -0.01, 0.01, -0.01, 0.01, 600, 600, NUM_NEUTRON, "test.txt", s);
}

//Main Function
int main() {
    //seed random generator (starts the generator)
    srand48((long)time(0));

    //Make a new Scene
    Scene s = makeScene();

    //Add Items and initialize Scene
    addItems(&s,10.0,0.068,4.2,1,0.7,3,0.001);
    initSimulation(&s);

    //Raytrace all particles through the Scene
    double i;
    for (i = 0; i < NUM_NEUTRON; i++) {
        Particle p = generateParticleFromSource(0.005, 0.02422, 4, NUM_NEUTRON);
        traceSingleNeutron(&p,s);
    }

    //Finish Simulation of the Scene
    finishSimulation(&s);

    return 0;
}
@endcode

*/

/**
    @file conic.h
    \brief General Framework for generating and raytracing geometries of the
    form @f$ r = k_1+k_2 z+k_3 z^2 @f$

    @author Giacomo Resta <gresta@mit.edu>
    @version 0.2

    @section LICENSE
    Permission is hereby granted, free of charge, to any person obtaining a
    copy of this software and associated documentation files (the "Software"),
    to deal in the Software without restriction, including without limitation
    the rights to use, copy, modify, merge, publish, distribute, sublicense,
    and/or sell copies of the Software, and to permit persons to whom the Software
    is furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
    INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
    PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE
    FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
    OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
    DEALINGS IN THE SOFTWARE.

    @section DESCRIPTION
     General Framework for generating and raytracing geometries of the form
     @f$ r = k_1+k_2 z+k_3 z^2 @f$
*/

/////////////////////////////////////
// Simulation
/////////////////////////////////////

#ifndef MIT_CONICS
#define MIT_CONICS

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/** @defgroup simgroup Simulator Internals
    Contains items general to the simulation
    @{
*/
//! Max number of ConicSurf allowed in a Scene
#define MAX_FLATSURF 100

//! Max number of ConicSurf allowed in a Scene
#define MAX_CONICSURF 100

//! Max number of Disks allowed in a Scene
#define MAX_DISK 100

//! Max number of Detectors allowed in a Scene
#define MAX_DETECTOR 10

//! If "1" simulator will record z location where neutron with greatest grazing angle reflected for each ConicSurf
/*! The information is stored in the max_ga and max_ga_z0 members of each ConicSurf, which are only present if
    this flag is 1. See source code for clarification.

@note You must use default traceNeutronConic function or implement saving routine yourself */
#define REC_MAX_GA 0

#define V2Q_conic 1.58825361e-3
#define Q2V_conic 629.622368
//! Stucture to represent a point
typedef struct {
    double x; //!< x-axis position of point
    double y; //!< y-axis position of point
    double z; //!< z-axis position of point
} Point;

//! Structure to represent a vector
typedef struct {
    double x; //!< x-axis length of vector
    double y; //!< y-axis length of vector
    double z; //!< z-axis length of vector
} Vec;

//! Structure to represent a particle
typedef struct {
    double _x; //!< x axis position of particle
    double _y; //!< y axis position of particle
    double _z; //!< z axis position of particle
    double _vx; //!< x axis components of velocity
    double _vy; //!< y axis components of velocity
    double _vz; //!< z axis components of velocity
    double _sx; //!< x spin vector components
    double _sy; //!< y spin vector components
    double _sz; //!< z spin vector components
    double w; //!< Weight of particle
    int silicon; //!< +1 if Particle is in silicon -1 if Particle is in air
    int absorb; //!< Absorb flag (0 is not absorbed)
    double _t; //!< Time of flight of particle
} Particle;

/*! \brief Function to make a point

@param x x-axis position
@param y y-axis position
@param z z-axis position
*/
Point makePoint(double x, double y, double z) {
    Point p;
    p.x = x;
    p.y = y;
    p.z = z;

    return p;
}

/*! \brief Function to make a vector

@param x x-axis length
@param y y-axis length
@param z z-axis length
*/
Vec makeVec(double x, double y, double z) {
    Vec p;
    p.x = x;
    p.y = y;
    p.z = z;

    return p;
}

//! Function to compute length of a vector
double getMagVec(Vec v) {
    return sqrt(v.x*v.x+v.y*v.y+v.z*v.z);
}

//! Function to compute dot product of two vectors
double dotVec(Vec v1, Vec v2) {
    return v1.x*v2.x+v1.y*v2.y+v1.z*v2.z;
}

/*! \brief Function to make a particle

@param x x-axis position
@param y y-axis position
@param z z-axis position
@param vx x-axis velocity
@param vy y-axis velocity
@param vz z-axis velocity
@param t time of flight of neutron
@param sx x-axis component of spin vector
@param sy y-axis component of spin vector
@param sz z-axis component of spin vector
@param w weight of particle
*/
Particle makeParticle(double x, double y, double z,
    double vx, double vy, double vz, double t,
    double sx, double sy, double sz, int silicon, double w) {

    Particle pa;

    pa._x = x;
    pa._y = y;
    pa._z = z;

    pa._vx = vx;
    pa._vy = vy;
    pa._vz = vz;

    pa._sx = sx;
    pa._sy = sy;
    pa._sz = sz;
    pa.silicon = silicon;
    pa.w = w;

    pa.absorb = 0;
    pa._t = t;

    return pa;
}

//! Function to get position of particle
Point getParticlePos(Particle p) {
    return makePoint(p._x,p._y,p._z);
}

//! Function to get velocity vector of particle
Vec getParticleVel(Particle p) {
    return makeVec(p._vx, p._vy, p._vz);
}

/*! \brief Function to move particle a specific time step.

Will not move particle if t < 0. Does not simulate
gravity.

@param t time step to move particle
@param p pointer to particle
*/
void moveParticleT(double t, Particle* p) {
    if (t < 0)
        return;
    if (p->silicon>0){
        p->_x = p->_x+p->_vx*t;
        p->_y = p->_y+p->_vy*t;
        p->_z = p->_z+p->_vz*t;
        p->_t = p->_t+t;
        p->w *= exp(-t*98900/52.338);//ugly hard coding in the penetration depth of silicon at 989 m/s
        //absorbParticle(p);
    }
    else{
        p->_x = p->_x+p->_vx*t;
        p->_y = p->_y+p->_vy*t;
        p->_z = p->_z+p->_vz*t;
        p->_t = p->_t+t;
    }
}

/*! \brief Function to move particle to position z.

Will not move particle if moving particle to z
position would require negative time.
Does not simulate gravity.

@param z z-position to move particle
@param p pointer to particle
*/
void moveParticleZ(double z, Particle* p) {
    double t = (z-p->_z)/p->_vz;
    moveParticleT(t, p);
}

/*! \brief Function to compute new position of particle
without modifying the position of the actual particle.

Will not move particle if t < 0. Does not simulate gravity.

@param t timestep to move particle
@param p particle
*/
Particle copyMoveParticleT(double t, Particle p) {
    Particle p2 = p;
    moveParticleT(t,&p2);
    return p2;
}

/*! \brief Function to move particle to position z
without modifying the position of the actual particle.

Will not move particle if moving particle to z
position would require negative time.
Does not simulate gravity.

@param z z-position to move particle
@param p pointer to particle
*/
Particle copyMoveParticleZ(double z, Particle p) {
    Particle p2 = p;
    moveParticleZ(z,&p2);
    return p2;
}

/*! \brief Mathematical Aid for Snell's Law for reflection.

Does not take into account grazing angle constrictions.
Only computes mathematical reflection.

@param n Normal vector
@param p Pointer to particle
*/
void reflectParticle(Vec n, Particle* p) {
    double vn = dotVec(getParticleVel(*p),n);

    p->_vx = p->_vx-2*vn*n.x;
    p->_vy = p->_vy-2*vn*n.y;
    p->_vz = p->_vz-2*vn*n.z;
}

/*! \brief Function to mark particle as absorbed

@param p Pointer to particle to be absorbed */
void absorbParticle(Particle* p)  {
    p->_vx = 0;
    p->_vy = 0;
    p->_vz = 0;
    p->w = 0;
    p->absorb = 1;
}

/*! \brief Function to set weight of particle.

Will set the weight of the particle to w.  */
void setWeightParticle(double w, Particle* pa) {
    pa->w = w;
}

/*! \brief Function to solve quadratic equations for smallest positive result.

If no positive result returns -1. Parameters are coefficents such that
@f$ 0=A z^2 + B z + C @f$

@return Will return smallest positive value or -1 if no smallest positive value
*/
double solveQuad(double A, double B, double C) {
    if (fabs(A) < 1e-11 && B != 0)           //FIXME: 1e-11 cutoff may cause problems
        return -C/B;
    else {
        double det = B*B - 4*A*C;
        if (det < 0)
            return -1;
        else {
            double sdet = sqrt(det);
            double s1 = (-B+sdet)/(2*A);
            double s2 = (-B-sdet)/(2*A);

            if (fabs(s1) < 1e-11) s1=0.0;     //FIXME: 1e-11 cutoff may cause problems
            if (fabs(s2) < 1e-11) s2=0.0;     //FIXME: 1e-11 cutoff may cause problems

            if (s1 > 0.0) {
                if (s2 > 0.0) {
                    if (s1 > s2)
                        return s2;
                    return s1;
                }
                else
                    return s1;
           }
           if (s2 > 0.0)
               return s2;
        }
    }
    return -1;
}

//! Returns sign of x, either -1 or 1
int sign(double x) {
    if (x < 0)
        return -1;
    return 1;
}

/** @} */ //end of simgroup

/*! @ingroup detectorgroup
\brief Structure for representing inline detectors

@warning Do not directly modify this structure*/
typedef struct {
    double z0; //!< z-axis position of detector
    double xmin; //!< Smallest x value to detect
    double xmax; //!< Largest x value to detect
    double xstep; //!< x size of subsampling pixel
    double ymin; //!< Smallest y value to detect
    double ymax; //!< Largest y value to detect
    double ystep; //!< y size of subsampling pixel
    int nrows;    //!< Number of pixels along y axis
    int ncols;    //!< Number of pixels along x axis
    double num_particles; //!< Number of particles being emitted from source
    double *num_count; //!< Pointer to the number of particles that hit detector
    double **data; //!< Pointer to 2d data, pixel_x_y = data[x][y]
    char* filename; //!< Name of output file of detector (should end in .txt)
} Detector;

/*! @ingroup diskgroup
\brief Structure for representing Disk geometry

 Creates a doughnut with inner radius r0 and outer radus r1 at position z0.
Neutrons between r0 and r1 are absorbed. */
typedef struct {
    double r0; //!< Inner radius of doughnut
    double r1; //!< Outer radius of doughnut
    double z0; //!< z-axis position of Disk
} Disk;

/*! @ingroup conicgroup */
enum ConicType {
    PARA,
    HYPER,
    ELLIP
};


/*! @ingroup conicgroup
\brief Structure to contain z-axis symetric conic sections

Contains any geometry that can be expressed as
@f$ r^2=k_1 + k_2 z + k_3 z^2 @f$

@warning Do not directly modify values in this structure directly */
typedef struct {
    double k1; //!< @f$ k_1 @f$ in equation below
    double k2; //!< @f$ k_2 @f$ in equation below
    double k3; //!< @f$ k_3 @f$ in equation below
    double zs; //!< z-axis position of start of mirror
    double ze; //!< z-axis position of end of mirror
    double m;  //!< m value for mirror (1.0 for Nickel)
    int doubleReflections; //!< 0 if reflections from the back of the surface cannot happen 1 otherwise
    //Only for reference
    double f1; //!< z-axis position of first focus
    double f2; //!< z-axis position of second focus, for paraboloid this is unassigned
    double a;  //!< Value of a, specific to geometry type
    double c;  //!< Value of c, for paraboloid this is unassigned
    enum ConicType type; //!< Type of mirror geometry

    #if REC_MAX_GA
    double max_ga; //!< Max Grazing Angle of Reflected Neutron (Exists only if REC_MAX_GA)
    double max_ga_z0; //!< Collision point of Max Grazing Neutron (Exists only if REC_MAX_GA)
    #endif

} ConicSurf;


/*! @ingroup flatgroup
\brief Structure to contain z-axis symetric flat sections

Contains flat geometries which can be expressed as
@f$ x = k_1 + k_2 z + k_3 z^2 @f$ or
@f$ y = k_1 + k_2 z + k_3 z^2 @f$

@warning Do not directly modify values in this structure directly */
typedef struct {
    double k1; //!< @f$ k_1 @f$ in equation below
    double k2; //!< @f$ k_2 @f$ in equation below
    double k3; //!< @f$ k_3 @f$ in equation below
    double zs; //!< z-axis position of start of mirror
    double ze; //!< z-axis position of end of mirror
    double ll; //!< left/lower limit of mirror along translational symmetry
    double rl; //!< right/upper limit of mirror along translational symmetry
    double m;  //!< m value for mirror (1.0 for Nickel)

    //Only for reference
    double f1; //!< z-axis position of first focus
    double f2; //!< z-axis position of second focus, for paraboloid this is unassigned
    double a;  //!< Value of a, specific to geometry type
    double c;  //!< Value of c, for paraboloid this is unassigned
    //enum FlatType type; //!< Type of mirror geometry
    int doubleReflections; // will determine whether the geometry allows double reflections
    #if REC_MAX_GA
    double max_ga; //!< Max Grazing Angle of Reflected Neutron (Exists only if REC_MAX_GA)
    double max_ga_z0; //!< Collision point of Max Grazing Neutron (Exists only if REC_MAX_GA)
    #endif

} FlatSurf;
/*! @ingroup simgroup
\brief Structure to hold all scene geometry

The number of possible ConicSurf, Disk and Detector in the Scene are
determined by MAX_CONICSURF, MAX_DISK and MAX_DETECTOR.
*/
typedef struct {
    FlatSurf f[MAX_FLATSURF]; //!< Array of all ConicSurf in Scene
    int num_f;                  //!< Number of ConicSurf in Scene

    ConicSurf c[MAX_CONICSURF]; //!< Array of all ConicSurf in Scene
    int num_c;                  //!< Number of ConicSurf in Scene

    Disk di[MAX_DISK];          //!< Array of all Disk in Scene
    int num_di;                 //!< Number of Disk in Scene

    Detector d[MAX_DETECTOR];  //!< Array of all Detector in Scene
    int num_d;                 //!< Number of Detector in Scene

    //! Function called to handle Neutron-Flat Interaction
    void (*traceNeutronFlat)(Particle*,FlatSurf);
    //! Function called to handle Neutron-Conic Interaction
    void (*traceNeutronConic)(Particle*,ConicSurf);
    //! Function called to handle Neutron-Disk Interaction
    void (*traceNeutronDisk)(Particle*,Disk);
    //! Function called to handle Neutron-Detector Interaction
    void (*traceNeutronDetector)(Particle*,Detector);

} Scene;

/////////////////////////////////////
// Inline Detector
/////////////////////////////////////

/** @defgroup detectorgroup Detector
    Contains code related to the inline detectors
    @{
*/

/*! \brief Function to make Detector

@param z0 z-axis position of detector
@param xmin Smallest x value to detect
@param xmax Largest x value to detect
@param ymin Smallest y value to detect
@param ymax Largest y value to detect
@param xres Number of pixels along x axis
@param yres Number of pixels along y axis
@param num_particles Total number of particles being emitted
@param filename Name of output file of detector (should end in .txt)
*/
Detector makeDetector(double z0,double xmin, double xmax, double ymin, double ymax, int xres,
    int yres, double num_particles, char* filename) {

    Detector d;
    d.z0 = z0;
    d.xmin = xmin;
    d.xmax = xmax;
    d.xstep = (xmax-xmin)/xres;

    d.ymin = ymin;
    d.ymax = ymax;
    d.ystep = (ymax-ymin)/yres;

    d.ncols = xres;
    d.nrows = yres;
    d.filename = filename;
    d.num_particles = num_particles;

    d.num_count = (double*)malloc(sizeof(double));
    if (d.num_count == NULL) {
        fprintf(stderr, "MEMORY ALLOCATION PROBLEM\n");
        exit(-1);
    }
    (*d.num_count) = 0;

    d.data = (double**)malloc(d.ncols*sizeof(double *));
    if (d.data == NULL) {
        fprintf(stderr, "MEMORY ALLOCATION PROBLEM\n");
        exit(-1);
    }
    int x;
    for(x = 0; x  < d.ncols; x++) {
        d.data[x] = (double*)malloc(d.ncols*sizeof(double));
        if (d.data[x] == NULL) {
            fprintf(stderr, "MEMORY ALLOCATION PROBLEM\n");
            exit(-1);
        }
        (*d.data[x]) = 0;
    }

    return d;
}

/*! \brief Function to make and add Detector

@param z0 z-axis position of detector
@param xmin Smallest x value to detect
@param xmax Largest x value to detect
@param ymin Smallest y value to detect
@param ymax Largest y value to detect
@param xres Number of pixels along x axis
@param yres Number of pixels along y axis
@param num_particles Total number of particles being emitted
@param filename Name of output file of detector (should end in .txt)
@param s Scene to add Detector to
*/
Detector* addDetector(double z0, double xmin, double xmax, double ymin, double ymax, double xres,
    double yres, double num_particles, char* filename, Scene* s) {
    if (s->num_d >= MAX_DETECTOR-1) {
        fprintf(stderr,"TOO MANY DETECTORS IN SCENE");
        exit(-1);
    }
    s->d[s->num_d] = makeDetector(z0,xmin,xmax,ymin,ymax,xres,yres,num_particles,filename);
    s->num_d++;
    return &s->d[s->num_d-1];
}

/*! \brief Function to compute time of first collision for a Detector.

@param p Particle to consider
@param d Detector to consider

@return Time until the propogation or -1 if particle will not hit detector
*/
double getTimeOfFirstCollisionDetector(Particle p, Detector d) {
    double t = (d.z0-p._z)/p._vz;
    if (t <= 0)
        return -1;
    Particle p2 = copyMoveParticleT(t,p);
    if (p2._x > d.xmax || p2._x < d.xmin || p2._y > d.ymax || p2._y < d.ymin)
        return -1;
    return t;
}

/*! \brief Function to raytrace Detector

@param p Pointer to particle to be traced
@param d Detector to be traced
*/
void traceNeutronDetector(Particle* p, Detector d) {
    double t = getTimeOfFirstCollisionDetector(*p, d);
    if (t < 0)
        return;
    moveParticleT(t,p);
    d.data[(int)floor((p->_x-d.xmin)/d.xstep)][(int)floor((p->_y-d.ymin)/d.ystep)] += p->w;
    (*d.num_count) += p->w;
}

/*! \brief Function to finalize detector

Will write data and free data array.

@param d Detector to finalize
*/
void finishDetector(Detector d) {
    int x,y;
    if (d.filename != "") {
        FILE *file;
        file = fopen(d.filename,"w");

        double intensity = (*d.num_count);
        fprintf(file, "#I=%e I_ERR=%e xmin=%f xmax=%f ymin=%f ymax=%f ncols=%i nrows=%i\n",
            intensity, sqrt(intensity/d.num_particles), d.xmin, d.xmax, d.ymin, d.ymax, d.ncols, d.nrows); //FIXME: check I_ERR sqrt(I/num_particles)

        //Write data
        for (x=0; x < d.ncols; x++) {
            for (y=0; y < d.nrows; y++)
                fprintf(file, "%e ", d.data[x][y]);
            fprintf(file, "\n");
        }
        fclose(file);
    }
    for (x=0; x < d.ncols; x++)
        free(d.data[x]);
    free(d.data);
    free(d.num_count);
}

/** @} */ //end of detectorgroup

/////////////////////////////////////
// Geometry Types
/////////////////////////////////////

/////////////////////////////////////
// Disks
/////////////////////////////////////

/** @defgroup diskgroup Disk
    Contains code related to Disks
    @{
*/

/*! \brief Function for creating a Disk structure

@param z0 z-axis position of Disk
@param r0 Inner radius of doughnut
@param r1 Outer radius of doughnut

@see Disk
*/
Disk makeDisk(double z0, double r0, double r1) {
    Disk d;

    d.r0 = r0;
    d.z0 = z0;
    d.r1 = r1;

    return d;
}

/*! \brief Function for making and adding Disk to Scene

@param z0 z-axis position of Disk
@param r0 Inner radius of doughnut
@param r1 Outer radius of doughnut
@param s Scene to add Disk to

@see Disk
*/
Disk* addDisk(double z0, double r0, double r1, Scene* s) {
    if (s->num_di >= MAX_DISK-1) {
        fprintf(stderr,"TOO MANY DISKS IN SCENE");
        exit(-1);
    }
    s->di[s->num_di] = makeDisk(z0, r0, r1);
    s->num_di++;
    return &s->di[s->num_di -1];
}

/*! \brief Function to compute time of first collision for a disk

@param p Particle to consider
@param d Disk to consider
@return Time until the propogation or -1 if particle will not hit disk
*/
double getTimeOfFirstCollisionDisk(Particle p, Disk d) {
    double tz = (d.z0-p._z)/p._vz;
    if (tz <= 0)
        return -1;
    Particle p2 = copyMoveParticleT(tz, p);
    double rp = sqrt(p2._x*p2._x+p2._y*p2._y);
    if (rp > d.r0 && rp < d.r1 && fabs(p2._z-d.z0) < 1e-11)
        return (d.z0-p._z)/p._vz;
    return -1;
}

/*! \brief Function to raytrace Disks

@param p Pointer to particle to be traced
@param d Disk to be traced
*/
void traceNeutronDisk(Particle* p, Disk d) {
    double t = getTimeOfFirstCollisionDisk(*p, d);

    if (t <= 0)
        return;

    moveParticleT(t, p);
    //absorbParticle(p); //Disk will only be used to propagate neutrons somewhere
}

/** @} */ //end of diskgroup

/////////////////////////////////////
// Z-Axis Symetric Conic Sections
/////////////////////////////////////

/** @defgroup conicgroup ConicSurf
    Contains code related to ConicSurfs
    @{
*/

/*! \brief Function to return radius of ConicSurf at a z-axis position.

Will return radius even if z is outside the bounds of zs and ze
for the particular ConicSurf.

@param z z-axis position to compute radius
@param s ConicSurf to compute radius of
*/
double rConic(double z, ConicSurf s) {
    return sqrt(s.k1+s.k2*z+s.k3*z*z);
}

/*! \brief Function for generating Hyperboloid ConicSurf.

@param f1 z position of focus closest to actual mirror surface
@param f2 z position of focus furthest from actual mirror surface
@param p A Point on the actual surface of the mirror
@param zstart z position of start of mirror surface
@param zend z position of end of mirror surface
@param m m value for reflectivity of the surface

@see ConicSurf
*/
ConicSurf makeHyperboloid(double f1, double f2, Point p,
   double zstart, double zend, double m) {
    ConicSurf s;
    s.zs = zstart;
    s.ze = zend;

    double r2 = p.x*p.x+p.y*p.y;
    double c = (f1-f2)/2;

    double u = p.z+c-f1;
    double a = sqrt(((u*u+c*c+r2)-sqrt(pow(u*u+c*c+r2,2)-4*c*c*u*u))/2);

    s.k3 = c*c/(a*a)-1;
    s.k2 = 2*s.k3*(c-f1);
    s.k1 = (s.k3)*(c-f1)*(c-f1)-c*c+a*a;

    s.m = m;
    s.f1 = f1;
    s.f2 = f2;
    s.a = a;
    s.c = c;

    s.type = HYPER;

    #if REC_MAX_GA
    s.max_ga = -1;
    s.max_ga_z0 = -1;
    #endif

    return s;
}

/*! \brief Function for generating Ellipsoid ConicSurf.

@param f1 z position of focus closest to actual mirror surface
@param f2 z position of focus furthest from actual mirror surface
@param p A Point on the actual surface of the mirror
@param zstart z position of start of mirror surface
@param zend z position of end of mirror surface
@param m m value for reflectivity of the surface

@see ConicSurf
*/
ConicSurf makeEllipsoid(double f1, double f2, Point p,
    double zstart, double zend, double m, int doubleReflections) {
    ConicSurf s;
    s.zs = zstart;
    s.ze = zend;
    s.doubleReflections = doubleReflections;
    double r2 = p.x*p.x+p.y*p.y;
    double c = (f1-f2)/2;

    double u = p.z+c-f1;
    double a = sqrt(((u*u+c*c+r2)+sqrt(pow(u*u+c*c+r2,2)-4*c*c*u*u))/2);

    s.k3 = c*c/(a*a)-1;
    s.k2 = 2*s.k3*(c-f1);
    s.k1 = (s.k3)*(c-f1)*(c-f1)-c*c+a*a;

    s.m = m;
    s.f1 = f1;
    s.f2 = f2;
    s.a = a;
    s.c = c;

    s.type = ELLIP;

    #if REC_MAX_GA
    s.max_ga = -1;
    s.max_ga_z0 = -1;
    #endif

    return s;
}


/*! \brief Function for generating Flat Ellipse with symmetry along the vertical y.

@param f1 z position of focus closest to actual mirror surface
@param f2 z position of focus furthest from actual mirror surface
@param b the short half axis of the ellipse, positive for translational symmetry along y, negative for translational symmetry along x
@param p A Point on the actual surface of the mirror
@param zstart z position of start of mirror surface
@param zend z position of end of mirror surface
@param ll the left/lower limit of the ellipse surface (in either y or x for b > 0 or b < 0
@param rl the right/upper limit of the ellipse surface (in either y or x for b > 0 or b < 0
@param m m value for reflectivity of the surface


@see ConicSurf
*/
FlatSurf makeFlatEllipse(
    double f1,
    double f2,
    Point p,
    double zstart,
    double zend,
    double ll,
    double rl,
    double m,
    int doubleReflections) {
        FlatSurf s;
        s.zs = zstart;
        s.ze = zend;
        s.doubleReflections = doubleReflections;
        s.ll = ll;
        s.rl = rl;
        double r2 = p.x*p.x + p.y*p.y;
        double c = (f1-f2)/2;
        double u = p.z+c-f1;
        double a = sqrt(((u*u+c*c+r2)+sqrt(pow(u*u+c*c+r2,2)-4*c*c*u*u))/2);

        s.k3 = c*c/(a*a)-1;
        s.k2 = 2*s.k3*(c-f1);
        s.k1 = (s.k3)*(c-f1)*(c-f1)-c*c+a*a;

        s.m = m;
        s.f1 = f1;
        s.f2 = f2;
        s.a = a;
        s.c = c;

        //s.type = ELLIP;

        #if REC_MAX_GA
        s.max_ga = -1;
        s.max_ga_z0 = -1;
        #endif

        return s;
}

/*! \brief Function for generating Paraboloid ConicSurf.

@param f z position of focus closest to actual mirror surface
@param p A Point on the actual surface of the mirror
@param zstart z position of start of mirror surface
@param zend z position of end of mirror surface
@param m m value for reflectivity of the surface

@see ConicSurf
*/
ConicSurf makeParaboloid(double f, Point p, double zstart,
    double zend, double m, int doubleReflections) {

    ConicSurf s;
    s.zs = zstart;
    s.ze = zend;
    s.doubleReflections = doubleReflections;
    double r2 = p.x*p.x+p.y*p.y;
    double a = (-(p.z-f)+sign(p.z-f)*sqrt((p.z-f)*(p.z-f)+r2))/2;

    s.k3 = 0.0;
    s.k2 = 4*a;
    s.k1 = s.k2*(a-f);

    s.m = m;
    s.f1 = f;
    s.a = a;

    s.type = PARA;

    #if REC_MAX_GA
    s.max_ga = -1;
    s.max_ga_z0 = -1;
    #endif

    return s;
}

/*! \brief Function for generating Flat Parabola for FlatSurf.

@param f z position of focus closest to actual mirror surface
@param p A Point on the actual surface of the mirror, putting one of x or y to 0 results in the surface being parallel to said coordinate
@param zstart z position of start of mirror surface
@param zend z position of end of mirror surface
@param ll the left/lower limit of the ellipse surface (in either y or x for b > 0 or b < 0
@param rl the right/upper limit of the ellipse surface (in either y or x for b > 0 or b < 0
@param m m value for reflectivity of the surface
@param doubleReflections wether double reflections are allowed

@see FlatSurf
*/
FlatSurf makeFlatparbola(
    double f,
    Point p,
    double zstart,
    double zend,
    double ll,
    double rl,
    double m,
    int doubleReflections) {

    FlatSurf s;
    s.zs = zstart;
    s.ze = zend;
    s.doubleReflections = doubleReflections;
    double r2 = p.x*p.x+p.y*p.y;
    double a = (-(p.z-f)+sign(p.z-f)*sqrt((p.z-f)*(p.z-f)+r2))/2;

    s.k3 = 0.0;
    s.k2 = 4*a;
    s.k1 = s.k2*(a-f);

    s.m = m;
    s.f1 = f;
    s.a = a;
    s.ll = ll;
    s.rl = rl;
    //s.type = PARA;

    #if REC_MAX_GA
    s.max_ga = -1;
    s.max_ga_z0 = -1;
    #endif

    return s;
}

/*! \brief Function for generating and adding Paraboloid ConicSurf.

@param f1 z position of focus closest to actual mirror surface
@param p A Point on the actual surface of the mirror
@param zstart z position of start of mirror surface
@param zend z position of end of mirror surface
@param m m value for reflectivity of the surface
@param s Scene to add Paraboloid to

@see ConicSurf
*/
ConicSurf* addParaboloid(double f1, Point p, double zstart, double zend,
    double m, int doubleReflections, Scene* s) {
    if (s->num_c >= MAX_CONICSURF-1) {
        fprintf(stderr,"TOO MANY CONICSURF IN SCENE");
        exit(-1);
    }
    s->c[s->num_c] = makeParaboloid(f1,p,zstart,zend,m, doubleReflections);
    s->num_c++;
    return &s->c[s->num_c-1];
}

/*! \brief Function for generating and adding a flat Parabolic FlatSurf.

@param f1 z position of focus closest to actual mirror surface
@param p A Point on the actual surface of the mirror
@param zstart z position of start of mirror surface
@param zend z position of end of mirror surface
@param ll the lower bound of the mirror
@param rl the upper bound of the mirror
@param m m value for reflectivity of the surface
@param s Scene to add Ellipsoid to
@param doubleReflections wether double reflections can occur
@see FlatSurf
*/
FlatSurf* addFlatParabola(
    double f,
    Point p,
    double zstart,
    double zend,
    double ll,
    double rl,
    double m,
    Scene* s,
    int doubleReflections) {
    if (s->num_f >= MAX_FLATSURF-1) {
        fprintf(stderr,"TOO MANY FLATSURF IN SCENE");
        exit(-1);
    }
    s->f[s->num_f] = makeFlatparbola(f,p,zstart,zend,ll,rl,m,doubleReflections);
    s->num_f++;
    return &s->f[s->num_f-1];
}

/*! \brief Function for generating and adding Hyperboloid ConicSurf.

@param f1 z position of focus closest to actual mirror surface
@param f2 z position of focus furthest from actual mirror surface
@param p A Point on the actual surface of the mirror
@param zstart z position of start of mirror surface
@param zend z position of end of mirror surface
@param m m value for reflectivity of the surface
@param s Scene to add Hyperboloid to

@see ConicSurf
*/
ConicSurf* addHyperboloid(double f1, double f2, Point p, double zstart,
    double zend, double m, Scene* s) {
    if (s->num_c >= MAX_CONICSURF-1) {
        fprintf(stderr,"TOO MANY CONICSURF IN SCENE");
        exit(-1);
    }
    s->c[s->num_c] = makeHyperboloid(f1,f2,p,zstart,zend,m);
    s->num_c++;
    return &s->c[s->num_c-1];
}

/*! \brief Function for generating and adding Ellipsoid ConicSurf.

@param f1 z position of focus closest to actual mirror surface
@param f2 z position of focus furthest from actual mirror surface
@param p A Point on the actual surface of the mirror
@param zstart z position of start of mirror surface
@param zend z position of end of mirror surface
@param m m value for reflectivity of the surface
@param s Scene to add Ellipsoid to

@see ConicSurf
*/
ConicSurf* addEllipsoid(double f1, double f2, Point p, double zstart,
    double zend, double m, int doubleReflections, Scene* s) {
    if (s->num_c >= MAX_CONICSURF-1) {
        fprintf(stderr,"TOO MANY CONICSURF IN SCENE");
        exit(-1);
    }
    s->c[s->num_c] = makeEllipsoid(f1,f2,p,zstart,zend,m,doubleReflections);
    s->num_c++;
    return &s->c[s->num_c-1];
}

/*! \brief Function for generating and adding a flat Ellipse FlatSurf.

@param f1 z position of focus closest to actual mirror surface
@param f2 z position of focus furthest from actual mirror surface
@param p A Point on the actual surface of the mirror
@param zstart z position of start of mirror surface
@param zend z position of end of mirror surface
@param m m value for reflectivity of the surface
@param s Scene to add Ellipsoid to

@see ConicSurf
*/
FlatSurf* addFlatEllipse(
    double f1,
    double f2,
    Point p,
    double zstart,
    double zend,
    double ll,
    double rl,
    double m,
    int doubleReflections,
    Scene* s) {
    if (s->num_f >= MAX_FLATSURF-1) {
        fprintf(stderr,"TOO MANY FLATSURF IN SCENE");
        exit(-1);
    }
    s->f[s->num_f] = makeFlatEllipse(f1,f2,p,zstart,zend,ll,rl,m,doubleReflections);
    s->num_f++;
    return &s->f[s->num_f-1];
}
//!TODO
double getGrazeAngleConic(Particle p, ConicSurf s) {
    /*
    double v = sqrt(dotVec(getParticleVel(p),getParticleVel(p)));
    double vn = dotVec(getParticleVel(p),n);
    return fabs(acos(vn/v)) - M_PI/2;
    */
}

/*! \brief Function for returning normal vector of ConicSurf at Point p

Will compute vector even if p is not on surface.
MAKE SURE p IS ON SURFACE

@param p Point to compute normal vector
@param s ConicSurf to compute normal vector of
*/
Vec getNormConic(Point p, ConicSurf s) {
    double det = s.k2*s.k2+4*s.k3*(p.x*p.x+p.y*p.y-s.k1);
    if (det <= 0.){

        return makeVec(-p.x/sqrt(p.x*p.x + p.y*p.y),-p.y/(p.x*p.x + p.y*p.y),0);
    }
    double den = sqrt(det);
    double nx = -2*p.x/den;
    double ny = -2*p.y/den;
    double nz = sign(2*s.k3*p.z+s.k2);
    double n = sqrt(nx*nx+ny*ny+nz*nz);
    //printf("%f,%f,%f \n", nx/n,ny/n,nz/n);
    return makeVec(nx/n,ny/n,nz/n);
}

/*! \brief Function for returning normal vector of FlatSurf at Point p

Will compute vector even if p is not on surface.
MAKE SURE p IS ON SURFACE

@param p Point to compute normal vector
@param s FlatSurf to compute normal vector of; for s.b > 0 surface posseses translation symmetry along y direction
*/
Vec getNormFlat(Point p, FlatSurf s) {
    double r;
    //if(s.b > 0){
    r = p.x;
    //else{
    //    r = p.y;
    //};
    double den;
    double det = s.k2*s.k2+4*s.k3*(r*r-s.k1);
    if (det > 0){
        den = sqrt(det);
    }
    else{
    return makeVec(-1,0,0); // if the neutron hits the apex of the ellipse we run into a divide by zero problem
    }
    double nx = -2*p.x/den;
    double ny = 0;
    double nz = sign(2*s.k3*p.z+s.k2);
    double n = sqrt(nx*nx+ny*ny+nz*nz);
    return makeVec(nx/n,ny/n,nz/n);
}

/*! \brief Function to compute time of first collision for a ConicSurf

@param p Particle to consider
@param s ConicSurf to consider
@return Time until the propogation or -1 if particle will not hit surface
*/
double getTimeOfFirstCollisionConic(Particle p, ConicSurf s) {
    double tz = (s.zs-p._z)/p._vz;
    if (tz < 0) {
       tz = 0;
       if (p._z > s.ze)
            return -1;
    }

    Particle p2 = copyMoveParticleT(tz,p);

    double A = p2._vx*p2._vx+p2._vy*p2._vy-s.k3*p2._vz*p2._vz;
    double B = 2*(p2._vx*p2._x+p2._vy*p2._y-s.k3*p2._vz*p2._z)-s.k2*p2._vz;
    double C = p2._x*p2._x+p2._y*p2._y-s.k3*p2._z*p2._z-s.k2*p2._z-s.k1;

    double t = solveQuad(A,B,C);

    if (t <= 0 || p2._vz*t+p2._z > s.ze || p2._vz*t+p2._z < s.zs)
        return -1;
    return t+tz;
}

/*! \brief Function to compute time of first collision for a FlatSurf

@param p Particle to consider
@param s FlatSurf to consider
@return Time until the propogation or -1 if particle will not hit surface
*/
//TODO
double getTimeOfFirstCollisionFlat(Particle p, FlatSurf s) {
    double tz = (s.zs-p._z)/p._vz;
    if (tz < 0) {
       tz = 0;
       if (p._z > s.ze)
            return -1;
    }

    Particle p2 = copyMoveParticleT(tz,p);
    double vs = 0;//the vector important for calculating the intersection with the ellipse
    double s0 = 0;
    double vt = 0;//the other component only important for testing whether the mirror is hit
    double t0 = 0;
    //if(s.b > 0){
    vs = p2._vx;
    s0 = p2._x;
    vt = p2._vy;
    t0 = p2._y;

    //}
    /*else{
    vs = p2._vy;
    s0 = p2._y;
    vt = p2._vx;
    t0 = p2._x;
    };
    */
    double A = vs*vs-s.k3*p2._vz*p2._vz;
    double B = 2*(vs*s0-s.k3*p2._vz*p2._z)-s.k2*p2._vz;
    double C = s0*s0-s.k3*p2._z*p2._z-s.k2*p2._z-s.k1;

    double t = solveQuad(A,B,C);

    if (t <= 0 || p2._vz*t+p2._z > s.ze || p2._vz*t+p2._z < s.zs||vt*t+t0 < s.ll||vt*t +t0 > s.rl)
        return -1;
    return t+tz;
}

/*! \brief Function to handle supermirror reflectivity copied from mcstas.
@note Uses only m-value for calculating reflectivity curve TODO more sophisticated formulae in the future

@param q k_i - k_f momentum transfer of the neutron at the super mirror surface
@param m supermirror m-value
@param R_0 low angle reflectivity
@param Q_c critical momentum transfer of the super mirror

@return p weight reduction of the neutron for further simulation
*/
double calcSupermirrorReflectivity(double q, double m, double R_0, double Q_c){
    double arg;
    double beta = 0;
    double alpha = 0;
    double W = 0;
    double weight = 1.0; //neutron weight to be transformed
    q = fabs(q);
    if (m >= 10){
        weight = 1.0;
        return weight;
    }
    if (W==0 && alpha==0) {
      m=m*0.9853+0.1978;
      W=-0.0002*m+0.0022;
      alpha=0.2304*m+5.0944;
      beta=-7.6251*m+68.1137;
      if (m<=3) {
	    alpha=m;
	    beta=0;
        }
    }
    arg = W > 0 ? (q - m*Q_c)/W : 11;
    if (arg > 10 || m <= 0 || Q_c <=0 || R_0 <= 0) {
      weight = 0.0;
      return weight;
    }

    if (m < 1) { Q_c *= m; m=1; }

    if(q <= Q_c) {
      weight = R_0;
      return weight;
    }


    weight = R_0*0.5*(1 - tanh(arg))*(1 - alpha*(q - Q_c) + beta*(q - Q_c)*(q - Q_c));
    return weight;
}

/*! \brief Function to handle reflection of neutron for a ConicSurf.

@note Uses step function for reflectivity

@warning Make sure particle has been moved to surface of mirror
before computing reflection

@param p Pointer of particle to reflect
@param s ConicSurf to use

@return Value of critical angle of the neutron or -1 if neutron is absorbed

@see traceNeutronConic()
*/
double reflectNeutronConic(Particle* p, ConicSurf s) {//TODO add super mirror reflectivity an passing neutrons
    Vec n = getNormConic(getParticlePos(*p),s);
    Vec pv = getParticleVel(*p);

    double v = getMagVec(pv);
    double vn = dotVec(pv,n);
    double weight = calcSupermirrorReflectivity(V2Q_conic*vn*2, s.m, 1.0, 0.0218);

    //Hitting shell from outside

    if (vn > 0 && !s.doubleReflections) {
        absorbParticle(p);
        return -1;
    }



    double ga = fabs(acos(vn/v)) - M_PI/2;
    double gc = 6.84459399932*s.m/v;

    if (weight <= 0) {
        printf("weight <0");
        absorbParticle(p);
        return -1;
    }
    else {
        p->_vx = p->_vx-2*vn*n.x;
        p->_vy = p->_vy-2*vn*n.y;
        p->_vz = p->_vz-2*vn*n.z;
        p->w *= weight;
    }
    return ga;
}

/*! \brief Function to handle reflection of neutron for a FlatSurf.

@note Uses step function for reflectivity

@warning Make sure particle has been moved to surface of mirror
before computing reflection

@param p Pointer of particle to reflect
@param s FlatSurf to use

@return Value of critical angle of the neutron or -1 if neutron is absorbed

@see traceNeutronConic()
*/
double reflectNeutronFlat(Particle* p, FlatSurf s) {
    Vec n = getNormFlat(getParticlePos(*p),s);
    Vec pv = getParticleVel(*p);

    double v = getMagVec(pv);
    double vn = dotVec(pv,n);
    //printf("before %f \n", p->w);
    //Hitting shell from outside For FlatSurface this has to be checked
    // make it able to reflect from the outside
    double ga = fabs(acos(vn/v)) - M_PI/2;
    double gc = 6.84459399932*s.m/v;


    double weight = 0;
    weight = calcSupermirrorReflectivity(V2Q_conic*2*vn, s.m, 0.995, 0.0218);

    if (vn > 0 && !s.doubleReflections) {
        absorbParticle(p);
        return -1;
    }


    if (weight < 0) {
        printf("this happens?");
        absorbParticle(p);
        return -1;
    } else {
        if (getRandom() <= weight){//to be updated to use the mcstas random function or quasoi deterministic model
            //printf("oh a reflections\n");
            p->_vx = p->_vx-2*vn*n.x;
            p->_vy = p->_vy-2*vn*n.y;
            p->_vz = p->_vz-2*vn*n.z;
        }
        else{
            //printf("no reflection \n");
            p->silicon *= -1;
            return ga;
        }
        //p->w *= weight;
        //printf("weight in reflectNeutron %f %f\n", p->w, weight);
    }
    return ga;
}


/*! \brief Function to handle raytracing of neutron for a ConicSurf.

@param p Pointer of particle to reflect
@param c ConicSurf to use

*/
void traceNeutronConic(Particle* p, ConicSurf c) {
    double t = getTimeOfFirstCollisionConic(*p, c);
    if (t < 0)
        return;
    else {
        moveParticleT(t, p);
        double ga = reflectNeutronConic(p, c);
#if REC_MAX_GA
        if (ga > c.max_ga) {
            c.max_ga = ga;
            c.max_ga_z0 = p->_z;
        }
#endif
    }
}

/*! \brief Function to handle raytracing of neutron for a FlatSurf.

@param p Pointer of particle to reflect
@param f FlatSurf to use

*/
void traceNeutronFlat(Particle* p, FlatSurf f) {
    double t = getTimeOfFirstCollisionFlat(*p, f);
    if (t < 0)
        return;
    else {

        moveParticleT(t, p);

        //printf("weight before reflect %f", p->w);
        double ga = reflectNeutronFlat(p, f);
        //printf("weight after reflect %f\n", p->w);
#if REC_MAX_GA
        if (ga > f.max_ga) {
            f.max_ga = ga;
            f.max_ga_z0 = p->_z;
        }
#endif
    }
}
/** @} */ //end of conicgroup
/////////////////////////////////////
// Scene Functions
/////////////////////////////////////
/** @ingroup simgroup
    @{
*/
enum GEO {
    NONE,
    DETECTOR,
    DISK,
    CONIC,
    FLAT
};

//! Function to generate an empty Scene
Scene makeScene() {
    Scene s;
    s.num_f = 0;
    s.num_c = 0;
    s.num_di = 0;
    s.num_d = 0;

    s.traceNeutronFlat = traceNeutronFlat;
    s.traceNeutronConic = traceNeutronConic;
    s.traceNeutronDisk = traceNeutronDisk;
    s.traceNeutronDetector = traceNeutronDetector;

    return s;
}

//! Function to init simulation items
/*! Should be called after all items
have been added to scene but before
neutrons are traced.

@param s Pointer of Scene to init
*/
void initSimulation(Scene* s) {
    //
}

/*! \brief Function to raytrace single neutron through geometries specified by d, di and c.

@param p Pointer of particle to trace
@param s Scene to trace
*/
void traceSingleNeutron(Particle* p, Scene s) {

    int contact = 1;
    do {
        double t;
        enum  GEO type = NONE;
        int index = -1;
        int i;

        for (i = 0; i < s.num_c; i++) {
            double t2 = getTimeOfFirstCollisionConic(*p,s.c[i]);

            if (t2 <= 0)
                continue;
            if (index == -1 || t2 < t) {
                type = CONIC;
                index = i;
                t = t2;
            }
        }

        for (i = 0; i < s.num_f; i++) {
            double t2 = getTimeOfFirstCollisionFlat(*p,s.f[i]);

            if (t2 <= 0)
                continue;
            if (index == -1 || t2 < t) {
                type = FLAT;
                index = i;
                t = t2;
            }
        }

        for (i = 0; i < s.num_di; i++)  {
            double t2 = getTimeOfFirstCollisionDisk(*p,s.di[i]);

            if (t2 <= 0)
                continue;
            else if (index == -1 || t2 < t) {
                type = DISK;
                index = i;
                t = t2;
            }
        }

        for (i = 0; i < s.num_d; i++) {
            double t2 = getTimeOfFirstCollisionDetector(*p,s.d[i]);

            if (t2 <= 0)
                continue;
            else if (index == -1 || t2 < t) {
                type = DETECTOR;
                index = i;
                t = t2;
            }
        }

        switch (type) {
            case DETECTOR:
                s.traceNeutronDetector(p, s.d[index]);
                break;
            case FLAT:
                s.traceNeutronFlat(p, s.f[index]);
                break;
            case DISK:
                s.traceNeutronDisk(p, s.di[index]);
                break;
            case CONIC:
                s.traceNeutronConic(p, s.c[index]);
                break;
            default:
                contact = 0;
                break;
        }
    } while (contact && !p->absorb);

}

//!Finishes tracing the scene
/*! This function should be called after all of the
particles have been raytraced.

@param s Pointer of Scene to finish tracing
*/
void finishSimulation(Scene* s) {
    int i;

    //Finish Detectors
    for (i=0; i < s->num_d; i++)
        finishDetector(s->d[i]);
}

/** @} */ //end of ingroup simgroup

#endif

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
double * get_r_at_z0(int number, double z_0, double r_0, double LStart, double LEnd, double lStart, double lEnd) {
    /*
    number: how many mirrors,
    z_0: z-position of returned points on mirror
    r_0: r-distance of returned points on mirror
    LStart: Position of the first focal point
    LEnd: Postion of the second focal point
    lStart: Beginning of the mirror
    lEnd: End of the mirror
    */
    int n = number;
    double *r_z0s = malloc(n*sizeof(double_t)); /* n is an array of 10 integers */
	r_z0s[0] = r_0;
    //helper variables as in conic_finite_mirror.h and explained in swissneutronics_??berlegungen
    double k1;
    double k2;
    double k3;
    double c;
    double u;
    double a;
    double r_lEnd;
    double r_lStart;
    //initial mirror is calculated from the initial point z0, r0
    c = (LEnd - LStart)/2;
    u = (z_0 + c - LEnd);
    a = sqrt((u*u+c*c+r_0*r_0+sqrt(pow(u*u+c*c+r_0*r_0, 2)-4*c*c*u*u))/2);
    k3 = c*c/(a*a)-1;
    k2 = 2*k3*(c-LEnd);
    k1 = k3*(c-LEnd)*(c-LEnd)-c*c+a*a;
    r_lEnd = sqrt(k1+ k2*lEnd + k3*lEnd*lEnd);
    r_lStart = r_lEnd*(lStart-LStart)/(lEnd-LStart);
    r_z0s[0] = r_0;
	//next mirror will be calculated with the point on the surface being lStart, r_lStart
	for( int k = 1; k < number;++k){
        c = (LEnd - LStart)/2;
        u = (lStart + c - LEnd);
        a = sqrt((u*u+c*c+r_lStart*r_lStart+sqrt(pow(u*u+c*c+r_lStart*r_lStart, 2)-4*c*c*u*u))/2);
        k3 = c*c/(a*a)-1;
        k2 = 2*k3*(c-LEnd);
        k1 = k3*(c-LEnd)*(c-LEnd)-c*c+a*a;
        r_lEnd = sqrt(k1+ k2*lEnd + k3*lEnd*lEnd);
        r_lStart = r_lEnd*(lStart-LStart)/(lEnd-LStart);
		if (r_z0s[k-1] > 0)
        {
		r_z0s[k] = sqrt(k1+ k2*z_0 + k3*z_0*z_0);
        }
        else
        {
        r_z0s[k] = -sqrt(k1+ k2*z_0 + k3*z_0*z_0);
        }
	};
   return r_z0s;
}




    //%include "w1_general.h"
    //%include "read_table-lib"
#line 7195 "./ess_extraction.c"

/* Shared user declarations for all components 'Guide_gravity'. */
#line 124 "/usr/share/mcstas/2.6.1/tools/Python/mcrun/../mccodelib/../../../optics/Guide_gravity.comp"
/*****************************************************************************
*
* McStas, neutron ray-tracing package
*         Copyright 1997-2006, All rights reserved
*         Risoe National Laboratory, Roskilde, Denmark
*         Institut Laue Langevin, Grenoble, France
*
* Library: share/ref-lib.h
*
* %Identification
* Written by: Peter Christiansen
* Date: August, 2006
* Origin: RISOE
* Release: McStas 1.10
* Version: $Revision$
*
* Commonly used reflection functions are declared in this file which
* are used by some guide and mirror components.
*
* Depends on read_table-lib
*
* Usage: within SHARE
* %include "ref-lib"
*
****************************************************************************/
/*******************************************************************************
*
* McStas, neutron ray-tracing package
*         Copyright 1997-2002, All rights reserved
*         Risoe National Laboratory, Roskilde, Denmark
*         Institut Laue Langevin, Grenoble, France
*
* Library: share/read_table-lib.h
*
* %Identification
* Written by: EF
* Date: Aug 28, 2002
* Origin: ILL
* Release: McStas 1.6
* Version: $Revision$
*
* This file is to be imported by components that may read data from table files
* It handles some shared functions.
*
* This library may be used directly as an external library. It has no dependency
*
* Usage: within SHARE
* %include "read_table-lib"
*
*******************************************************************************/

#ifndef READ_TABLE_LIB_H
#define READ_TABLE_LIB_H "$Revision$"

#define READ_TABLE_STEPTOL  0.04 /* tolerancy for constant step approx */

#ifndef MC_PATHSEP_C
#ifdef WIN32
#define MC_PATHSEP_C '\\'
#define MC_PATHSEP_S "\\"
#else  /* !WIN32 */
#ifdef MAC
#define MC_PATHSEP_C ':'
#define MC_PATHSEP_S ":"
#else  /* !MAC */
#define MC_PATHSEP_C '/'
#define MC_PATHSEP_S "/"
#endif /* !MAC */
#endif /* !WIN32 */
#endif /* !MC_PATHSEP_C */

#ifndef MCSTAS
#ifdef WIN32
#define MCSTAS "C:\\mcstas\\lib"
#else  /* !WIN32 */
#ifdef MAC
#define MCSTAS ":mcstas:lib" /* ToDo: What to put here? */
#else  /* !MAC */
#define MCSTAS "/usr/local/lib/mcstas"
#endif /* !MAC */
#endif /* !WIN32 */
#endif /* !MCSTAS */

#include <sys/stat.h>
#include <stdio.h>
#include <stdlib.h>

  typedef struct struct_table
  {
    char    filename[1024];
    long    filesize;
    char   *header;  /* text header, e.g. comments */
    double *data;    /* vector { x[0], y[0], ... x[n-1], y[n-1]... } */
    double  min_x;   /* min value of first column */
    double  max_x;   /* max value of first column */
    double  step_x;  /* minimal step value of first column */
    long    rows;    /* number of rows in matrix block */
    long    columns; /* number of columns in matrix block */

    long    begin;   /* start fseek index of block */
    long    end;     /* stop  fseek index of block */
    long    block_number;  /* block index. 0 is catenation of all */
    long    array_length;  /* number of elements in the t_Table array */
    char    monotonic;     /* true when 1st column/vector data is monotonic */
    char    constantstep;  /* true when 1st column/vector data has constant step */
    char    method[32];    /* interpolation method: nearest, linear */
  } t_Table;

/*maximum number of rows to rebin a table = 1M*/
enum { mcread_table_rebin_maxsize = 1000000 };

typedef struct t_Read_table_file_item {
    int ref_count;
    t_Table *table_ref;
} t_Read_table_file_item;

typedef enum enum_Read_table_file_actions {STORE,FIND,GC}  t_Read_table_file_actions;

/* read_table-lib function prototypes */
/* ========================================================================= */

/* 'public' functions */
long     Table_Read              (t_Table *Table, char *File, long block_number);
long     Table_Read_Offset       (t_Table *Table, char *File, long block_number,
                                  long *offset, long max_lines);
long     Table_Read_Offset_Binary(t_Table *Table, char *File, char *Type,
                                  long *Offset, long Rows, long Columns);
long     Table_Rebin(t_Table *Table); /* rebin table with regular 1st column and interpolate all columns 2:end */
long     Table_Info (t_Table Table);
double   Table_Index(t_Table Table,   long i, long j); /* get indexed value */
double   Table_Value(t_Table Table, double X, long j); /* search X in 1st column and return interpolated value in j-column */
t_Table *Table_Read_Array(char *File, long *blocks);
void     Table_Free_Array(t_Table *Table);
long     Table_Info_Array(t_Table *Table);
int      Table_SetElement(t_Table *Table, long i, long j, double value);
long     Table_Init(t_Table *Table, long rows, long columns); /* create a Table */
double   Table_Value2d(t_Table Table, double X, double Y);    /* same as Table_Index with non-integer indices and 2d interpolation */
MCDETECTOR Table_Write(t_Table Table, char*file, char*xl, char*yl, 
           double x1, double x2, double y1, double y2); /* write Table to disk */
void * Table_File_List_Handler(t_Read_table_file_actions action, void *item, void *item_modifier);
t_Table *Table_File_List_find(char *name, int block, int offset);
int Table_File_List_gc(t_Table *tab);
void *Table_File_List_store(t_Table *tab);

#define Table_ParseHeader(header, ...) \
  Table_ParseHeader_backend(header,__VA_ARGS__,NULL);

char **Table_ParseHeader_backend(char *header, ...);

/* private functions */
void Table_Free(t_Table *Table);
long Table_Read_Handle(t_Table *Table, FILE *fid, long block_number, long max_lines, char *name);
static void Table_Stat(t_Table *Table);
double Table_Interp1d(double x, double x1, double y1, double x2, double y2);
double Table_Interp1d_nearest(double x, double x1, double y1, double x2, double y2);
double Table_Interp2d(double x, double y, double x1, double y1, double x2, double y2,
double z11, double z12, double z21, double z22);

#endif

/* end of read_table-lib.h */
/*******************************************************************************
*
* McStas, neutron ray-tracing package
*         Copyright (C) 1997-2009, All rights reserved
*         Risoe National Laboratory, Roskilde, Denmark
*         Institut Laue Langevin, Grenoble, France
*
* Library: share/read_table-lib.c
*
* %Identification
* Written by: EF
* Date: Aug 28, 2002
* Origin: ILL
* Release: McStas CVS_090504
* Version: $Revision: 5052 $
*
* This file is to be imported by components that may read data from table files
* It handles some shared functions. Embedded within instrument in runtime mode.
*
* Usage: within SHARE
* %include "read_table-lib"
*
*******************************************************************************/

#ifndef READ_TABLE_LIB_H
#include "read_table-lib.h"
#endif


/*******************************************************************************
 * void *Table_File_List_Handler(action, item, item_modifier)
 *   ACTION: handle file entries in the read_table-lib file list. If a file is read - it is supposed to be
 *   stored in a list such that we can avoid reading the same file many times.
 *   input  action: FIND, STORE, GC. check if file exists in the list, store an item in the list, or check if it can be garbage collected.
 *   input item: depends on the action.
 *    FIND)  item is a filename, and item_modifier is the block number
 *    STORE) item is the Table to store - item_modifier is ignored
 *    GC)    item is the Table to check. If it has a ref_count >1 then this is simply decremented.
 *   return  depends on the action
 *    FIND)  return a reference to a table+ref_count item if found - NULL otherwise. I.e. NULL means the file has not been read before and must be read again.
 *    STORE) return NULL always
 *    GC)    return NULL if no garbage collection is needed, return an adress to the t_Table which should be garbage collected. 0x1 is returned if
 *           the item is not found in the list
*******************************************************************************/
void * Table_File_List_Handler(t_Read_table_file_actions action, void *item, void *item_modifier){

    /* logic here is Read_Table should include a call to FIND. If found the return value should just be used as
     * if the table had been read from disk. If not found then read the table and STORE.
     * Table_Free should include a call to GC. If this returns non-NULL then we should proceed with freeing the memory
     * associated with the table item - otherwise only decrement the reference counter since there are more references
     * that may need it.*/

    static t_Read_table_file_item read_table_file_list[1024];  
    static int read_table_file_count=0;

    t_Read_table_file_item *tr;
    switch(action){
        case FIND:
            /*interpret data item as a filename, if it is found return a pointer to the table and increment refcount.
             * if not found return the item itself*/
            tr=read_table_file_list;
            while ( tr->table_ref!=NULL ){
                int i=*((int*) item_modifier);
                int j=*( ((int*) item_modifier)+1);
                if ( !strcmp(tr->table_ref->filename,(char *) item) &&
                        tr->table_ref->block_number==i && tr->table_ref->begin==j ){
                    tr->ref_count++;
                    return (void *) tr;
                }
                tr++;
            }
            return NULL;
        case STORE:
            /*find an available slot and store references to table there*/
            tr=&(read_table_file_list[read_table_file_count++]);
            tr->table_ref = ((t_Table *) item);
            tr->ref_count++;
            return NULL;
        case GC:
            /* Should this item be garbage collected (freed) - if so scratch the entry and return the address of the item - 
             * else decrement ref_count and return NULL.
             * A non-NULL return expects the item to actually be freed afterwards.*/
            tr=read_table_file_list;
            while ( tr->table_ref!=NULL ){
                if ( tr->table_ref->data ==((t_Table *)item)->data && 
                        tr->table_ref->block_number == ((t_Table *)item)->block_number){
                    /*matching item found*/
                    if (tr->ref_count>1){
                        /*the item is found and no garbage collection needed*/
                        tr->ref_count--;
                        return NULL;
                    }else{
                        /* The item is found and the reference counter is 1.
                         * This means we should garbage collect. Move remaining list items up one slot,
                         * and return the table for garbage collection by caller*/
                        while (tr->table_ref!=NULL){
                            *tr=*(tr+1);
                            tr++;
                        }
                        read_table_file_count--;
                        return (t_Table *) item;
                    }
                }
                tr++;
            }
            /* item not found, and so should be garbage collected. This could be the case if freeing a
             * Table that has been constructed from code - not read from file. Return 0x1 to flag it for
             * collection.*/
            return (void *) 0x1 ;
    }
}

/* Access functions to the handler*/

/********************************************
 * t_Table *Table_File_List_find(char *name, int block, int offset)
 * input name: filename to search for in the file list
 * input block: data block in the file as each file may contain more than 1 data block.
 * return a ref. to a table if it is found (you may use this pointer and skip reading the file), NULL otherwise (i.e. go ahead and read the file)
*********************************************/
t_Table *Table_File_List_find(char *name, int block, int offset){
    int vars[2]={block,offset};
    t_Read_table_file_item *item = Table_File_List_Handler(FIND,name, vars);
    if (item == NULL){
        return NULL;
    }else{
        return item->table_ref;
    }
}
/********************************************
 * int Table_File_List_gc(t_Table *tab)
 * input tab: the table to check for references.
 * return 0: no garbage collection needed
 *        1: Table's data and header (at least) should be freed.
*********************************************/
int Table_File_List_gc(t_Table *tab){
    void *rval=Table_File_List_Handler(GC,tab,0);
    if (rval==NULL) return 0;
    else return 1;
}


/*****************************************************************************
 * void *Table_File_List_store(t_Table *tab)
 * input tab: pointer to table to store.
 * return None. 
*******************************************************************************/
void *Table_File_List_store(t_Table *tab){
    return Table_File_List_Handler(STORE,tab,0);
}


/*******************************************************************************
* FILE *Open_File(char *name, char *Mode, char *path)
*   ACTION: search for a file and open it. Optionally return the opened path.
*   input   name:  file name from which table should be extracted
*           mode: "r", "w", "a" or any valid fopen mode
*           path:  NULL or a pointer to at least 1024 allocated chars
*   return  initialized file handle or NULL in case of error
*******************************************************************************/

  FILE *Open_File(char *File, const char *Mode, char *Path)
  {
    char path[1024];
    FILE *hfile = NULL;
    
    if (!File || File[0]=='\0')                     return(NULL);
    if (!strcmp(File,"NULL") || !strcmp(File,"0"))  return(NULL);
    
    /* search in current or full path */
    strncpy(path, File, 1024);
    hfile = fopen(path, Mode);
    if(!hfile)
    {
      char dir[1024];

      if (!hfile && mcinstrument_source[0] != '\0' && strlen(mcinstrument_source)) /* search in instrument source location */
      {
        char *path_pos   = NULL;
        /* extract path: searches for last file separator */
        path_pos    = strrchr(mcinstrument_source, MC_PATHSEP_C);  /* last PATHSEP */
        if (path_pos) {
          long path_length = path_pos +1 - mcinstrument_source;  /* from start to path+sep */
          if (path_length) {
            strncpy(dir, mcinstrument_source, path_length);
            dir[path_length] = '\0';
            snprintf(path, 1024, "%s%c%s", dir, MC_PATHSEP_C, File);
            hfile = fopen(path, Mode);
          }
        }
      }
      if (!hfile && mcinstrument_exe[0] != '\0' && strlen(mcinstrument_exe)) /* search in PWD instrument executable location */
      {
        char *path_pos   = NULL;
        /* extract path: searches for last file separator */
        path_pos    = strrchr(mcinstrument_exe, MC_PATHSEP_C);  /* last PATHSEP */
        if (path_pos) {
          long path_length = path_pos +1 - mcinstrument_exe;  /* from start to path+sep */
          if (path_length) {
            strncpy(dir, mcinstrument_exe, path_length);
            dir[path_length] = '\0';
            snprintf(path, 1024, "%s%c%s", dir, MC_PATHSEP_C, File);
            hfile = fopen(path, Mode);
          }
        }
      }
      if (!hfile) /* search in HOME or . */
      {
        strcpy(dir, getenv("HOME") ? getenv("HOME") : ".");
        snprintf(path, 1024, "%s%c%s", dir, MC_PATHSEP_C, File);
        hfile = fopen(path, Mode);
      }
      if (!hfile) /* search in MCSTAS/data */
      {
        strcpy(dir, getenv(FLAVOR_UPPER) ? getenv(FLAVOR_UPPER) : MCSTAS);
        snprintf(path, 1024, "%s%c%s%c%s", dir, MC_PATHSEP_C, "data", MC_PATHSEP_C, File);
        hfile = fopen(path, Mode);
      }
      if (!hfile) /* search in MVCSTAS/contrib */
      {
        strcpy(dir, getenv(FLAVOR_UPPER) ? getenv(FLAVOR_UPPER) : MCSTAS);
        snprintf(path, 1024, "%s%c%s%c%s", dir, MC_PATHSEP_C, "contrib", MC_PATHSEP_C, File);
        hfile = fopen(path, Mode);
      }
      if(!hfile)
      {
        fprintf(stderr, "Error: Could not open input file '%s' (Open_File)\n", File);
        return (NULL);
      }
    }
    if (Path) strncpy(Path, path, 1024);
    return(hfile);
  } /* end Open_File */

/*******************************************************************************
* long Read_Table(t_Table *Table, char *name, int block_number)
*   ACTION: read a single Table from a text file
*   input   Table: pointer to a t_Table structure
*           name:  file name from which table should be extracted
*           block_number: if the file does contain more than one
*                 data block, then indicates which one to get (from index 1)
*                 a 0 value means append/catenate all
*   return  initialized single Table t_Table structure containing data, header, ...
*           number of read elements (-1: error, 0:header only)
* The routine stores any line starting with '#', '%' and ';' into the header
* File is opened, read and closed
* Other lines are interpreted as numerical data, and stored.
* Data block should be a rectangular matrix or vector.
* Data block may be rebinned with Table_Rebin (also sort in ascending order)
*******************************************************************************/
  long Table_Read(t_Table *Table, char *File, long block_number)
  { /* reads all or a single data block from 'file' and returns a Table structure  */
    return(Table_Read_Offset(Table, File, block_number, NULL, 0));
  } /* end Table_Read */

/*******************************************************************************
* long Table_Read_Offset(t_Table *Table, char *name, int block_number, long *offset
*                        long max_rows)
*   ACTION: read a single Table from a text file, starting at offset
*     Same as Table_Read(..) except:
*   input   offset:    pointer to an offset (*offset should be 0 at start)
*           max_rows: max number of data rows to read from file (0 means all)
*   return  initialized single Table t_Table structure containing data, header, ...
*           number of read elements (-1: error, 0:header only)
*           updated *offset position (where end of reading occured)
*******************************************************************************/
  long Table_Read_Offset(t_Table *Table, char *File,
                         long block_number, long *offset,
                         long max_rows)
  { /* reads all/a data block in 'file' and returns a Table structure  */
    FILE *hfile;
    long  nelements=0;
    long  begin=0;
    long  filesize=0;
    char  name[1024];
    char  path[1024];
    struct stat stfile;

    /*Need to be able to store the pointer*/
    if (!Table) return(-1);
    
    //if (offset && *offset) snprintf(name, 1024, "%s@%li", File, *offset);
    //else                   
    strncpy(name, File, 1024);
    if(offset && *offset){
        begin=*offset;
    }
    /* Check if the table has already been read from file.
     * If so just reuse the table, if not (this is flagged by returning NULL
     * set up a new table and read the data into it */
    t_Table *tab_p= Table_File_List_find(name,block_number,begin);
    if ( tab_p!=NULL ){
        /*table was found in the Table_File_List*/
        // printf("Reusing input file '%s' (Table_Read_Offset)\n", name);
        *Table=*tab_p;
        return Table->rows*Table->columns;
    }

    /* open the file */
    hfile = Open_File(File, "r", path);
    if (!hfile) return(-1);
    else {
      MPI_MASTER(
      printf("Opening input file '%s' (Table_Read_Offset)\n", path);
      );
    }
    
    /* read file state */
    stat(path,&stfile); filesize = stfile.st_size;
    if (offset && *offset) fseek(hfile, *offset, SEEK_SET);
    begin     = ftell(hfile);
    
    Table_Init(Table, 0, 0);

    /* read file content and set the Table */
    nelements = Table_Read_Handle(Table, hfile, block_number, max_rows, name);
    Table->begin = begin;
    Table->end   = ftell(hfile);
    Table->filesize = (filesize>0 ? filesize : 0);
    Table_Stat(Table);
    
    Table_File_List_store(Table);

    if (offset) *offset=Table->end;
    fclose(hfile);
    return(nelements);

  } /* end Table_Read_Offset */

/*******************************************************************************
* long Table_Read_Offset_Binary(t_Table *Table, char *File, char *type,
*                               long *offset, long rows, long columns)
*   ACTION: read a single Table from a binary file, starting at offset
*     Same as Table_Read_Offset(..) except that it handles binary files.
*   input   type: may be "float"/NULL or "double"
*           offset: pointer to an offset (*offset should be 0 at start)
*           rows   : number of rows (0 means read all)
*           columns: number of columns
*   return  initialized single Table t_Table structure containing data, header, ...
*           number of read elements (-1: error, 0:header only)
*           updated *offset position (where end of reading occured)
*******************************************************************************/
  long Table_Read_Offset_Binary(t_Table *Table, char *File, char *type,
                                long *offset, long rows, long columns)
  { /* reads all/a data block in binary 'file' and returns a Table structure  */
    long    nelements, sizeofelement;
    long    filesize;
    FILE   *hfile;
    char    path[1024];
    struct stat stfile;
    double *data;
    long    i;
    long    begin;

    if (!Table) return(-1);

    Table_Init(Table, 0, 0);
    
    /* open the file */
    hfile = Open_File(File, "r", path);
    if (!hfile) return(-1);
    else {
      MPI_MASTER(
      printf("Opening input file '%s' (Table_Read, Binary)\n", path);
      );
    }
    
    /* read file state */
    stat(File,&stfile);
    filesize = stfile.st_size;
    Table->filesize=filesize;
    
    /* read file content */
    if (type && !strcmp(type,"double")) sizeofelement = sizeof(double);
    else  sizeofelement = sizeof(float);
    if (offset && *offset) fseek(hfile, *offset, SEEK_SET);
    begin     = ftell(hfile);
    if (rows && filesize > sizeofelement*columns*rows)
      nelements = columns*rows;
    else nelements = (long)(filesize/sizeofelement);
    if (!nelements || filesize <= *offset) return(0);
    data    = (double*)malloc(nelements*sizeofelement);
    if (!data) {
      fprintf(stderr,"Error: allocating %ld elements for %s file '%s'. Too big (Table_Read_Offset_Binary).\n", nelements, type, File);
      exit(-1);
    }
    nelements = fread(data, sizeofelement, nelements, hfile);

    if (!data || !nelements)
    {
      fprintf(stderr,"Error: reading %ld elements from %s file '%s' (Table_Read_Offset_Binary)\n", nelements, type, File);
      exit(-1);
    }
    Table->begin   = begin;
    Table->end     = ftell(hfile);
    if (offset) *offset=Table->end;
    fclose(hfile);
    data = (double*)realloc(data, (double)nelements*sizeofelement);
    /* copy file data into Table */
    if (type && !strcmp(type,"double")) Table->data = data;
    else {
      float  *s;
      double *dataf;
      s     = (float*)data;
      dataf = (double*)malloc(sizeof(double)*nelements);
      for (i=0; i<nelements; i++)
        dataf[i]=s[i];
      free(data);
      Table->data = dataf;
    }
    strncpy(Table->filename, File, 1024);
    Table->rows    = nelements/columns;
    Table->columns = columns;
    Table->array_length = 1;
    Table->block_number = 1;

    Table_Stat(Table);

    return(nelements);
  } /* end Table_Read_Offset_Binary */

/*******************************************************************************
* long Table_Read_Handle(t_Table *Table, FILE *fid, int block_number, long max_rows, char *name)
*   ACTION: read a single Table from a text file handle (private)
*   input   Table:pointer to a t_Table structure
*           fid:  pointer to FILE handle
*           block_number: if the file does contain more than one
*                 data block, then indicates which one to get (from index 1)
*                 a 0 value means append/catenate all
*           max_rows: if non 0, only reads that number of lines
*   return  initialized single Table t_Table structure containing data, header, ...
*           modified Table t_Table structure containing data, header, ...
*           number of read elements (-1: error, 0:header only)
* The routine stores any line starting with '#', '%' and ';' into the header
* Other lines are interpreted as numerical data, and stored.
* Data block should be a rectangular matrix or vector.
* Data block may be rebined with Table_Rebin (also sort in ascending order)
*******************************************************************************/
  long Table_Read_Handle(t_Table *Table, FILE *hfile,
                         long block_number, long max_rows, char *name)
  { /* reads all/a data block from 'file' handle and returns a Table structure  */
    double *Data;
    char *Header              = NULL;
    long  malloc_size         = CHAR_BUF_LENGTH;
    long  malloc_size_h       = 4096;
    long  Rows = 0,   Columns = 0;
    long  count_in_array      = 0;
    long  count_in_header     = 0;
    long  block_Current_index = 0;
    char  flag_End_row_loop   = 0;

    if (!Table) return(-1);
    Table_Init(Table, 0, 0);
    if (name && name[0]!='\0') strncpy(Table->filename, name, 1024);

    if(!hfile) {
       fprintf(stderr, "Error: File handle is NULL (Table_Read_Handle).\n");
       return (-1);
    }
    Header = (char*)  calloc(malloc_size_h, sizeof(char));
    Data   = (double*)calloc(malloc_size,   sizeof(double));
    if ((Header == NULL) || (Data == NULL)) {
       fprintf(stderr, "Error: Could not allocate Table and Header (Table_Read_Handle).\n");
       return (-1);
    }

    int flag_In_array = 0;
    do { /* while (!flag_End_row_loop) */
      char  line[1024*CHAR_BUF_LENGTH];
      long  back_pos=0;   /* ftell start of line */

      back_pos = ftell(hfile);
      if (fgets(line, 1024*CHAR_BUF_LENGTH, hfile) != NULL) { /* analyse line */
        /* first skip blank and tabulation characters */
        int i = strspn(line, " \t");

        /* handle comments: stored in header */
        if (NULL != strchr("#%;/", line[i]))
        { /* line is a comment */
          count_in_header += strlen(line);
          if (count_in_header >= malloc_size_h) {
            /* if succeed and in array : add (and realloc if necessary) */
            malloc_size_h = count_in_header+4096;
            Header        = (char*)realloc(Header, malloc_size_h*sizeof(char));
          }
          strncat(Header, line, 4096);
          flag_In_array=0;
          /* exit line and file if passed desired block */
          if (block_number > 0 && block_number == block_Current_index) {
            flag_End_row_loop = 1;
          }

          /* Continue with next line */
          continue;
        }

        /* get the number of columns splitting line with strtok */
        char  *lexeme;
        char  flag_End_Line = 0;
        long  block_Num_Columns = 0;
        const char seps[] = " ,;\t\n\r";

        lexeme = strtok(line, seps);
        while (!flag_End_Line) {
          if ((lexeme != NULL) && (lexeme[0] != '\0')) {
            /* reading line: the token is not empty */
            double X;
            int    count=1;
            /* test if we have 'NaN','Inf' */
            if (!strncasecmp(lexeme,"NaN",3))
              X = 0;
            else if (!strncasecmp(lexeme,"Inf",3) || !strncasecmp(lexeme,"+Inf",4))
              X = FLT_MAX;
            else if (!strncasecmp(lexeme,"-Inf",4))
              X = -FLT_MAX;
            else
              count = sscanf(lexeme,"%lg",&X);
            if (count == 1) {
              /* reading line: the token is a number in the line */
              if (!flag_In_array) {
                /* reading num: not already in a block: starts a new data block */
                block_Current_index++;
                flag_In_array    = 1;
                block_Num_Columns= 0;
                if (block_number > 0) {
                  /* initialise a new data block */
                  Rows = 0;
                  count_in_array = 0;
                } /* else append */
              }
              /* reading num: all blocks or selected block */
              if (flag_In_array && (block_number == 0 ||
                  block_number == block_Current_index)) {
                /* starting block: already the desired number of rows ? */
                if (block_Num_Columns == 0 &&
                    max_rows > 0 && Rows >= max_rows) {
                  flag_End_Line      = 1;
                  flag_End_row_loop  = 1;
                  flag_In_array      = 0;
                  /* reposition to begining of line (ignore line) */
                  fseek(hfile, back_pos, SEEK_SET);
                } else { /* store into data array */
                  if (count_in_array >= malloc_size) {
                    /* realloc data buffer if necessary */
                    malloc_size = count_in_array*1.5;
                    Data = (double*) realloc(Data, malloc_size*sizeof(double));
                    if (Data == NULL) {
                      fprintf(stderr, "Error: Can not re-allocate memory %li (Table_Read_Handle).\n",
                              malloc_size*sizeof(double));
                      return (-1);
                    }
                  }
                  if (0 == block_Num_Columns) Rows++;
                  Data[count_in_array] = X;
                  count_in_array++;
                  block_Num_Columns++;
                }
              } /* reading num: end if flag_In_array */
            } /* end reading num: end if sscanf lexeme -> numerical */
            else {
              /* reading line: the token is not numerical in that line. end block */
              if (block_Current_index == block_number) {
                flag_End_Line = 1;
                flag_End_row_loop = 1;
              } else {
                flag_In_array = 0;
                flag_End_Line = 1;
              }
            }
          }
          else {
            /* no more tokens in line */
            flag_End_Line = 1;
            if (block_Num_Columns > 0) Columns = block_Num_Columns;
          }

          // parse next token
          lexeme = strtok(NULL, seps);

        } /* while (!flag_End_Line) */
      } /* end: if fgets */
      else flag_End_row_loop = 1; /* else fgets : end of file */

    } while (!flag_End_row_loop); /* end while flag_End_row_loop */

    Table->block_number = block_number;
    Table->array_length = 1;

    // shrink header to actual size (plus terminating 0-byte)
    if (count_in_header) {
      Header = (char*)realloc(Header, count_in_header*sizeof(char) + 1);
    }
    Table->header = Header;

    if (count_in_array*Rows*Columns == 0)
    {
      Table->rows         = 0;
      Table->columns      = 0;
      free(Data);
      return (0);
    }
    if (Rows * Columns != count_in_array)
    {
      fprintf(stderr, "Warning: Read_Table :%s %s Data has %li values that should be %li x %li\n",
        (Table->filename[0] != '\0' ? Table->filename : ""),
        (!block_number ? " catenated" : ""),
        count_in_array, Rows, Columns);
      Columns = count_in_array; Rows = 1;
    }
    Data     = (double*)realloc(Data, count_in_array*sizeof(double));
    Table->data         = Data;
    Table->rows         = Rows;
    Table->columns      = Columns;

    return (count_in_array);

  } /* end Table_Read_Handle */

/*******************************************************************************
* long Table_Rebin(t_Table *Table)
*   ACTION: rebin a single Table, sorting 1st column in ascending order
*   input   Table: single table containing data.
*                  The data block is reallocated in this process
*   return  updated Table with increasing, evenly spaced first column (index 0)
*           number of data elements (-1: error, 0:empty data)
*******************************************************************************/
  long Table_Rebin(t_Table *Table)
  {
    double new_step=0;
    long   i;
    /* performs linear interpolation on X axis (0-th column) */

    if (!Table) return(-1);
    if (!Table->data 
    || Table->rows*Table->columns == 0 || !Table->step_x)
      return(0);
    Table_Stat(Table); /* recompute statitstics and minimal step */
    new_step = Table->step_x; /* minimal step in 1st column */

    if (!(Table->constantstep)) /* not already evenly spaced */
    {
      long Length_Table;
      double *New_Table;

      Length_Table = ceil(fabs(Table->max_x - Table->min_x)/new_step)+1;
      /*return early if the rebinned table will become too large*/
      if (Length_Table > mcread_table_rebin_maxsize){
        fprintf(stderr,"WARNING: (Table_Rebin): Rebinning table from %s would exceed 1M rows. Skipping.\n", Table->filename); 
        return(Table->rows*Table->columns);
      }
      New_Table    = (double*)malloc(Length_Table*Table->columns*sizeof(double));

      for (i=0; i < Length_Table; i++)
      {
        long   j;
        double X;
        X = Table->min_x + i*new_step;
        New_Table[i*Table->columns] = X;
        for (j=1; j < Table->columns; j++)
          New_Table[i*Table->columns+j]
                = Table_Value(*Table, X, j);
      } /* end for i */

      Table->rows = Length_Table;
      Table->step_x = new_step;
      Table->max_x = Table->min_x + (Length_Table-1)*new_step; 
      /*max might not be the same anymore
       * Use Length_Table -1 since the first and laset rows are the limits of the defined interval.*/
      free(Table->data);
      Table->data = New_Table;
      Table->constantstep=1;
    } /* end else (!constantstep) */
    return (Table->rows*Table->columns);
  } /* end Table_Rebin */

/*******************************************************************************
* double Table_Index(t_Table Table, long i, long j)
*   ACTION: read an element [i,j] of a single Table
*   input   Table: table containing data
*           i : index of row      (0:Rows-1)
*           j : index of column   (0:Columns-1)
*   return  Value = data[i][j]
* Returns Value from the i-th row, j-th column of Table
* Tests are performed on indexes i,j to avoid errors
*******************************************************************************/

#ifndef MIN
#define MIN(a, b)  (((a) < (b)) ? (a) : (b))
#endif
#ifndef MAX
#define MAX(a, b)  (((a) > (b)) ? (a) : (b))
#endif

double Table_Index(t_Table Table, long i, long j)
{
  long AbsIndex;

  if (Table.rows == 1 || Table.columns == 1) {
    /* vector */
    j = MIN(MAX(0, i+j), Table.columns*Table.rows - 1);
    i = 0;
  } else {
    /* matrix */
    i = MIN(MAX(0, i), Table.rows - 1);
    j = MIN(MAX(0, j), Table.columns - 1);
  }

  /* handle vectors specifically */
  AbsIndex = i*(Table.columns)+j;

  if (Table.data != NULL)
    return (Table.data[AbsIndex]);
  else
    return 0;
} /* end Table_Index */

/*******************************************************************************
* void Table_SetElement(t_Table *Table, long i, long j, double value)
*   ACTION: set an element [i,j] of a single Table
*   input   Table: table containing data
*           i : index of row      (0:Rows-1)
*           j : index of column   (0:Columns-1)
*           value = data[i][j]
* Returns 0 in case of error
* Tests are performed on indexes i,j to avoid errors
*******************************************************************************/
int Table_SetElement(t_Table *Table, long i, long j,
                     double value)
{
  long AbsIndex;

  if (Table->rows == 1 || Table->columns == 1) {
    /* vector */
    j = MIN(MAX(0, i+j), Table->columns*Table->rows - 1); i=0;
  } else {
    /* matrix */
    i = MIN(MAX(0, i), Table->rows - 1);
    j = MIN(MAX(0, j), Table->columns - 1);
  }

  AbsIndex = i*(Table->columns)+j;
  if (Table->data != NULL) {
    Table->data[AbsIndex] = value;
    return 1;
  }

  return 0;
} /* end Table_SetElement */

/*******************************************************************************
* double Table_Value(t_Table Table, double X, long j)
*   ACTION: read column [j] of a single Table at row which 1st column is X
*   input   Table: table containing data.
*           X : data value in the first column (index 0)
*           j : index of column from which is extracted the Value (0:Columns-1)
*   return  Value = data[index for X][j] with linear interpolation
* Returns Value from the j-th column of Table corresponding to the
* X value for the 1st column (index 0)
* Tests are performed (within Table_Index) on indexes i,j to avoid errors
* NOTE: data should rather be monotonic, and evenly sampled.
*******************************************************************************/
double Table_Value(t_Table Table, double X, long j)
{
  long   Index = -1;
  double X1=0, Y1=0, X2=0, Y2=0;
  double ret=0;

  if (X > Table.max_x) return Table_Index(Table,Table.rows-1  ,j);
  if (X < Table.min_x) return Table_Index(Table,0  ,j);

  // Use constant-time lookup when possible
  if(Table.constantstep) {
    Index = (long)floor(
              (X - Table.min_x) / (Table.max_x - Table.min_x) * (Table.rows-1));
    X1 = Table_Index(Table,Index  ,0);
    X2 = Table_Index(Table,Index+1,0);
  }
  // Use binary search on large, monotonic tables
  else if(Table.monotonic && Table.rows > 100) {
    long left = Table.min_x;
    long right = Table.max_x;

    while (!((X1 <= X) && (X < X2)) && (right - left > 1)) {
      Index = (left + right) / 2;

      X1 = Table_Index(Table, Index-1, 0);
      X2 = Table_Index(Table, Index,   0);

      if (X < X1) {
        right = Index;
      } else {
        left  = Index;
      }
    }
  }

  // Fall back to linear search, if no-one else has set X1, X2 correctly
  if (!((X1 <= X) && (X < X2))) {
    /* look for index surrounding X in the table -> Index */
    for (Index=1; Index <= Table.rows-1; Index++) {
        X1 = Table_Index(Table, Index-1,0);
        X2 = Table_Index(Table, Index  ,0);
        if ((X1 <= X) && (X < X2)) break;
      } /* end for Index */
  }

  Y1 = Table_Index(Table,Index-1,j);
  Y2 = Table_Index(Table,Index  ,j);

  if (!strcmp(Table.method,"linear")) {
    ret = Table_Interp1d(X, X1,Y1, X2,Y2);
  }
  else if (!strcmp(Table.method,"nearest")) {
    ret = Table_Interp1d_nearest(X, X1,Y1, X2,Y2);
  }

  return ret;
} /* end Table_Value */

/*******************************************************************************
* double Table_Value2d(t_Table Table, double X, double Y)
*   ACTION: read element [X,Y] of a matrix Table
*   input   Table: table containing data.
*           X : row index, may be non integer
*           Y : column index, may be non integer
*   return  Value = data[index X][index Y] with bi-linear interpolation
* Returns Value for the indices [X,Y]
* Tests are performed (within Table_Index) on indexes i,j to avoid errors
* NOTE: data should rather be monotonic, and evenly sampled.
*******************************************************************************/
  double Table_Value2d(t_Table Table, double X, double Y)
  {
    long   x1,x2,y1,y2;
    double z11,z12,z21,z22;
    double ret=0;

    x1 = (long)floor(X);
    y1 = (long)floor(Y);

    if (x1 > Table.rows-1 || x1 < 0) {
      x2 = x1;
    } else {
      x2 = x1 + 1;
    }

    if (y1 > Table.columns-1 || y1 < 0) {
      y2 = y1;
    } else {
      y2 = y1 + 1;
    }

    z11 = Table_Index(Table, x1, y1);

    if (y2 != y1) z12=Table_Index(Table, x1, y2); else z12 = z11;
    if (x2 != x1) z21=Table_Index(Table, x2, y1); else z21 = z11;
    if (y2 != y1) z22=Table_Index(Table, x2, y2); else z22 = z21;

    if (!strcmp(Table.method,"linear"))
      ret = Table_Interp2d(X,Y, x1,y1,x2,y2, z11,z12,z21,z22);
    else {
      if (fabs(X-x1) < fabs(X-x2)) {
        if (fabs(Y-y1) < fabs(Y-y2)) ret = z11; else ret = z12;
      } else {
        if (fabs(Y-y1) < fabs(Y-y2)) ret = z21; else ret = z22;
      }
    }
    return ret;
  } /* end Table_Value2d */


/*******************************************************************************
* void Table_Free(t_Table *Table)
*   ACTION: free a single Table. First Call Table_File_list_gc. If this returns
*   non-zero it means there are more refernces to the table, and so the table
*   should not bee freed.
*   return: empty Table
*******************************************************************************/
  void Table_Free(t_Table *Table)
  {
    if( !Table_File_List_gc(Table) ){
       return;
    } 
    if (!Table) return;
    if (Table->data   != NULL) free(Table->data);
    if (Table->header != NULL) free(Table->header);
    Table->data   = NULL;
    Table->header = NULL;
  } /* end Table_Free */

/******************************************************************************
* void Table_Info(t_Table Table)
*    ACTION: print informations about a single Table
*******************************************************************************/
  long Table_Info(t_Table Table)
  {
    char buffer[256];
    long ret=0;

    if (!Table.block_number) strcpy(buffer, "catenated");
    else sprintf(buffer, "block %li", Table.block_number);
    printf("Table from file '%s' (%s)",
      Table.filename[0] != '\0' ? Table.filename : "", buffer);
    if ((Table.data != NULL) && (Table.rows*Table.columns))
    {
      printf(" is %li x %li ", Table.rows, Table.columns);
      if (Table.rows*Table.columns > 1)
        printf("(x=%g:%g)", Table.min_x, Table.max_x);
      else printf("(x=%g) ", Table.min_x);
      ret = Table.rows*Table.columns;
      if (Table.monotonic)    printf(", monotonic");
      if (Table.constantstep) printf(", constant step");
      printf(". interpolation: %s\n", Table.method);
    }
    else printf(" is empty.\n");

    if (Table.header && strlen(Table.header)) {
      char *header;
      int  i;
      header = malloc(80);
      if (!header) return(ret);
      for (i=0; i<80; header[i++]=0);
      strncpy(header, Table.header, 75);
      if (strlen(Table.header) > 75) {
        strcat( header, " ...");
      }
      for (i=0; i<strlen(header); i++)
        if (header[i] == '\n' || header[i] == '\r') header[i] = ';';
      printf("  '%s'\n", header);
      free(header);
    }

    return(ret);
  } /* end Table_Info */

/******************************************************************************
* long Table_Init(t_Table *Table, m, n)
*   ACTION: initialise a Table to empty m by n table
*   return: empty Table
******************************************************************************/
long Table_Init(t_Table *Table, long rows, long columns)
{
  double *data=NULL;
  long   i;

  if (!Table) return(0);

  Table->header  = NULL;
  Table->filename[0]= '\0';
  Table->filesize= 0;
  Table->min_x   = 0;
  Table->max_x   = 0;
  Table->step_x  = 0;
  Table->block_number = 0;
  Table->array_length = 0;
  Table->monotonic    = 0;
  Table->constantstep = 0;
  Table->begin   = 0;
  Table->end     = 0;
  strcpy(Table->method,"linear");

  if (rows*columns >= 1) {
    data    = (double*)malloc(rows*columns*sizeof(double));
    if (data) for (i=0; i < rows*columns; data[i++]=0);
    else {
      fprintf(stderr,"Error: allocating %ld double elements."
                     "Too big (Table_Init).\n", rows*columns);
      rows = columns = 0;
    }
  }
  Table->rows    = (rows >= 1 ? rows : 0);
  Table->columns = (columns >= 1 ? columns : 0);
  Table->data    = data;
  return(Table->rows*Table->columns);
} /* end Table_Init */

/******************************************************************************
* long Table_Write(t_Table Table, char *file, x1,x2, y1,y2)
*   ACTION: write a Table to disk (ascii).
*     when x1=x2=0 or y1=y2=0, the table default limits are used.
*   return: 0=all is fine, non-0: error
*******************************************************************************/
MCDETECTOR Table_Write(t_Table Table, char *file, char *xl, char *yl, 
  double x1, double x2, double y1, double y2)
{
  long    i =0;
  MCDETECTOR detector;

  if ((Table.data == NULL) && (Table.rows*Table.columns)) {
    detector.m = 0;
    return(detector); /* Table is empty - nothing to do */
  }
  if (!x1 && !x2) {
    x1 = Table.min_x;
    x2 = Table.max_x;
  }
  if (!y1 && !y2) {
    y1 = 1;
    y2 = Table.columns;
  }

  /* transfer content of the Table into a 2D detector */
  Coords coords = { 0, 0, 0};

  if (Table.rows == 1 || Table.columns == 1) {
    detector = mcdetector_out_1D(Table.filename,
                      xl ? xl : "", yl ? yl : "",
                      "x", x1, x2,
                      Table.rows * Table.columns,
                      NULL, Table.data, NULL,
                      file, file, coords);
  } else {
    detector = mcdetector_out_2D(Table.filename,
                      xl ? xl : "", yl ? yl : "",
                      x1, x2, y1, y2,
                      Table.rows, Table.columns,
                      NULL, Table.data, NULL,
                      file, file, coords);
  }
  return(detector);
}

/******************************************************************************
* void Table_Stat(t_Table *Table)
*   ACTION: computes min/max/mean step of 1st column for a single table (private)
*   return: updated Table
*******************************************************************************/
  static void Table_Stat(t_Table *Table)
  {
    long   i;
    double max_x, min_x;
    double row=1;
    char   monotonic=1;
    char   constantstep=1;
    double step=0;
    long n;

    if (!Table) return;
    if (!Table->rows || !Table->columns) return;
    if (Table->rows == 1) row=0; // single row
    max_x = -FLT_MAX;
    min_x =  FLT_MAX;
    n     = (row ? Table->rows : Table->columns);
    /* get min and max of first column/vector */
    for (i=0; i < n; i++)
    {
      double X;
      X = (row ? Table_Index(*Table,i  ,0)
                               : Table_Index(*Table,0, i));
      if (X < min_x) min_x = X;
      if (X > max_x) max_x = X;
    } /* for */
    
    /* test for monotonicity and constant step if the table is an XY or single vector */
    if (n > 1) {
      /* mean step */
      step = (max_x - min_x)/(n-1);
      /* now test if table is monotonic on first column, and get minimal step size */
      for (i=0; i < n-1; i++) {
        double X, diff;;
        X    = (row ? Table_Index(*Table,i  ,0)
                    : Table_Index(*Table,0,  i));
        diff = (row ? Table_Index(*Table,i+1,0)
                    : Table_Index(*Table,0,  i+1)) - X;
        if (diff && fabs(diff) < fabs(step)) step = diff;
        /* change sign ? */
        if ((max_x - min_x)*diff < 0 && monotonic)
          monotonic = 0;
      } /* end for */
      
      /* now test if steps are constant within READ_TABLE_STEPTOL */
      if(!step){
        /*means there's a disconitnuity -> not constantstep*/
        constantstep=0;
      }else if (monotonic) {
        for (i=0; i < n-1; i++) {
          double X, diff;
          X    = (row ? Table_Index(*Table,i  ,0)
              : Table_Index(*Table,0,  i));
          diff = (row ? Table_Index(*Table,i+1,0)
              : Table_Index(*Table,0,  i+1)) - X;
          if ( fabs(step)*(1+READ_TABLE_STEPTOL) < fabs(diff) ||
                fabs(diff) < fabs(step)*(1-READ_TABLE_STEPTOL) )
          { constantstep = 0; break; }
        }
      }

    }
    Table->step_x= step;
    Table->max_x = max_x;
    Table->min_x = min_x;
    Table->monotonic = monotonic;
    Table->constantstep = constantstep;
  } /* end Table_Stat */

/******************************************************************************
* t_Table *Table_Read_Array(char *File, long *blocks)
*   ACTION: read as many data blocks as available, iteratively from file
*   return: initialized t_Table array, last element is an empty Table.
*           the number of extracted blocks in non NULL pointer *blocks
*******************************************************************************/
  t_Table *Table_Read_Array(char *File, long *blocks)
  {
    t_Table *Table_Array=NULL;
    long offset=0;
    long block_number=0;
    long allocated=256;
    long nelements=1;

    /* first allocate an initial empty t_Table array */
    Table_Array = (t_Table *)malloc(allocated*sizeof(t_Table));
    if (!Table_Array) {
      fprintf(stderr, "Error: Can not allocate memory %li (Table_Read_Array).\n",
         allocated*sizeof(t_Table));
      *blocks = 0;
      return (NULL);
    }

    while (nelements > 0)
    {
      t_Table Table;

      /* if ok, set t_Table block number else exit loop */
      block_number++;
      Table.block_number = block_number;
      
      /* access file at offset and get following block. Block number is from the set offset
       * hence the hardcoded 1 - i.e. the next block counted from offset.*/
      nelements = Table_Read_Offset(&Table, File, 1, &offset,0);
      /*if the block is empty - don't store it*/
      if (nelements>0){
          /* if t_Table array is not long enough, expand and realocate */
          if (block_number >= allocated-1) {
              allocated += 256;
              Table_Array = (t_Table *)realloc(Table_Array,
                      allocated*sizeof(t_Table));
              if (!Table_Array) {
                  fprintf(stderr, "Error: Can not re-allocate memory %li (Table_Read_Array).\n",
                          allocated*sizeof(t_Table));
                  *blocks = 0;
                  return (NULL);
              }
          }
          /* store it into t_Table array */
          //snprintf(Table.filename, 1024, "%s#%li", File, block_number-1);
          Table_Array[block_number-1] = Table;
      }
      /* continues until we find an empty block */
    }
    /* send back number of extracted blocks */
    if (blocks) *blocks = block_number-1;

    /* now store total number of elements in Table array */
    for (offset=0; offset < block_number;
      Table_Array[offset++].array_length = block_number-1);

    return(Table_Array);
  } /* end Table_Read_Array */
/*******************************************************************************
* void Table_Free_Array(t_Table *Table)
*   ACTION: free a Table array
*******************************************************************************/
  void Table_Free_Array(t_Table *Table)
  {
    long index;
    if (!Table) return;
    for (index=0;index < Table[0].array_length; index++){
            Table_Free(&Table[index]);
    }
    free(Table);
  } /* end Table_Free_Array */

/******************************************************************************
* long Table_Info_Array(t_Table *Table)
*    ACTION: print informations about a Table array
*    return: number of elements in the Table array
*******************************************************************************/
  long Table_Info_Array(t_Table *Table)
  {
    long index=0;

    if (!Table) return(-1);
    while (index < Table[index].array_length
       && (Table[index].data || Table[index].header)
       && (Table[index].rows*Table[index].columns) ) {
      Table_Info(Table[index]);
      index++;
    }
    printf("This Table array contains %li elements\n", index);
    return(index);
  } /* end Table_Info_Array */

/******************************************************************************
* char **Table_ParseHeader(char *header, symbol1, symbol2, ..., NULL)
*    ACTION: search for char* symbols in header and return their value or NULL
*            the search is not case sensitive.
*            Last argument MUST be NULL
*    return: array of char* with line following each symbol, or NULL if not found
*******************************************************************************/
#ifndef MyNL_ARGMAX
#define MyNL_ARGMAX 50
#endif

char **Table_ParseHeader_backend(char *header, ...){
  va_list ap;
  char exit_flag=0;
  int counter   =0;
  char **ret    =NULL;
  if (!header || header[0]=='\0') return(NULL);

  ret = (char**)calloc(MyNL_ARGMAX, sizeof(char*));
  if (!ret) {
    printf("Table_ParseHeader: Cannot allocate %i values array for Parser (Table_ParseHeader).\n",
      MyNL_ARGMAX);
    return(NULL);
  }
  for (counter=0; counter < MyNL_ARGMAX; ret[counter++] = NULL);
  counter=0;

  va_start(ap, header);
  while(!exit_flag && counter < MyNL_ARGMAX-1)
  {
    char *arg_char=NULL;
    char *pos     =NULL;
    /* get variable argument value as a char */
    arg_char = va_arg(ap, char *);
    if (!arg_char || arg_char[0]=='\0'){
      exit_flag = 1; break;
    }
    /* search for the symbol in the header */
    pos = (char*)strcasestr(header, arg_char);
    if (pos) {
      char *eol_pos;
      eol_pos = strchr(pos+strlen(arg_char), '\n');
      if (!eol_pos)
        eol_pos = strchr(pos+strlen(arg_char), '\r');
      if (!eol_pos)
        eol_pos = pos+strlen(pos)-1;
      ret[counter] = (char*)malloc(eol_pos - pos);
      if (!ret[counter]) {
        printf("Table_ParseHeader: Cannot allocate value[%i] array for Parser searching for %s (Table_ParseHeader).\n",
          counter, arg_char);
        exit_flag = 1; break;
      }
      strncpy(ret[counter], pos+strlen(arg_char), eol_pos - pos - strlen(arg_char));
      ret[counter][eol_pos - pos - strlen(arg_char)]='\0';
    }
    counter++;
  }
  va_end(ap);
  return(ret);
} /* Table_ParseHeader */

/******************************************************************************
* double Table_Interp1d(x, x1, y1, x2, y2)
*    ACTION: interpolates linearly at x between y1=f(x1) and y2=f(x2)
*    return: y=f(x) value
*******************************************************************************/
double Table_Interp1d(double x,
  double x1, double y1,
  double x2, double y2)
{
  double slope;
  if (x2 == x1) return (y1+y2)/2;
  if (y1 == y2) return  y1;
  slope = (y2 - y1)/(x2 - x1);
  return y1+slope*(x - x1);
} /* Table_Interp1d */

/******************************************************************************
* double Table_Interp1d_nearest(x, x1, y1, x2, y2)
*    ACTION: table lookup with nearest method at x between y1=f(x1) and y2=f(x2)
*    return: y=f(x) value
*******************************************************************************/
double Table_Interp1d_nearest(double x,
  double x1, double y1,
  double x2, double y2)
{
  if (fabs(x-x1) < fabs(x-x2)) return (y1);
  else return(y2);
} /* Table_Interp1d_nearest */

/******************************************************************************
* double Table_Interp2d(x,y, x1,y1, x2,y2, z11,z12,z21,z22)
*    ACTION: interpolates bi-linearly at (x,y) between z1=f(x1,y1) and z2=f(x2,y2)
*    return: z=f(x,y) value
*    x,y |   x1   x2
*    ----------------
*     y1 |   z11  z21
*     y2 |   z12  z22
*******************************************************************************/
double Table_Interp2d(double x, double y,
  double x1, double y1,
  double x2, double y2,
  double z11, double z12, double z21, double z22)
{
  double ratio_x, ratio_y;
  if (x2 == x1) return Table_Interp1d(y, y1,z11, y2,z12);
  if (y1 == y2) return Table_Interp1d(x, x1,z11, x2,z21);

  ratio_y = (y - y1)/(y2 - y1);
  ratio_x = (x - x1)/(x2 - x1);
  return (1-ratio_x)*(1-ratio_y)*z11 + ratio_x*(1-ratio_y)*z21
    + ratio_x*ratio_y*z22         + (1-ratio_x)*ratio_y*z12;
} /* Table_Interp2d */

/* end of read_table-lib.c */


#ifndef REF_LIB_H
#define REF_LIB_H "$Revision$"

void StdReflecFunc(double, double*, double*);
void TableReflecFunc(double, t_Table*, double*);

#endif

/* end of ref-lib.h */
/****************************************************************************
*
* McStas, neutron ray-tracing package
*         Copyright 1997-2006, All rights reserved
*         Risoe National Laboratory, Roskilde, Denmark
*         Institut Laue Langevin, Grenoble, France
*
* Library: share/ref-lib.c
*
* %Identification
* Written by: Peter Christiansen
* Date: August, 2006
* Origin: RISOE
* Release: McStas 1.10
* Version: $Revision$
*
* Commonly used reflection functions are declared in this file which
* are used by some guide and mirror components.
*
* Variable names have prefix 'mc_ref_' for 'McStas Reflection' 
* to avoid conflicts
*
* Usage: within SHARE
* %include "ref-lib"
*
****************************************************************************/

#ifndef REF_LIB_H
#include "ref-lib.h"
#endif

#ifndef READ_TABLE_LIB_H
#include "read_table-lib.h"
#include "read_table-lib.c"
#endif

/****************************************************************************
* void StdReflecFunc(double q, double *par, double *r)
* 
* The McStas standard analytic parametrization of the reflectivity.
* The parameters are:
* R0:      [1]    Low-angle reflectivity
* Qc:      [AA-1] Critical scattering vector
* alpha:   [AA]   Slope of reflectivity
* m:       [1]    m-value of material. Zero means completely absorbing.
* W:       [AA-1] Width of supermirror cut-off
*****************************************************************************/
void StdReflecFunc(double mc_pol_q, double *mc_pol_par, double *mc_pol_r) {
    double R0    = mc_pol_par[0];
    double Qc    = mc_pol_par[1];
    double alpha = mc_pol_par[2];
    double m     = mc_pol_par[3];
    double W     = mc_pol_par[4];
    double beta  = 0;
    mc_pol_q     = fabs(mc_pol_q);
    double arg;
        
    /* Simpler parametrization from Henrik Jacobsen uses these values that depend on m only.
       double m_value=m*0.9853+0.1978;
       double W=-0.0002*m_value+0.0022;
       double alpha=0.2304*m_value+5.0944;
       double beta=-7.6251*m_value+68.1137; 
       If W and alpha are set to 0, use Henrik's approach for estimating these parameters
       and apply the formulation:
       arg = R0*0.5*(1-tanh(arg))*(1-alpha*(q-Qc)+beta*(q-Qc)*(q-Qc));
    */  
    if (W==0 && alpha==0) {
      m=m*0.9853+0.1978;
      W=-0.0002*m+0.0022;
      alpha=0.2304*m+5.0944;
      beta=-7.6251*m+68.1137;
      if (m<=3) {
	alpha=m;
	beta=0;
      }
    }
    
    arg = W > 0 ? (mc_pol_q - m*Qc)/W : 11;

    if (arg > 10 || m <= 0 || Qc <=0 || R0 <= 0) {
      *mc_pol_r = 0;
      return;
    }
    
    if (m < 1) { Qc *= m; m=1; }
    
    if(mc_pol_q <= Qc) {      
      *mc_pol_r = R0;
      return;
    }
    
    
    *mc_pol_r = R0*0.5*(1 - tanh(arg))*(1 - alpha*(mc_pol_q - Qc) + beta*(mc_pol_q - Qc)*(mc_pol_q - Qc));
    
    return;
  }

/****************************************************************************
* void TableReflecFunc(double q, t_Table *par, double *r) {
* 
* Looks up the reflectivity in a table using the routines in read_table-lib.
*****************************************************************************/
void TableReflecFunc(double mc_pol_q, t_Table *mc_pol_par, double *mc_pol_r) {
    
  *mc_pol_r = Table_Value(*mc_pol_par, mc_pol_q, 1);
  if(*mc_pol_r>1)
    *mc_pol_r = 1;
  return;
}

/* end of ref-lib.c */

#ifndef Gravity_guide_Version
#define Gravity_guide_Version "$Revision$"

#ifndef PROP_GRAV_DT
#error McStas : You need PROP_GRAV_DT (McStas >= 1.4.3) to run this component
#endif

/*
* G:       (m/s^2) Gravitation acceleration along y axis [-9.81]
* Gx:      (m/s^2) Gravitation acceleration along x axis [0]
* Gy:      (m/s^2) Gravitation acceleration along y axis [-9.81]
* Gz:      (m/s^2) Gravitation acceleration along z axis [0]
* mh:      (1)    m-value of material for left/right vert. mirrors
* mv:      (1)    m-value of material for top/bottom horz. mirrors
* mx:      (1)    m-value of material for left/right vert. mirrors
* my:      (1)    m-value of material for top/bottom horz. mirrors
*/

  typedef struct Gravity_guide_Vars
  {
    double gx;
    double gy;
    double gz;
    double nx[6], ny[6], nz[6];
    double wx[6], wy[6], wz[6];
    double A[6], norm_n2[6], norm_n[6];
    long   N_reflection[7];
    double w1c, h1c;
    double w2c, h2c;
    double M[5];
    double Alpha[5];
    double nzC[5], norm_n2xy[5], Axy[5];
    double wav_lr, wav_tb, wav_z;
    double chamfer_z, chamfer_lr, chamfer_tb;
    char   compcurname[256];
    double fc_freq, fc_phase;
    double warnings;
  } Gravity_guide_Vars_type;

  void Gravity_guide_Init(Gravity_guide_Vars_type *aVars,
    MCNUM a_w1, MCNUM a_h1, MCNUM a_w2, MCNUM a_h2, MCNUM a_l, MCNUM a_R0,
    MCNUM a_Qc, MCNUM a_alpha, MCNUM a_m, MCNUM a_W, MCNUM a_nslit, MCNUM a_d,
    MCNUM a_Gx, MCNUM a_Gy, MCNUM a_Gz,
    MCNUM a_mleft, MCNUM a_mright, MCNUM a_mtop, MCNUM a_mbottom, MCNUM a_nhslit,
    MCNUM a_wavy_lr, MCNUM a_wavy_tb, MCNUM a_wavy_z, MCNUM a_wavy,
    MCNUM a_chamfers_z, MCNUM a_chamfers_lr, MCNUM a_chamfers_tb, MCNUM a_chamfers,
    MCNUM a_nu, MCNUM a_phase, MCNUM a_aleft, MCNUM a_aright, MCNUM a_atop, MCNUM a_abottom)
  {
    int i;

    for (i=0; i<7; aVars->N_reflection[i++] = 0);
    for (i=0; i<5; aVars->M[i++] = 0);
    for (i=0; i<5; aVars->Alpha[i++] = 0);

    aVars->gx = a_Gx; /* The gravitation vector in the current component axis system */
    aVars->gy = a_Gy;
    aVars->gz = a_Gz;
    aVars->warnings=0;

    if (a_nslit <= 0 || a_nhslit <= 0) { fprintf(stderr,"%s: Fatal: no channel in this guide (nhslit or nslit=0).\n", aVars->compcurname); exit(-1); }
    if (a_d < 0) { fprintf(stderr,"%s: Fatal: subdividing walls have negative thickness in this guide (d<0).\n", aVars->compcurname); exit(-1); }
    aVars->w1c = (a_w1 - (a_nslit-1) *a_d)/(double)a_nslit;
    aVars->w2c = (a_w2 - (a_nslit-1) *a_d)/(double)a_nslit;
    aVars->h1c = (a_h1 - (a_nhslit-1)*a_d)/(double)a_nhslit;
    aVars->h2c = (a_h2 - (a_nhslit-1)*a_d)/(double)a_nhslit;

    for (i=0; i <= 4;   aVars->M[i++]=a_m);
    for (i=0; i <= 4;   aVars->Alpha[i++]=a_alpha);
    if (a_mleft   >= 0) aVars->M[1] =a_mleft  ;
    if (a_mright  >= 0) aVars->M[2] =a_mright ;
    if (a_mtop    >= 0) aVars->M[3] =a_mtop   ;
    if (a_mbottom >= 0) aVars->M[4] =a_mbottom;
    if (a_aleft   >= 0) aVars->Alpha[1] =a_aleft  ;
    if (a_aright  >= 0) aVars->Alpha[2] =a_aright ;
    if (a_atop    >= 0) aVars->Alpha[3] =a_atop   ;
    if (a_abottom >= 0) aVars->Alpha[4] =a_abottom;

    /* n: normal vectors to surfaces */
    aVars->nx[1] =  a_l; aVars->ny[1] =  0;   aVars->nz[1] =  0.5*(aVars->w2c-aVars->w1c);  /* 1:+X left       */
    aVars->nx[2] = -a_l; aVars->ny[2] =  0;   aVars->nz[2] = -aVars->nz[1];             /* 2:-X right      */
    aVars->nx[3] =  0;   aVars->ny[3] =  a_l; aVars->nz[3] =  0.5*(aVars->h2c-aVars->h1c);  /* 3:+Y top        */
    aVars->nx[4] =  0;   aVars->ny[4] = -a_l; aVars->nz[4] = -aVars->nz[3];             /* 4:-Y bottom     */
    aVars->nx[5] =  0;   aVars->ny[5] =  0;   aVars->nz[5] =  a_l;                      /* 5:+Z exit       */
    aVars->nx[0] =  0;   aVars->ny[0] =  0;   aVars->nz[0] = -a_l;                      /* 0:Z0 input      */
    /* w: a point on these surfaces */
    aVars->wx[1] = +(aVars->w1c)/2; aVars->wy[1] =  0;              aVars->wz[1] = 0;   /* 1:+X left       */
    aVars->wx[2] = -(aVars->w1c)/2; aVars->wy[2] =  0;              aVars->wz[2] = 0;   /* 2:-X right      */
    aVars->wx[3] =  0;              aVars->wy[3] = +(aVars->h1c)/2; aVars->wz[3] = 0;   /* 3:+Y top        */
    aVars->wx[4] =  0;              aVars->wy[4] = -(aVars->h1c)/2; aVars->wz[4] = 0;   /* 4:-Y bottom     */
    aVars->wx[5] =  0;              aVars->wy[5] =  0;              aVars->wz[5] = a_l; /* 5:+Z exit       */
    aVars->wx[0] =  0;              aVars->wy[0] =  0;              aVars->wz[0] = 0;   /* 0:Z0 input      */

    for (i=0; i <= 5; i++)
    {
      aVars->A[i] = scalar_prod(aVars->nx[i], aVars->ny[i], aVars->nz[i], aVars->gx, aVars->gy, aVars->gz)/2;
      aVars->norm_n2[i] = aVars->nx[i]*aVars->nx[i] + aVars->ny[i]*aVars->ny[i] + aVars->nz[i]*aVars->nz[i];
      if (aVars->norm_n2[i] <= 0)
        { fprintf(stderr,"%s: Fatal: normal vector norm %i is null/negative ! check guide dimensions.\n", aVars->compcurname, i); exit(-1); } /* should never occur */
      else
        aVars->norm_n[i] = sqrt(aVars->norm_n2[i]);
    }
    /* partial computations for l/r/t/b sides, to save computing time */
    for (i=1; i <= 4; i++)
    { /* stores nz that changes in case non box element (focus/defocus) */
      aVars->nzC[i]      =  aVars->nz[i]; /* partial xy terms */
      aVars->norm_n2xy[i]=  aVars->nx[i]*aVars->nx[i] + aVars->ny[i]*aVars->ny[i];
      aVars->Axy[i]      = (aVars->nx[i]*aVars->gx    + aVars->ny[i]*aVars->gy)/2;
    }
    /* handle waviness init */
    if (a_wavy && (!a_wavy_tb && !a_wavy_lr && !a_wavy_z))
    { aVars->wav_tb=aVars->wav_lr=aVars->wav_z=a_wavy; }
    else
    { aVars->wav_tb=a_wavy_tb; aVars->wav_lr=a_wavy_lr; aVars->wav_z=a_wavy_z; }
    aVars->wav_tb *= DEG2RAD/(sqrt(8*log(2)));   /* Convert from deg FWHM to rad Gaussian sigma */
    aVars->wav_lr *= DEG2RAD/(sqrt(8*log(2)));
    aVars->wav_z  *= DEG2RAD/(sqrt(8*log(2)));
    /* handle chamfers init */
    if (a_chamfers && (!a_chamfers_z && !a_chamfers_lr && !a_chamfers_tb))
    { aVars->chamfer_z=aVars->chamfer_lr=aVars->chamfer_tb=a_chamfers; }
    else
    {
      aVars->chamfer_z=a_chamfers_z;
      aVars->chamfer_lr=a_chamfers_lr;
      aVars->chamfer_tb=a_chamfers_tb;
    }

    aVars->fc_freq  = a_nu;
    aVars->fc_phase = a_phase;
  }

  int Gravity_guide_Trace(double *dt,
        Gravity_guide_Vars_type *aVars,
        double cx, double cy, double cz,
        double cvx, double cvy, double cvz,
        double cxnum, double cxk, double cynum, double cyk,
        double *cnx, double *cny,double *cnz)
  {
    double B, C;
    int    ret=0;
    int    side=0;
    double n1;
    double dt0, dt_min=0;
    int    i;
    double loc_num, loc_nslit;
    int    i_slope=3;

    /* look if there is a previous intersection with guide sides */
    /* A = 0.5 n.g; B = n.v; C = n.(r-W); */
    /* 5=+Z side: n=(0, 0, -l) ; W = (0, 0, l) (at z=l, guide exit)*/
    B = aVars->nz[5]*cvz; C = aVars->nz[5]*(cz - aVars->wz[5]);
    ret = solve_2nd_order(&dt0, NULL, aVars->A[5], B, C);
    if (ret && dt0>1e-10) { dt_min = dt0; side=5; }

    loc_num = cynum; loc_nslit = cyk;
    for (i=4; i>0; i--)
    {
      if (i == 2) { i_slope=1; loc_num = cxnum; loc_nslit = cxk; }

      if (aVars->nzC[i_slope] != 0) {
        n1 = loc_nslit - 2*(loc_num);  /* slope of l/r/u/d sides depends on the channel ! */
        loc_num++; /* use partial computations to alter nz and A */
        aVars->nz[i]= aVars->nzC[i]*n1;
        aVars->A[i] = aVars->Axy[i] + aVars->nz[i]*aVars->gz/2;
      }
      if (i < 3)
      {      B = aVars->nx[i]*cvx + aVars->nz[i]*cvz; C = aVars->nx[i]*(cx-aVars->wx[i]) + aVars->nz[i]*cz; }
      else { B = aVars->ny[i]*cvy + aVars->nz[i]*cvz; C = aVars->ny[i]*(cy-aVars->wy[i]) + aVars->nz[i]*cz; }
      ret = solve_2nd_order(&dt0, NULL, aVars->A[i], B, C);
      if (ret && dt0>1e-10 && (dt0<dt_min || !dt_min))
      { dt_min = dt0; side=i;
        if (aVars->nzC[i] != 0)
        { aVars->norm_n2[i] = aVars->norm_n2xy[i] + aVars->nz[i]*aVars->nz[i];
          aVars->norm_n[i]  = sqrt(aVars->norm_n2[i]); }
      }
     }

    *dt = dt_min;
    /* handles waviness: rotate n vector */
    if (side > 0 && side < 5 && (aVars->wav_z || aVars->wav_lr || aVars->wav_tb))
    {
      double nt_x, nt_y, nt_z;  /* transverse vector */
      double nn_x, nn_y, nn_z;  /* normal vector (tmp) */
      double phi;
      /* normal vector n_z = [ 0,0,1], n_t = n x n_z; */
      vec_prod(nt_x,nt_y,nt_z, aVars->nx[side],aVars->ny[side],aVars->nz[side], 0,0,1);
      /* rotate n with angle wavy_z around n_t -> nn */
      if (aVars->wav_z) {
        phi = aVars->wav_z;
        rotate(nn_x,nn_y,nn_z, aVars->nx[side],aVars->ny[side],aVars->nz[side], aVars->wav_z*randnorm(), nt_x,nt_y,nt_z);
      } else { nn_x=aVars->nx[side]; nn_y=aVars->ny[side]; nn_z=aVars->nz[side]; }
      /* rotate n with angle wavy_{x|y} around n_z -> nt */
      phi = (side <=2) ? aVars->wav_lr : aVars->wav_tb;
      if (phi) {
        rotate(nt_x,nt_y,nt_z, nn_x,nn_y,nn_z, phi*randnorm(), 0,0,1);
      } else { nt_x=nn_x; nt_y=nn_y; nt_z=nn_z; }
      *cnx=nt_x; *cny=nt_y; *cnz=nt_z;
    } else
    { *cnx=aVars->nx[side]; *cny=aVars->ny[side]; *cnz=aVars->nz[side]; }
    return (side);
  }



#endif
#line 8993 "./ess_extraction.c"

/* Instrument parameters. */
MCNUM mcipdet_width;
MCNUM mcipdet_width_focus;
MCNUM mcipsource_width;
MCNUM mcipguide_width;
MCNUM mcipL_source;
MCNUM mcipguide_length;
MCNUM mcipdL;
MCNUM mcipL_min;
MCNUM mcipL_max;
MCNUM mcipsource_divergence;
MCNUM mcipdivergence_max;
MCNUM mcipfocal_length;
MCNUM mcipmirrors;
MCNUM mcipincoming_length;
MCNUM mcipg;
MCNUM mcipmax_div;
MCNUM mcippixels;
MCNUM mcipflux;
MCNUM mcipplaceholder;

#define mcNUMIPAR 19
int mcnumipar = 19;
struct mcinputtable_struct mcinputtable[mcNUMIPAR+1] = {
  "det_width", &mcipdet_width, instr_type_double, "0.2", 
  "det_width_focus", &mcipdet_width_focus, instr_type_double, "0.02", 
  "source_width", &mcipsource_width, instr_type_double, "0.02", 
  "guide_width", &mcipguide_width, instr_type_double, "0.1", 
  "L_source", &mcipL_source, instr_type_double, "4", 
  "guide_length", &mcipguide_length, instr_type_double, "160", 
  "dL", &mcipdL, instr_type_double, "1", 
  "L_min", &mcipL_min, instr_type_double, "3", 
  "L_max", &mcipL_max, instr_type_double, "5", 
  "source_divergence", &mcipsource_divergence, instr_type_double, "0.8", 
  "divergence_max", &mcipdivergence_max, instr_type_double, "1", 
  "focal_length", &mcipfocal_length, instr_type_double, "3", 
  "mirrors", &mcipmirrors, instr_type_double, "20", 
  "incoming_length", &mcipincoming_length, instr_type_double, "5000", 
  "g", &mcipg, instr_type_double, "-9.81", 
  "max_div", &mcipmax_div, instr_type_double, "10", 
  "pixels", &mcippixels, instr_type_double, "50", 
  "flux", &mcipflux, instr_type_double, "1", 
  "placeholder", &mcipplaceholder, instr_type_double, "1", 
  NULL, NULL, instr_type_double, ""
};

/* User declarations from instrument definition. */
#define mccompcurname  template_simple
#define mccompcurtype  INSTRUMENT
#define mccompcurindex 0
#define mcposatemplate_simple coords_set(0,0,0)
#define det_width mcipdet_width
#define det_width_focus mcipdet_width_focus
#define source_width mcipsource_width
#define guide_width mcipguide_width
#define L_source mcipL_source
#define guide_length mcipguide_length
#define dL mcipdL
#define L_min mcipL_min
#define L_max mcipL_max
#define source_divergence mcipsource_divergence
#define divergence_max mcipdivergence_max
#define focal_length mcipfocal_length
#define mirrors mcipmirrors
#define incoming_length mcipincoming_length
#define g mcipg
#define max_div mcipmax_div
#define pixels mcippixels
#define flux mcipflux
#define placeholder mcipplaceholder
#undef placeholder
#undef flux
#undef pixels
#undef max_div
#undef g
#undef incoming_length
#undef mirrors
#undef focal_length
#undef divergence_max
#undef source_divergence
#undef L_max
#undef L_min
#undef dL
#undef guide_length
#undef L_source
#undef guide_width
#undef source_width
#undef det_width_focus
#undef det_width
#undef mcposatemplate_simple
#undef mccompcurindex
#undef mccompcurtype
#undef mccompcurname

/* neutron state table at each component input (local coords) */
/* [x, y, z, vx, vy, vz, t, sx, sy, sz, p] */
MCNUM mccomp_storein[11*22];
/* Components position table (absolute and relative coords) */
Coords mccomp_posa[22];
Coords mccomp_posr[22];
/* Counter for each comp to check for inactive ones */
MCNUM  mcNCounter[22];
MCNUM  mcPCounter[22];
MCNUM  mcP2Counter[22];
#define mcNUMCOMP 21 /* number of components */
/* Counter for PROP ABSORB */
MCNUM  mcAbsorbProp[22];
/* Flag true when previous component acted on the neutron (SCATTER) */
MCNUM mcScattered=0;
/* Flag true when neutron should be restored (RESTORE) */
MCNUM mcRestore=0;
/* Declarations of component definition and setting parameters. */

/* Setting parameters for component 'origin' [1]. */
char mccorigin_profile[16384];
MCNUM mccorigin_percent;
MCNUM mccorigin_flag_save;
MCNUM mccorigin_minutes;

/* Setting parameters for component 'source_div' [3]. */
MCNUM mccsource_div_xwidth;
MCNUM mccsource_div_yheight;
MCNUM mccsource_div_focus_aw;
MCNUM mccsource_div_focus_ah;
MCNUM mccsource_div_E0;
MCNUM mccsource_div_dE;
MCNUM mccsource_div_lambda0;
MCNUM mccsource_div_dlambda;
MCNUM mccsource_div_gauss;
MCNUM mccsource_div_flux;

/* Setting parameters for component 'psd_monitor_source' [4]. */
int mccpsd_monitor_source_nx;
int mccpsd_monitor_source_ny;
char mccpsd_monitor_source_filename[16384];
MCNUM mccpsd_monitor_source_xmin;
MCNUM mccpsd_monitor_source_xmax;
MCNUM mccpsd_monitor_source_ymin;
MCNUM mccpsd_monitor_source_ymax;
MCNUM mccpsd_monitor_source_xwidth;
MCNUM mccpsd_monitor_source_yheight;
MCNUM mccpsd_monitor_source_restore_neutron;

/* Definition parameters for component 'parabolic_optic_before_guide_v' [5]. */
#define mccparabolic_optic_before_guide_v_reflect "supermirror_m3.rfl" /* declared as a string. May produce warnings at compile */
/* Setting parameters for component 'parabolic_optic_before_guide_v' [5]. */
MCNUM mccparabolic_optic_before_guide_v_sourceDist;
MCNUM mccparabolic_optic_before_guide_v_LStart;
MCNUM mccparabolic_optic_before_guide_v_LEnd;
MCNUM mccparabolic_optic_before_guide_v_lStart;
MCNUM mccparabolic_optic_before_guide_v_lEnd;
MCNUM mccparabolic_optic_before_guide_v_r_0;
MCNUM mccparabolic_optic_before_guide_v_nummirror;
MCNUM mccparabolic_optic_before_guide_v_mf;
MCNUM mccparabolic_optic_before_guide_v_mb;
MCNUM mccparabolic_optic_before_guide_v_mirror_width;
MCNUM mccparabolic_optic_before_guide_v_doubleReflections;

/* Definition parameters for component 'parabolic_optic_before_guide_h' [6]. */
#define mccparabolic_optic_before_guide_h_reflect "supermirror_m3.rfl" /* declared as a string. May produce warnings at compile */
/* Setting parameters for component 'parabolic_optic_before_guide_h' [6]. */
MCNUM mccparabolic_optic_before_guide_h_sourceDist;
MCNUM mccparabolic_optic_before_guide_h_LStart;
MCNUM mccparabolic_optic_before_guide_h_LEnd;
MCNUM mccparabolic_optic_before_guide_h_lStart;
MCNUM mccparabolic_optic_before_guide_h_lEnd;
MCNUM mccparabolic_optic_before_guide_h_r_0;
MCNUM mccparabolic_optic_before_guide_h_nummirror;
MCNUM mccparabolic_optic_before_guide_h_mf;
MCNUM mccparabolic_optic_before_guide_h_mb;
MCNUM mccparabolic_optic_before_guide_h_mirror_width;
MCNUM mccparabolic_optic_before_guide_h_doubleReflections;

/* Setting parameters for component 'psd_monitor_afteropticsource' [8]. */
int mccpsd_monitor_afteropticsource_nx;
int mccpsd_monitor_afteropticsource_ny;
char mccpsd_monitor_afteropticsource_filename[16384];
MCNUM mccpsd_monitor_afteropticsource_xmin;
MCNUM mccpsd_monitor_afteropticsource_xmax;
MCNUM mccpsd_monitor_afteropticsource_ymin;
MCNUM mccpsd_monitor_afteropticsource_ymax;
MCNUM mccpsd_monitor_afteropticsource_xwidth;
MCNUM mccpsd_monitor_afteropticsource_yheight;
MCNUM mccpsd_monitor_afteropticsource_restore_neutron;

/* Setting parameters for component 'guide_gravity_1' [9]. */
MCNUM mccguide_gravity_1_w1;
MCNUM mccguide_gravity_1_h1;
MCNUM mccguide_gravity_1_w2;
MCNUM mccguide_gravity_1_h2;
MCNUM mccguide_gravity_1_l;
MCNUM mccguide_gravity_1_R0;
MCNUM mccguide_gravity_1_Qc;
MCNUM mccguide_gravity_1_alpha;
MCNUM mccguide_gravity_1_m;
MCNUM mccguide_gravity_1_W;
MCNUM mccguide_gravity_1_nslit;
MCNUM mccguide_gravity_1_d;
MCNUM mccguide_gravity_1_mleft;
MCNUM mccguide_gravity_1_mright;
MCNUM mccguide_gravity_1_mtop;
MCNUM mccguide_gravity_1_mbottom;
MCNUM mccguide_gravity_1_nhslit;
MCNUM mccguide_gravity_1_G;
MCNUM mccguide_gravity_1_aleft;
MCNUM mccguide_gravity_1_aright;
MCNUM mccguide_gravity_1_atop;
MCNUM mccguide_gravity_1_abottom;
MCNUM mccguide_gravity_1_wavy;
MCNUM mccguide_gravity_1_wavy_z;
MCNUM mccguide_gravity_1_wavy_tb;
MCNUM mccguide_gravity_1_wavy_lr;
MCNUM mccguide_gravity_1_chamfers;
MCNUM mccguide_gravity_1_chamfers_z;
MCNUM mccguide_gravity_1_chamfers_lr;
MCNUM mccguide_gravity_1_chamfers_tb;
MCNUM mccguide_gravity_1_nelements;
MCNUM mccguide_gravity_1_nu;
MCNUM mccguide_gravity_1_phase;
char mccguide_gravity_1_reflect[16384];

/* Definition parameters for component 'divpos_monitor' [11]. */
#define mccdivpos_monitor_nh 100
#define mccdivpos_monitor_ndiv 100
/* Setting parameters for component 'divpos_monitor' [11]. */
char mccdivpos_monitor_filename[16384];
MCNUM mccdivpos_monitor_xmin;
MCNUM mccdivpos_monitor_xmax;
MCNUM mccdivpos_monitor_ymin;
MCNUM mccdivpos_monitor_ymax;
MCNUM mccdivpos_monitor_xwidth;
MCNUM mccdivpos_monitor_yheight;
MCNUM mccdivpos_monitor_maxdiv_h;
MCNUM mccdivpos_monitor_restore_neutron;
MCNUM mccdivpos_monitor_nx;
MCNUM mccdivpos_monitor_ny;
MCNUM mccdivpos_monitor_nz;
int mccdivpos_monitor_nowritefile;

/* Setting parameters for component 'psd_monitor_g1' [12]. */
int mccpsd_monitor_g1_nx;
int mccpsd_monitor_g1_ny;
char mccpsd_monitor_g1_filename[16384];
MCNUM mccpsd_monitor_g1_xmin;
MCNUM mccpsd_monitor_g1_xmax;
MCNUM mccpsd_monitor_g1_ymin;
MCNUM mccpsd_monitor_g1_ymax;
MCNUM mccpsd_monitor_g1_xwidth;
MCNUM mccpsd_monitor_g1_yheight;
MCNUM mccpsd_monitor_g1_restore_neutron;

/* Definition parameters for component 'divhlambda_monitor_g1' [13]. */
#define mccdivhlambda_monitor_g1_nL 100
#define mccdivhlambda_monitor_g1_nh 100
/* Setting parameters for component 'divhlambda_monitor_g1' [13]. */
char mccdivhlambda_monitor_g1_filename[16384];
MCNUM mccdivhlambda_monitor_g1_xmin;
MCNUM mccdivhlambda_monitor_g1_xmax;
MCNUM mccdivhlambda_monitor_g1_ymin;
MCNUM mccdivhlambda_monitor_g1_ymax;
MCNUM mccdivhlambda_monitor_g1_xwidth;
MCNUM mccdivhlambda_monitor_g1_yheight;
MCNUM mccdivhlambda_monitor_g1_maxdiv_h;
MCNUM mccdivhlambda_monitor_g1_Lmin;
MCNUM mccdivhlambda_monitor_g1_Lmax;
MCNUM mccdivhlambda_monitor_g1_restore_neutron;
MCNUM mccdivhlambda_monitor_g1_nx;
MCNUM mccdivhlambda_monitor_g1_ny;
MCNUM mccdivhlambda_monitor_g1_nz;
int mccdivhlambda_monitor_g1_nowritefile;

/* Definition parameters for component 'parabolic_optic1' [14]. */
#define mccparabolic_optic1_reflect "supermirror_m3.rfl" /* declared as a string. May produce warnings at compile */
/* Setting parameters for component 'parabolic_optic1' [14]. */
MCNUM mccparabolic_optic1_sourceDist;
MCNUM mccparabolic_optic1_LStart;
MCNUM mccparabolic_optic1_LEnd;
MCNUM mccparabolic_optic1_lStart;
MCNUM mccparabolic_optic1_lEnd;
MCNUM mccparabolic_optic1_r_0;
MCNUM mccparabolic_optic1_nummirror;
MCNUM mccparabolic_optic1_mf;
MCNUM mccparabolic_optic1_mb;
MCNUM mccparabolic_optic1_mirror_width;
MCNUM mccparabolic_optic1_doubleReflections;

/* Definition parameters for component 'parabolic_optic2' [15]. */
#define mccparabolic_optic2_reflect "supermirror_m3.rfl" /* declared as a string. May produce warnings at compile */
/* Setting parameters for component 'parabolic_optic2' [15]. */
MCNUM mccparabolic_optic2_sourceDist;
MCNUM mccparabolic_optic2_LStart;
MCNUM mccparabolic_optic2_LEnd;
MCNUM mccparabolic_optic2_lStart;
MCNUM mccparabolic_optic2_lEnd;
MCNUM mccparabolic_optic2_r_0;
MCNUM mccparabolic_optic2_nummirror;
MCNUM mccparabolic_optic2_mf;
MCNUM mccparabolic_optic2_mb;
MCNUM mccparabolic_optic2_mirror_width;
MCNUM mccparabolic_optic2_doubleReflections;

/* Setting parameters for component 'psd_monitor_f_zoom' [17]. */
int mccpsd_monitor_f_zoom_nx;
int mccpsd_monitor_f_zoom_ny;
char mccpsd_monitor_f_zoom_filename[16384];
MCNUM mccpsd_monitor_f_zoom_xmin;
MCNUM mccpsd_monitor_f_zoom_xmax;
MCNUM mccpsd_monitor_f_zoom_ymin;
MCNUM mccpsd_monitor_f_zoom_ymax;
MCNUM mccpsd_monitor_f_zoom_xwidth;
MCNUM mccpsd_monitor_f_zoom_yheight;
MCNUM mccpsd_monitor_f_zoom_restore_neutron;

/* Setting parameters for component 'psd_monitor_f' [18]. */
int mccpsd_monitor_f_nx;
int mccpsd_monitor_f_ny;
char mccpsd_monitor_f_filename[16384];
MCNUM mccpsd_monitor_f_xmin;
MCNUM mccpsd_monitor_f_xmax;
MCNUM mccpsd_monitor_f_ymin;
MCNUM mccpsd_monitor_f_ymax;
MCNUM mccpsd_monitor_f_xwidth;
MCNUM mccpsd_monitor_f_yheight;
MCNUM mccpsd_monitor_f_restore_neutron;

/* Definition parameters for component 'f_divpos' [19]. */
#define mccf_divpos_nh 100
#define mccf_divpos_ndiv 100
/* Setting parameters for component 'f_divpos' [19]. */
char mccf_divpos_filename[16384];
MCNUM mccf_divpos_xmin;
MCNUM mccf_divpos_xmax;
MCNUM mccf_divpos_ymin;
MCNUM mccf_divpos_ymax;
MCNUM mccf_divpos_xwidth;
MCNUM mccf_divpos_yheight;
MCNUM mccf_divpos_maxdiv_h;
MCNUM mccf_divpos_restore_neutron;
MCNUM mccf_divpos_nx;
MCNUM mccf_divpos_ny;
MCNUM mccf_divpos_nz;
int mccf_divpos_nowritefile;

/* Definition parameters for component 'divhlambda_monitor_f' [20]. */
#define mccdivhlambda_monitor_f_nL 100
#define mccdivhlambda_monitor_f_nh 100
/* Setting parameters for component 'divhlambda_monitor_f' [20]. */
char mccdivhlambda_monitor_f_filename[16384];
MCNUM mccdivhlambda_monitor_f_xmin;
MCNUM mccdivhlambda_monitor_f_xmax;
MCNUM mccdivhlambda_monitor_f_ymin;
MCNUM mccdivhlambda_monitor_f_ymax;
MCNUM mccdivhlambda_monitor_f_xwidth;
MCNUM mccdivhlambda_monitor_f_yheight;
MCNUM mccdivhlambda_monitor_f_maxdiv_h;
MCNUM mccdivhlambda_monitor_f_Lmin;
MCNUM mccdivhlambda_monitor_f_Lmax;
MCNUM mccdivhlambda_monitor_f_restore_neutron;
MCNUM mccdivhlambda_monitor_f_nx;
MCNUM mccdivhlambda_monitor_f_ny;
MCNUM mccdivhlambda_monitor_f_nz;
int mccdivhlambda_monitor_f_nowritefile;

/* User component declarations. */

/* User declarations for component 'origin' [1]. */
#define mccompcurname  origin
#define mccompcurtype  Progress_bar
#define mccompcurindex 1
#define IntermediateCnts mccorigin_IntermediateCnts
#define StartTime mccorigin_StartTime
#define EndTime mccorigin_EndTime
#define CurrentTime mccorigin_CurrentTime
#define profile mccorigin_profile
#define percent mccorigin_percent
#define flag_save mccorigin_flag_save
#define minutes mccorigin_minutes
#line 44 "/usr/share/mcstas/2.6.1/tools/Python/mcrun/../mccodelib/../../../misc/Progress_bar.comp"
#ifndef PROGRESS_BAR
#define PROGRESS_BAR
#else
#error Only one Progress_bar component may be used in an instrument definition.
#endif

double IntermediateCnts;
time_t StartTime;
time_t EndTime;
time_t CurrentTime;
#line 9383 "./ess_extraction.c"
#undef minutes
#undef flag_save
#undef percent
#undef profile
#undef CurrentTime
#undef EndTime
#undef StartTime
#undef IntermediateCnts
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'source' [2]. */
#define mccompcurname  source
#define mccompcurtype  Arm
#define mccompcurindex 2
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'source_div' [3]. */
#define mccompcurname  source_div
#define mccompcurtype  Source_div
#define mccompcurindex 3
#define thetah mccsource_div_thetah
#define thetav mccsource_div_thetav
#define sigmah mccsource_div_sigmah
#define sigmav mccsource_div_sigmav
#define tan_h mccsource_div_tan_h
#define tan_v mccsource_div_tan_v
#define p_init mccsource_div_p_init
#define dist mccsource_div_dist
#define focus_xw mccsource_div_focus_xw
#define focus_yh mccsource_div_focus_yh
#define xwidth mccsource_div_xwidth
#define yheight mccsource_div_yheight
#define focus_aw mccsource_div_focus_aw
#define focus_ah mccsource_div_focus_ah
#define E0 mccsource_div_E0
#define dE mccsource_div_dE
#define lambda0 mccsource_div_lambda0
#define dlambda mccsource_div_dlambda
#define gauss mccsource_div_gauss
#define flux mccsource_div_flux
#line 69 "/usr/share/mcstas/2.6.1/tools/Python/mcrun/../mccodelib/../../../sources/Source_div.comp"
double thetah, thetav, sigmah, sigmav, tan_h, tan_v, p_init, dist, focus_xw, focus_yh;
#line 9430 "./ess_extraction.c"
#undef flux
#undef gauss
#undef dlambda
#undef lambda0
#undef dE
#undef E0
#undef focus_ah
#undef focus_aw
#undef yheight
#undef xwidth
#undef focus_yh
#undef focus_xw
#undef dist
#undef p_init
#undef tan_v
#undef tan_h
#undef sigmav
#undef sigmah
#undef thetav
#undef thetah
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'psd_monitor_source' [4]. */
#define mccompcurname  psd_monitor_source
#define mccompcurtype  PSD_monitor
#define mccompcurindex 4
#define PSD_N mccpsd_monitor_source_PSD_N
#define PSD_p mccpsd_monitor_source_PSD_p
#define PSD_p2 mccpsd_monitor_source_PSD_p2
#define nx mccpsd_monitor_source_nx
#define ny mccpsd_monitor_source_ny
#define filename mccpsd_monitor_source_filename
#define xmin mccpsd_monitor_source_xmin
#define xmax mccpsd_monitor_source_xmax
#define ymin mccpsd_monitor_source_ymin
#define ymax mccpsd_monitor_source_ymax
#define xwidth mccpsd_monitor_source_xwidth
#define yheight mccpsd_monitor_source_yheight
#define restore_neutron mccpsd_monitor_source_restore_neutron
#line 62 "/usr/share/mcstas/2.6.1/tools/Python/mcrun/../mccodelib/../../../monitors/PSD_monitor.comp"
  DArray2d PSD_N;
  DArray2d PSD_p;
  DArray2d PSD_p2;
#line 9476 "./ess_extraction.c"
#undef restore_neutron
#undef yheight
#undef xwidth
#undef ymax
#undef ymin
#undef xmax
#undef xmin
#undef filename
#undef ny
#undef nx
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'parabolic_optic_before_guide_v' [5]. */
#define mccompcurname  parabolic_optic_before_guide_v
#define mccompcurtype  FlatEllipse_finite_mirror
#define mccompcurindex 5
#define reflect mccparabolic_optic_before_guide_v_reflect
#define s mccparabolic_optic_before_guide_v_s
#define pTable mccparabolic_optic_before_guide_v_pTable
#define R0 mccparabolic_optic_before_guide_v_R0
#define Qc mccparabolic_optic_before_guide_v_Qc
#define W mccparabolic_optic_before_guide_v_W
#define alpha mccparabolic_optic_before_guide_v_alpha
#define transmit mccparabolic_optic_before_guide_v_transmit
#define sourceDist mccparabolic_optic_before_guide_v_sourceDist
#define LStart mccparabolic_optic_before_guide_v_LStart
#define LEnd mccparabolic_optic_before_guide_v_LEnd
#define lStart mccparabolic_optic_before_guide_v_lStart
#define lEnd mccparabolic_optic_before_guide_v_lEnd
#define r_0 mccparabolic_optic_before_guide_v_r_0
#define nummirror mccparabolic_optic_before_guide_v_nummirror
#define mf mccparabolic_optic_before_guide_v_mf
#define mb mccparabolic_optic_before_guide_v_mb
#define mirror_width mccparabolic_optic_before_guide_v_mirror_width
#define doubleReflections mccparabolic_optic_before_guide_v_doubleReflections
#line 61 "FlatEllipse_finite_mirror.comp"
    //Scene where all geometry is added to
    Scene s;
    //Variables for Reflectivity from McStas Tables
    //t_Table pTable;
    //double R0 = 0.99;
    //double Qc = 0.021;
    //double W = 0.003;
    //double alpha = 6.07;
    //double transmit = 0;
    double *pointer_lStart;// lStart has to be used in Trace later, this requires a pointer
    //Function to handle Conic-Neutron collisions with reflectivity from McStas Tables
    void traceNeutronConicWithTables(Particle* p, ConicSurf c);
    double *rs;
    double dt;
    int silicon; // +1: neutron in silicon, -1: neutron in air, 0: mirrorwidth is 0; neutron cannot be in silicon
#line 9533 "./ess_extraction.c"
#undef doubleReflections
#undef mirror_width
#undef mb
#undef mf
#undef nummirror
#undef r_0
#undef lEnd
#undef lStart
#undef LEnd
#undef LStart
#undef sourceDist
#undef transmit
#undef alpha
#undef W
#undef Qc
#undef R0
#undef pTable
#undef s
#undef reflect
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'parabolic_optic_before_guide_h' [6]. */
#define mccompcurname  parabolic_optic_before_guide_h
#define mccompcurtype  FlatEllipse_finite_mirror
#define mccompcurindex 6
#define reflect mccparabolic_optic_before_guide_h_reflect
#define s mccparabolic_optic_before_guide_h_s
#define pTable mccparabolic_optic_before_guide_h_pTable
#define R0 mccparabolic_optic_before_guide_h_R0
#define Qc mccparabolic_optic_before_guide_h_Qc
#define W mccparabolic_optic_before_guide_h_W
#define alpha mccparabolic_optic_before_guide_h_alpha
#define transmit mccparabolic_optic_before_guide_h_transmit
#define sourceDist mccparabolic_optic_before_guide_h_sourceDist
#define LStart mccparabolic_optic_before_guide_h_LStart
#define LEnd mccparabolic_optic_before_guide_h_LEnd
#define lStart mccparabolic_optic_before_guide_h_lStart
#define lEnd mccparabolic_optic_before_guide_h_lEnd
#define r_0 mccparabolic_optic_before_guide_h_r_0
#define nummirror mccparabolic_optic_before_guide_h_nummirror
#define mf mccparabolic_optic_before_guide_h_mf
#define mb mccparabolic_optic_before_guide_h_mb
#define mirror_width mccparabolic_optic_before_guide_h_mirror_width
#define doubleReflections mccparabolic_optic_before_guide_h_doubleReflections
#line 61 "FlatEllipse_finite_mirror.comp"
    //Scene where all geometry is added to
    Scene s;
    //Variables for Reflectivity from McStas Tables
    //t_Table pTable;
    //double R0 = 0.99;
    //double Qc = 0.021;
    //double W = 0.003;
    //double alpha = 6.07;
    //double transmit = 0;
    double *pointer_lStart;// lStart has to be used in Trace later, this requires a pointer
    //Function to handle Conic-Neutron collisions with reflectivity from McStas Tables
    void traceNeutronConicWithTables(Particle* p, ConicSurf c);
    double *rs;
    double dt;
    int silicon; // +1: neutron in silicon, -1: neutron in air, 0: mirrorwidth is 0; neutron cannot be in silicon
#line 9596 "./ess_extraction.c"
#undef doubleReflections
#undef mirror_width
#undef mb
#undef mf
#undef nummirror
#undef r_0
#undef lEnd
#undef lStart
#undef LEnd
#undef LStart
#undef sourceDist
#undef transmit
#undef alpha
#undef W
#undef Qc
#undef R0
#undef pTable
#undef s
#undef reflect
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'after_optic_source' [7]. */
#define mccompcurname  after_optic_source
#define mccompcurtype  Arm
#define mccompcurindex 7
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'psd_monitor_afteropticsource' [8]. */
#define mccompcurname  psd_monitor_afteropticsource
#define mccompcurtype  PSD_monitor
#define mccompcurindex 8
#define PSD_N mccpsd_monitor_afteropticsource_PSD_N
#define PSD_p mccpsd_monitor_afteropticsource_PSD_p
#define PSD_p2 mccpsd_monitor_afteropticsource_PSD_p2
#define nx mccpsd_monitor_afteropticsource_nx
#define ny mccpsd_monitor_afteropticsource_ny
#define filename mccpsd_monitor_afteropticsource_filename
#define xmin mccpsd_monitor_afteropticsource_xmin
#define xmax mccpsd_monitor_afteropticsource_xmax
#define ymin mccpsd_monitor_afteropticsource_ymin
#define ymax mccpsd_monitor_afteropticsource_ymax
#define xwidth mccpsd_monitor_afteropticsource_xwidth
#define yheight mccpsd_monitor_afteropticsource_yheight
#define restore_neutron mccpsd_monitor_afteropticsource_restore_neutron
#line 62 "/usr/share/mcstas/2.6.1/tools/Python/mcrun/../mccodelib/../../../monitors/PSD_monitor.comp"
  DArray2d PSD_N;
  DArray2d PSD_p;
  DArray2d PSD_p2;
#line 9649 "./ess_extraction.c"
#undef restore_neutron
#undef yheight
#undef xwidth
#undef ymax
#undef ymin
#undef xmax
#undef xmin
#undef filename
#undef ny
#undef nx
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'guide_gravity_1' [9]. */
#define mccompcurname  guide_gravity_1
#define mccompcurtype  Guide_gravity
#define mccompcurindex 9
#define GVars mccguide_gravity_1_GVars
#define pTable mccguide_gravity_1_pTable
#define w1 mccguide_gravity_1_w1
#define h1 mccguide_gravity_1_h1
#define w2 mccguide_gravity_1_w2
#define h2 mccguide_gravity_1_h2
#define l mccguide_gravity_1_l
#define R0 mccguide_gravity_1_R0
#define Qc mccguide_gravity_1_Qc
#define alpha mccguide_gravity_1_alpha
#define m mccguide_gravity_1_m
#define W mccguide_gravity_1_W
#define nslit mccguide_gravity_1_nslit
#define d mccguide_gravity_1_d
#define mleft mccguide_gravity_1_mleft
#define mright mccguide_gravity_1_mright
#define mtop mccguide_gravity_1_mtop
#define mbottom mccguide_gravity_1_mbottom
#define nhslit mccguide_gravity_1_nhslit
#define G mccguide_gravity_1_G
#define aleft mccguide_gravity_1_aleft
#define aright mccguide_gravity_1_aright
#define atop mccguide_gravity_1_atop
#define abottom mccguide_gravity_1_abottom
#define wavy mccguide_gravity_1_wavy
#define wavy_z mccguide_gravity_1_wavy_z
#define wavy_tb mccguide_gravity_1_wavy_tb
#define wavy_lr mccguide_gravity_1_wavy_lr
#define chamfers mccguide_gravity_1_chamfers
#define chamfers_z mccguide_gravity_1_chamfers_z
#define chamfers_lr mccguide_gravity_1_chamfers_lr
#define chamfers_tb mccguide_gravity_1_chamfers_tb
#define nelements mccguide_gravity_1_nelements
#define nu mccguide_gravity_1_nu
#define phase mccguide_gravity_1_phase
#define reflect mccguide_gravity_1_reflect
#line 334 "/usr/share/mcstas/2.6.1/tools/Python/mcrun/../mccodelib/../../../optics/Guide_gravity.comp"
  Gravity_guide_Vars_type GVars;
  t_Table pTable;
#line 9710 "./ess_extraction.c"
#undef reflect
#undef phase
#undef nu
#undef nelements
#undef chamfers_tb
#undef chamfers_lr
#undef chamfers_z
#undef chamfers
#undef wavy_lr
#undef wavy_tb
#undef wavy_z
#undef wavy
#undef abottom
#undef atop
#undef aright
#undef aleft
#undef G
#undef nhslit
#undef mbottom
#undef mtop
#undef mright
#undef mleft
#undef d
#undef nslit
#undef W
#undef m
#undef alpha
#undef Qc
#undef R0
#undef l
#undef h2
#undef w2
#undef h1
#undef w1
#undef pTable
#undef GVars
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'monitor_1' [10]. */
#define mccompcurname  monitor_1
#define mccompcurtype  Arm
#define mccompcurindex 10
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'divpos_monitor' [11]. */
#define mccompcurname  divpos_monitor
#define mccompcurtype  DivPos_monitor
#define mccompcurindex 11
#define nh mccdivpos_monitor_nh
#define ndiv mccdivpos_monitor_ndiv
#define Div_N mccdivpos_monitor_Div_N
#define Div_p mccdivpos_monitor_Div_p
#define Div_p2 mccdivpos_monitor_Div_p2
#define filename mccdivpos_monitor_filename
#define xmin mccdivpos_monitor_xmin
#define xmax mccdivpos_monitor_xmax
#define ymin mccdivpos_monitor_ymin
#define ymax mccdivpos_monitor_ymax
#define xwidth mccdivpos_monitor_xwidth
#define yheight mccdivpos_monitor_yheight
#define maxdiv_h mccdivpos_monitor_maxdiv_h
#define restore_neutron mccdivpos_monitor_restore_neutron
#define nx mccdivpos_monitor_nx
#define ny mccdivpos_monitor_ny
#define nz mccdivpos_monitor_nz
#define nowritefile mccdivpos_monitor_nowritefile
#line 62 "/usr/share/mcstas/2.6.1/tools/Python/mcrun/../mccodelib/../../../monitors/DivPos_monitor.comp"
double Div_N[nh][ndiv];
double Div_p[nh][ndiv];
double Div_p2[nh][ndiv];
#line 9785 "./ess_extraction.c"
#undef nowritefile
#undef nz
#undef ny
#undef nx
#undef restore_neutron
#undef maxdiv_h
#undef yheight
#undef xwidth
#undef ymax
#undef ymin
#undef xmax
#undef xmin
#undef filename
#undef Div_p2
#undef Div_p
#undef Div_N
#undef ndiv
#undef nh
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'psd_monitor_g1' [12]. */
#define mccompcurname  psd_monitor_g1
#define mccompcurtype  PSD_monitor
#define mccompcurindex 12
#define PSD_N mccpsd_monitor_g1_PSD_N
#define PSD_p mccpsd_monitor_g1_PSD_p
#define PSD_p2 mccpsd_monitor_g1_PSD_p2
#define nx mccpsd_monitor_g1_nx
#define ny mccpsd_monitor_g1_ny
#define filename mccpsd_monitor_g1_filename
#define xmin mccpsd_monitor_g1_xmin
#define xmax mccpsd_monitor_g1_xmax
#define ymin mccpsd_monitor_g1_ymin
#define ymax mccpsd_monitor_g1_ymax
#define xwidth mccpsd_monitor_g1_xwidth
#define yheight mccpsd_monitor_g1_yheight
#define restore_neutron mccpsd_monitor_g1_restore_neutron
#line 62 "/usr/share/mcstas/2.6.1/tools/Python/mcrun/../mccodelib/../../../monitors/PSD_monitor.comp"
  DArray2d PSD_N;
  DArray2d PSD_p;
  DArray2d PSD_p2;
#line 9829 "./ess_extraction.c"
#undef restore_neutron
#undef yheight
#undef xwidth
#undef ymax
#undef ymin
#undef xmax
#undef xmin
#undef filename
#undef ny
#undef nx
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'divhlambda_monitor_g1' [13]. */
#define mccompcurname  divhlambda_monitor_g1
#define mccompcurtype  DivLambda_monitor
#define mccompcurindex 13
#define nL mccdivhlambda_monitor_g1_nL
#define nh mccdivhlambda_monitor_g1_nh
#define Div_N mccdivhlambda_monitor_g1_Div_N
#define Div_p mccdivhlambda_monitor_g1_Div_p
#define Div_p2 mccdivhlambda_monitor_g1_Div_p2
#define filename mccdivhlambda_monitor_g1_filename
#define xmin mccdivhlambda_monitor_g1_xmin
#define xmax mccdivhlambda_monitor_g1_xmax
#define ymin mccdivhlambda_monitor_g1_ymin
#define ymax mccdivhlambda_monitor_g1_ymax
#define xwidth mccdivhlambda_monitor_g1_xwidth
#define yheight mccdivhlambda_monitor_g1_yheight
#define maxdiv_h mccdivhlambda_monitor_g1_maxdiv_h
#define Lmin mccdivhlambda_monitor_g1_Lmin
#define Lmax mccdivhlambda_monitor_g1_Lmax
#define restore_neutron mccdivhlambda_monitor_g1_restore_neutron
#define nx mccdivhlambda_monitor_g1_nx
#define ny mccdivhlambda_monitor_g1_ny
#define nz mccdivhlambda_monitor_g1_nz
#define nowritefile mccdivhlambda_monitor_g1_nowritefile
#line 62 "/usr/share/mcstas/2.6.1/tools/Python/mcrun/../mccodelib/../../../monitors/DivLambda_monitor.comp"
double Div_N[nL][nh];
double Div_p[nL][nh];
double Div_p2[nL][nh];
#line 9875 "./ess_extraction.c"
#undef nowritefile
#undef nz
#undef ny
#undef nx
#undef restore_neutron
#undef Lmax
#undef Lmin
#undef maxdiv_h
#undef yheight
#undef xwidth
#undef ymax
#undef ymin
#undef xmax
#undef xmin
#undef filename
#undef Div_p2
#undef Div_p
#undef Div_N
#undef nh
#undef nL
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'parabolic_optic1' [14]. */
#define mccompcurname  parabolic_optic1
#define mccompcurtype  FlatEllipse_finite_mirror
#define mccompcurindex 14
#define reflect mccparabolic_optic1_reflect
#define s mccparabolic_optic1_s
#define pTable mccparabolic_optic1_pTable
#define R0 mccparabolic_optic1_R0
#define Qc mccparabolic_optic1_Qc
#define W mccparabolic_optic1_W
#define alpha mccparabolic_optic1_alpha
#define transmit mccparabolic_optic1_transmit
#define sourceDist mccparabolic_optic1_sourceDist
#define LStart mccparabolic_optic1_LStart
#define LEnd mccparabolic_optic1_LEnd
#define lStart mccparabolic_optic1_lStart
#define lEnd mccparabolic_optic1_lEnd
#define r_0 mccparabolic_optic1_r_0
#define nummirror mccparabolic_optic1_nummirror
#define mf mccparabolic_optic1_mf
#define mb mccparabolic_optic1_mb
#define mirror_width mccparabolic_optic1_mirror_width
#define doubleReflections mccparabolic_optic1_doubleReflections
#line 61 "FlatEllipse_finite_mirror.comp"
    //Scene where all geometry is added to
    Scene s;
    //Variables for Reflectivity from McStas Tables
    //t_Table pTable;
    //double R0 = 0.99;
    //double Qc = 0.021;
    //double W = 0.003;
    //double alpha = 6.07;
    //double transmit = 0;
    double *pointer_lStart;// lStart has to be used in Trace later, this requires a pointer
    //Function to handle Conic-Neutron collisions with reflectivity from McStas Tables
    void traceNeutronConicWithTables(Particle* p, ConicSurf c);
    double *rs;
    double dt;
    int silicon; // +1: neutron in silicon, -1: neutron in air, 0: mirrorwidth is 0; neutron cannot be in silicon
#line 9939 "./ess_extraction.c"
#undef doubleReflections
#undef mirror_width
#undef mb
#undef mf
#undef nummirror
#undef r_0
#undef lEnd
#undef lStart
#undef LEnd
#undef LStart
#undef sourceDist
#undef transmit
#undef alpha
#undef W
#undef Qc
#undef R0
#undef pTable
#undef s
#undef reflect
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'parabolic_optic2' [15]. */
#define mccompcurname  parabolic_optic2
#define mccompcurtype  FlatEllipse_finite_mirror
#define mccompcurindex 15
#define reflect mccparabolic_optic2_reflect
#define s mccparabolic_optic2_s
#define pTable mccparabolic_optic2_pTable
#define R0 mccparabolic_optic2_R0
#define Qc mccparabolic_optic2_Qc
#define W mccparabolic_optic2_W
#define alpha mccparabolic_optic2_alpha
#define transmit mccparabolic_optic2_transmit
#define sourceDist mccparabolic_optic2_sourceDist
#define LStart mccparabolic_optic2_LStart
#define LEnd mccparabolic_optic2_LEnd
#define lStart mccparabolic_optic2_lStart
#define lEnd mccparabolic_optic2_lEnd
#define r_0 mccparabolic_optic2_r_0
#define nummirror mccparabolic_optic2_nummirror
#define mf mccparabolic_optic2_mf
#define mb mccparabolic_optic2_mb
#define mirror_width mccparabolic_optic2_mirror_width
#define doubleReflections mccparabolic_optic2_doubleReflections
#line 61 "FlatEllipse_finite_mirror.comp"
    //Scene where all geometry is added to
    Scene s;
    //Variables for Reflectivity from McStas Tables
    //t_Table pTable;
    //double R0 = 0.99;
    //double Qc = 0.021;
    //double W = 0.003;
    //double alpha = 6.07;
    //double transmit = 0;
    double *pointer_lStart;// lStart has to be used in Trace later, this requires a pointer
    //Function to handle Conic-Neutron collisions with reflectivity from McStas Tables
    void traceNeutronConicWithTables(Particle* p, ConicSurf c);
    double *rs;
    double dt;
    int silicon; // +1: neutron in silicon, -1: neutron in air, 0: mirrorwidth is 0; neutron cannot be in silicon
#line 10002 "./ess_extraction.c"
#undef doubleReflections
#undef mirror_width
#undef mb
#undef mf
#undef nummirror
#undef r_0
#undef lEnd
#undef lStart
#undef LEnd
#undef LStart
#undef sourceDist
#undef transmit
#undef alpha
#undef W
#undef Qc
#undef R0
#undef pTable
#undef s
#undef reflect
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'monitor_2' [16]. */
#define mccompcurname  monitor_2
#define mccompcurtype  Arm
#define mccompcurindex 16
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'psd_monitor_f_zoom' [17]. */
#define mccompcurname  psd_monitor_f_zoom
#define mccompcurtype  PSD_monitor
#define mccompcurindex 17
#define PSD_N mccpsd_monitor_f_zoom_PSD_N
#define PSD_p mccpsd_monitor_f_zoom_PSD_p
#define PSD_p2 mccpsd_monitor_f_zoom_PSD_p2
#define nx mccpsd_monitor_f_zoom_nx
#define ny mccpsd_monitor_f_zoom_ny
#define filename mccpsd_monitor_f_zoom_filename
#define xmin mccpsd_monitor_f_zoom_xmin
#define xmax mccpsd_monitor_f_zoom_xmax
#define ymin mccpsd_monitor_f_zoom_ymin
#define ymax mccpsd_monitor_f_zoom_ymax
#define xwidth mccpsd_monitor_f_zoom_xwidth
#define yheight mccpsd_monitor_f_zoom_yheight
#define restore_neutron mccpsd_monitor_f_zoom_restore_neutron
#line 62 "/usr/share/mcstas/2.6.1/tools/Python/mcrun/../mccodelib/../../../monitors/PSD_monitor.comp"
  DArray2d PSD_N;
  DArray2d PSD_p;
  DArray2d PSD_p2;
#line 10055 "./ess_extraction.c"
#undef restore_neutron
#undef yheight
#undef xwidth
#undef ymax
#undef ymin
#undef xmax
#undef xmin
#undef filename
#undef ny
#undef nx
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'psd_monitor_f' [18]. */
#define mccompcurname  psd_monitor_f
#define mccompcurtype  PSD_monitor
#define mccompcurindex 18
#define PSD_N mccpsd_monitor_f_PSD_N
#define PSD_p mccpsd_monitor_f_PSD_p
#define PSD_p2 mccpsd_monitor_f_PSD_p2
#define nx mccpsd_monitor_f_nx
#define ny mccpsd_monitor_f_ny
#define filename mccpsd_monitor_f_filename
#define xmin mccpsd_monitor_f_xmin
#define xmax mccpsd_monitor_f_xmax
#define ymin mccpsd_monitor_f_ymin
#define ymax mccpsd_monitor_f_ymax
#define xwidth mccpsd_monitor_f_xwidth
#define yheight mccpsd_monitor_f_yheight
#define restore_neutron mccpsd_monitor_f_restore_neutron
#line 62 "/usr/share/mcstas/2.6.1/tools/Python/mcrun/../mccodelib/../../../monitors/PSD_monitor.comp"
  DArray2d PSD_N;
  DArray2d PSD_p;
  DArray2d PSD_p2;
#line 10094 "./ess_extraction.c"
#undef restore_neutron
#undef yheight
#undef xwidth
#undef ymax
#undef ymin
#undef xmax
#undef xmin
#undef filename
#undef ny
#undef nx
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'f_divpos' [19]. */
#define mccompcurname  f_divpos
#define mccompcurtype  DivPos_monitor
#define mccompcurindex 19
#define nh mccf_divpos_nh
#define ndiv mccf_divpos_ndiv
#define Div_N mccf_divpos_Div_N
#define Div_p mccf_divpos_Div_p
#define Div_p2 mccf_divpos_Div_p2
#define filename mccf_divpos_filename
#define xmin mccf_divpos_xmin
#define xmax mccf_divpos_xmax
#define ymin mccf_divpos_ymin
#define ymax mccf_divpos_ymax
#define xwidth mccf_divpos_xwidth
#define yheight mccf_divpos_yheight
#define maxdiv_h mccf_divpos_maxdiv_h
#define restore_neutron mccf_divpos_restore_neutron
#define nx mccf_divpos_nx
#define ny mccf_divpos_ny
#define nz mccf_divpos_nz
#define nowritefile mccf_divpos_nowritefile
#line 62 "/usr/share/mcstas/2.6.1/tools/Python/mcrun/../mccodelib/../../../monitors/DivPos_monitor.comp"
double Div_N[nh][ndiv];
double Div_p[nh][ndiv];
double Div_p2[nh][ndiv];
#line 10138 "./ess_extraction.c"
#undef nowritefile
#undef nz
#undef ny
#undef nx
#undef restore_neutron
#undef maxdiv_h
#undef yheight
#undef xwidth
#undef ymax
#undef ymin
#undef xmax
#undef xmin
#undef filename
#undef Div_p2
#undef Div_p
#undef Div_N
#undef ndiv
#undef nh
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'divhlambda_monitor_f' [20]. */
#define mccompcurname  divhlambda_monitor_f
#define mccompcurtype  DivLambda_monitor
#define mccompcurindex 20
#define nL mccdivhlambda_monitor_f_nL
#define nh mccdivhlambda_monitor_f_nh
#define Div_N mccdivhlambda_monitor_f_Div_N
#define Div_p mccdivhlambda_monitor_f_Div_p
#define Div_p2 mccdivhlambda_monitor_f_Div_p2
#define filename mccdivhlambda_monitor_f_filename
#define xmin mccdivhlambda_monitor_f_xmin
#define xmax mccdivhlambda_monitor_f_xmax
#define ymin mccdivhlambda_monitor_f_ymin
#define ymax mccdivhlambda_monitor_f_ymax
#define xwidth mccdivhlambda_monitor_f_xwidth
#define yheight mccdivhlambda_monitor_f_yheight
#define maxdiv_h mccdivhlambda_monitor_f_maxdiv_h
#define Lmin mccdivhlambda_monitor_f_Lmin
#define Lmax mccdivhlambda_monitor_f_Lmax
#define restore_neutron mccdivhlambda_monitor_f_restore_neutron
#define nx mccdivhlambda_monitor_f_nx
#define ny mccdivhlambda_monitor_f_ny
#define nz mccdivhlambda_monitor_f_nz
#define nowritefile mccdivhlambda_monitor_f_nowritefile
#line 62 "/usr/share/mcstas/2.6.1/tools/Python/mcrun/../mccodelib/../../../monitors/DivLambda_monitor.comp"
double Div_N[nL][nh];
double Div_p[nL][nh];
double Div_p2[nL][nh];
#line 10189 "./ess_extraction.c"
#undef nowritefile
#undef nz
#undef ny
#undef nx
#undef restore_neutron
#undef Lmax
#undef Lmin
#undef maxdiv_h
#undef yheight
#undef xwidth
#undef ymax
#undef ymin
#undef xmax
#undef xmin
#undef filename
#undef Div_p2
#undef Div_p
#undef Div_N
#undef nh
#undef nL
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

Coords mcposaorigin, mcposrorigin;
Rotation mcrotaorigin, mcrotrorigin;
Coords mcposasource, mcposrsource;
Rotation mcrotasource, mcrotrsource;
Coords mcposasource_div, mcposrsource_div;
Rotation mcrotasource_div, mcrotrsource_div;
Coords mcposapsd_monitor_source, mcposrpsd_monitor_source;
Rotation mcrotapsd_monitor_source, mcrotrpsd_monitor_source;
Coords mcposaparabolic_optic_before_guide_v, mcposrparabolic_optic_before_guide_v;
Rotation mcrotaparabolic_optic_before_guide_v, mcrotrparabolic_optic_before_guide_v;
Coords mcposaparabolic_optic_before_guide_h, mcposrparabolic_optic_before_guide_h;
Rotation mcrotaparabolic_optic_before_guide_h, mcrotrparabolic_optic_before_guide_h;
Coords mcposaafter_optic_source, mcposrafter_optic_source;
Rotation mcrotaafter_optic_source, mcrotrafter_optic_source;
Coords mcposapsd_monitor_afteropticsource, mcposrpsd_monitor_afteropticsource;
Rotation mcrotapsd_monitor_afteropticsource, mcrotrpsd_monitor_afteropticsource;
Coords mcposaguide_gravity_1, mcposrguide_gravity_1;
Rotation mcrotaguide_gravity_1, mcrotrguide_gravity_1;
Coords mcposamonitor_1, mcposrmonitor_1;
Rotation mcrotamonitor_1, mcrotrmonitor_1;
Coords mcposadivpos_monitor, mcposrdivpos_monitor;
Rotation mcrotadivpos_monitor, mcrotrdivpos_monitor;
Coords mcposapsd_monitor_g1, mcposrpsd_monitor_g1;
Rotation mcrotapsd_monitor_g1, mcrotrpsd_monitor_g1;
Coords mcposadivhlambda_monitor_g1, mcposrdivhlambda_monitor_g1;
Rotation mcrotadivhlambda_monitor_g1, mcrotrdivhlambda_monitor_g1;
Coords mcposaparabolic_optic1, mcposrparabolic_optic1;
Rotation mcrotaparabolic_optic1, mcrotrparabolic_optic1;
Coords mcposaparabolic_optic2, mcposrparabolic_optic2;
Rotation mcrotaparabolic_optic2, mcrotrparabolic_optic2;
Coords mcposamonitor_2, mcposrmonitor_2;
Rotation mcrotamonitor_2, mcrotrmonitor_2;
Coords mcposapsd_monitor_f_zoom, mcposrpsd_monitor_f_zoom;
Rotation mcrotapsd_monitor_f_zoom, mcrotrpsd_monitor_f_zoom;
Coords mcposapsd_monitor_f, mcposrpsd_monitor_f;
Rotation mcrotapsd_monitor_f, mcrotrpsd_monitor_f;
Coords mcposaf_divpos, mcposrf_divpos;
Rotation mcrotaf_divpos, mcrotrf_divpos;
Coords mcposadivhlambda_monitor_f, mcposrdivhlambda_monitor_f;
Rotation mcrotadivhlambda_monitor_f, mcrotrdivhlambda_monitor_f;

MCNUM mcnx, mcny, mcnz, mcnvx, mcnvy, mcnvz, mcnt, mcnsx, mcnsy, mcnsz, mcnp;

/* end declare */

void mcinit(void) {
#define mccompcurname  template_simple
#define mccompcurtype  INSTRUMENT
#define mccompcurindex 0
#define mcposatemplate_simple coords_set(0,0,0)
#define det_width mcipdet_width
#define det_width_focus mcipdet_width_focus
#define source_width mcipsource_width
#define guide_width mcipguide_width
#define L_source mcipL_source
#define guide_length mcipguide_length
#define dL mcipdL
#define L_min mcipL_min
#define L_max mcipL_max
#define source_divergence mcipsource_divergence
#define divergence_max mcipdivergence_max
#define focal_length mcipfocal_length
#define mirrors mcipmirrors
#define incoming_length mcipincoming_length
#define g mcipg
#define max_div mcipmax_div
#define pixels mcippixels
#define flux mcipflux
#define placeholder mcipplaceholder
#undef placeholder
#undef flux
#undef pixels
#undef max_div
#undef g
#undef incoming_length
#undef mirrors
#undef focal_length
#undef divergence_max
#undef source_divergence
#undef L_max
#undef L_min
#undef dL
#undef guide_length
#undef L_source
#undef guide_width
#undef source_width
#undef det_width_focus
#undef det_width
#undef mcposatemplate_simple
#undef mccompcurindex
#undef mccompcurtype
#undef mccompcurname
  /* Computation of coordinate transformations. */
  {
    Coords mctc1, mctc2, mcLastComp;
    Rotation mctr1;
    double mcAccumulatedILength = 0;
    /* Initialize "last" component origin as (0,0,0) */
    mcLastComp = coords_set(0,0,0);

    mcDEBUG_INSTR()
  /* Component initializations. */
    /* Component origin. */
  /* Setting parameters for component origin. */
  SIG_MESSAGE("origin (Init:SetPar)");
#line 39 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  if("NULL") strncpy(mccorigin_profile, "NULL" ? "NULL" : "", 16384); else mccorigin_profile[0]='\0';
#line 39 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccorigin_percent = 10;
#line 39 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccorigin_flag_save = 0;
#line 39 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccorigin_minutes = 0;
#line 10327 "./ess_extraction.c"

  SIG_MESSAGE("origin (Init:Place/Rotate)");
  rot_set_rotation(mcrotaorigin,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 10334 "./ess_extraction.c"
  rot_copy(mcrotrorigin, mcrotaorigin);
  mcposaorigin = coords_set(
#line 59 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
    0,
#line 59 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
    0,
#line 59 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
    0);
#line 10343 "./ess_extraction.c"
  mctc1 = coords_neg(mcposaorigin);
  mcposrorigin = rot_apply(mcrotaorigin, mctc1);
  mcDEBUG_COMPONENT("origin", mcposaorigin, mcrotaorigin)
  mccomp_posa[1] = mcposaorigin;
  mccomp_posr[1] = mcposrorigin;
  mcNCounter[1]  = mcPCounter[1] = mcP2Counter[1] = 0;
  mcAbsorbProp[1]= 0;
    /* Component source. */
  /* Setting parameters for component source. */
  SIG_MESSAGE("source (Init:SetPar)");

  SIG_MESSAGE("source (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 10360 "./ess_extraction.c"
  rot_mul(mctr1, mcrotaorigin, mcrotasource);
  rot_transpose(mcrotaorigin, mctr1);
  rot_mul(mcrotasource, mctr1, mcrotrsource);
  mctc1 = coords_set(
#line 63 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
    0,
#line 63 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
    0,
#line 63 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
    0);
#line 10371 "./ess_extraction.c"
  rot_transpose(mcrotaorigin, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposasource = coords_add(mcposaorigin, mctc2);
  mctc1 = coords_sub(mcposaorigin, mcposasource);
  mcposrsource = rot_apply(mcrotasource, mctc1);
  mcDEBUG_COMPONENT("source", mcposasource, mcrotasource)
  mccomp_posa[2] = mcposasource;
  mccomp_posr[2] = mcposrsource;
  mcNCounter[2]  = mcPCounter[2] = mcP2Counter[2] = 0;
  mcAbsorbProp[2]= 0;
    /* Component source_div. */
  /* Setting parameters for component source_div. */
  SIG_MESSAGE("source_div (Init:SetPar)");
#line 67 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccsource_div_xwidth = mcipsource_width;
#line 66 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccsource_div_yheight = mcipsource_width;
#line 69 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccsource_div_focus_aw = mcipsource_divergence;
#line 70 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccsource_div_focus_ah = mcipsource_divergence;
#line 64 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccsource_div_E0 = 0.0;
#line 64 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccsource_div_dE = 0.0;
#line 71 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccsource_div_lambda0 = mcipL_source;
#line 73 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccsource_div_dlambda = mcipdL;
#line 64 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccsource_div_gauss = 0;
#line 72 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccsource_div_flux = mcipflux;
#line 10405 "./ess_extraction.c"

  SIG_MESSAGE("source_div (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 10412 "./ess_extraction.c"
  rot_mul(mctr1, mcrotasource, mcrotasource_div);
  rot_transpose(mcrotasource, mctr1);
  rot_mul(mcrotasource_div, mctr1, mcrotrsource_div);
  mctc1 = coords_set(
#line 74 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
    0,
#line 74 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
    0,
#line 74 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
    0);
#line 10423 "./ess_extraction.c"
  rot_transpose(mcrotasource, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposasource_div = coords_add(mcposasource, mctc2);
  mctc1 = coords_sub(mcposasource, mcposasource_div);
  mcposrsource_div = rot_apply(mcrotasource_div, mctc1);
  mcDEBUG_COMPONENT("source_div", mcposasource_div, mcrotasource_div)
  mccomp_posa[3] = mcposasource_div;
  mccomp_posr[3] = mcposrsource_div;
  mcNCounter[3]  = mcPCounter[3] = mcP2Counter[3] = 0;
  mcAbsorbProp[3]= 0;
    /* Component psd_monitor_source. */
  /* Setting parameters for component psd_monitor_source. */
  SIG_MESSAGE("psd_monitor_source (Init:SetPar)");
#line 80 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccpsd_monitor_source_nx = mcippixels;
#line 81 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccpsd_monitor_source_ny = mcippixels;
#line 77 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  if("source_psd.dat") strncpy(mccpsd_monitor_source_filename, "source_psd.dat" ? "source_psd.dat" : "", 16384); else mccpsd_monitor_source_filename[0]='\0';
#line 50 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccpsd_monitor_source_xmin = -0.05;
#line 50 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccpsd_monitor_source_xmax = 0.05;
#line 50 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccpsd_monitor_source_ymin = -0.05;
#line 50 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccpsd_monitor_source_ymax = 0.05;
#line 78 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccpsd_monitor_source_xwidth = mcipdet_width;
#line 79 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccpsd_monitor_source_yheight = mcipdet_width;
#line 82 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccpsd_monitor_source_restore_neutron = 1;
#line 10457 "./ess_extraction.c"

  SIG_MESSAGE("psd_monitor_source (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 10464 "./ess_extraction.c"
  rot_mul(mctr1, mcrotasource, mcrotapsd_monitor_source);
  rot_transpose(mcrotasource_div, mctr1);
  rot_mul(mcrotapsd_monitor_source, mctr1, mcrotrpsd_monitor_source);
  mctc1 = coords_set(
#line 83 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
    0,
#line 83 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
    0,
#line 83 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
    mcipfocal_length -1);
#line 10475 "./ess_extraction.c"
  rot_transpose(mcrotasource, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposapsd_monitor_source = coords_add(mcposasource, mctc2);
  mctc1 = coords_sub(mcposasource_div, mcposapsd_monitor_source);
  mcposrpsd_monitor_source = rot_apply(mcrotapsd_monitor_source, mctc1);
  mcDEBUG_COMPONENT("psd_monitor_source", mcposapsd_monitor_source, mcrotapsd_monitor_source)
  mccomp_posa[4] = mcposapsd_monitor_source;
  mccomp_posr[4] = mcposrpsd_monitor_source;
  mcNCounter[4]  = mcPCounter[4] = mcP2Counter[4] = 0;
  mcAbsorbProp[4]= 0;
    /* Component parabolic_optic_before_guide_v. */
  /* Setting parameters for component parabolic_optic_before_guide_v. */
  SIG_MESSAGE("parabolic_optic_before_guide_v (Init:SetPar)");
#line 88 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccparabolic_optic_before_guide_v_sourceDist = - mcipfocal_length + 1;
#line 89 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccparabolic_optic_before_guide_v_LStart = - mcipfocal_length + 1;
#line 90 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccparabolic_optic_before_guide_v_LEnd = mcipincoming_length + 1;
#line 91 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccparabolic_optic_before_guide_v_lStart = 0;
#line 92 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccparabolic_optic_before_guide_v_lEnd = 1;
#line 93 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccparabolic_optic_before_guide_v_r_0 = 0.0457293;
#line 95 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccparabolic_optic_before_guide_v_nummirror = mcipmirrors;
#line 97 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccparabolic_optic_before_guide_v_mf = 4.1;
#line 98 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccparabolic_optic_before_guide_v_mb = 0;
#line 94 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccparabolic_optic_before_guide_v_mirror_width = 0;
#line 96 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccparabolic_optic_before_guide_v_doubleReflections = 1;
#line 10511 "./ess_extraction.c"

  SIG_MESSAGE("parabolic_optic_before_guide_v (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 102 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
    (0)*DEG2RAD,
#line 102 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
    (0)*DEG2RAD,
#line 102 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
    (90)*DEG2RAD);
#line 10521 "./ess_extraction.c"
  rot_mul(mctr1, mcrotasource, mcrotaparabolic_optic_before_guide_v);
  rot_transpose(mcrotapsd_monitor_source, mctr1);
  rot_mul(mcrotaparabolic_optic_before_guide_v, mctr1, mcrotrparabolic_optic_before_guide_v);
  mctc1 = coords_set(
#line 101 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
    0,
#line 101 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
    0,
#line 101 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
    mcipfocal_length -1);
#line 10532 "./ess_extraction.c"
  rot_transpose(mcrotasource, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaparabolic_optic_before_guide_v = coords_add(mcposasource, mctc2);
  mctc1 = coords_sub(mcposapsd_monitor_source, mcposaparabolic_optic_before_guide_v);
  mcposrparabolic_optic_before_guide_v = rot_apply(mcrotaparabolic_optic_before_guide_v, mctc1);
  mcDEBUG_COMPONENT("parabolic_optic_before_guide_v", mcposaparabolic_optic_before_guide_v, mcrotaparabolic_optic_before_guide_v)
  mccomp_posa[5] = mcposaparabolic_optic_before_guide_v;
  mccomp_posr[5] = mcposrparabolic_optic_before_guide_v;
  mcNCounter[5]  = mcPCounter[5] = mcP2Counter[5] = 0;
  mcAbsorbProp[5]= 0;
    /* Component parabolic_optic_before_guide_h. */
  /* Setting parameters for component parabolic_optic_before_guide_h. */
  SIG_MESSAGE("parabolic_optic_before_guide_h (Init:SetPar)");
#line 106 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccparabolic_optic_before_guide_h_sourceDist = - mcipfocal_length;
#line 107 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccparabolic_optic_before_guide_h_LStart = - mcipfocal_length;
#line 108 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccparabolic_optic_before_guide_h_LEnd = mcipincoming_length;
#line 109 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccparabolic_optic_before_guide_h_lStart = 0;
#line 110 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccparabolic_optic_before_guide_h_lEnd = 1;
#line 111 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccparabolic_optic_before_guide_h_r_0 = 0.048502;
#line 113 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccparabolic_optic_before_guide_h_nummirror = mcipmirrors;
#line 115 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccparabolic_optic_before_guide_h_mf = 4.1;
#line 116 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccparabolic_optic_before_guide_h_mb = 0;
#line 112 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccparabolic_optic_before_guide_h_mirror_width = 0;
#line 114 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccparabolic_optic_before_guide_h_doubleReflections = 1;
#line 10568 "./ess_extraction.c"

  SIG_MESSAGE("parabolic_optic_before_guide_h (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 120 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
    (0)*DEG2RAD,
#line 120 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
    (0)*DEG2RAD,
#line 120 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
    (0)*DEG2RAD);
#line 10578 "./ess_extraction.c"
  rot_mul(mctr1, mcrotasource, mcrotaparabolic_optic_before_guide_h);
  rot_transpose(mcrotaparabolic_optic_before_guide_v, mctr1);
  rot_mul(mcrotaparabolic_optic_before_guide_h, mctr1, mcrotrparabolic_optic_before_guide_h);
  mctc1 = coords_set(
#line 119 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
    0,
#line 119 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
    0,
#line 119 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
    mcipfocal_length);
#line 10589 "./ess_extraction.c"
  rot_transpose(mcrotasource, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaparabolic_optic_before_guide_h = coords_add(mcposasource, mctc2);
  mctc1 = coords_sub(mcposaparabolic_optic_before_guide_v, mcposaparabolic_optic_before_guide_h);
  mcposrparabolic_optic_before_guide_h = rot_apply(mcrotaparabolic_optic_before_guide_h, mctc1);
  mcDEBUG_COMPONENT("parabolic_optic_before_guide_h", mcposaparabolic_optic_before_guide_h, mcrotaparabolic_optic_before_guide_h)
  mccomp_posa[6] = mcposaparabolic_optic_before_guide_h;
  mccomp_posr[6] = mcposrparabolic_optic_before_guide_h;
  mcNCounter[6]  = mcPCounter[6] = mcP2Counter[6] = 0;
  mcAbsorbProp[6]= 0;
    /* Component after_optic_source. */
  /* Setting parameters for component after_optic_source. */
  SIG_MESSAGE("after_optic_source (Init:SetPar)");

  SIG_MESSAGE("after_optic_source (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 10609 "./ess_extraction.c"
  rot_mul(mctr1, mcrotasource, mcrotaafter_optic_source);
  rot_transpose(mcrotaparabolic_optic_before_guide_h, mctr1);
  rot_mul(mcrotaafter_optic_source, mctr1, mcrotrafter_optic_source);
  mctc1 = coords_set(
#line 125 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
    0,
#line 125 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
    0,
#line 125 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
    mcipfocal_length + 1);
#line 10620 "./ess_extraction.c"
  rot_transpose(mcrotasource, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaafter_optic_source = coords_add(mcposasource, mctc2);
  mctc1 = coords_sub(mcposaparabolic_optic_before_guide_h, mcposaafter_optic_source);
  mcposrafter_optic_source = rot_apply(mcrotaafter_optic_source, mctc1);
  mcDEBUG_COMPONENT("after_optic_source", mcposaafter_optic_source, mcrotaafter_optic_source)
  mccomp_posa[7] = mcposaafter_optic_source;
  mccomp_posr[7] = mcposrafter_optic_source;
  mcNCounter[7]  = mcPCounter[7] = mcP2Counter[7] = 0;
  mcAbsorbProp[7]= 0;
    /* Component psd_monitor_afteropticsource. */
  /* Setting parameters for component psd_monitor_afteropticsource. */
  SIG_MESSAGE("psd_monitor_afteropticsource (Init:SetPar)");
#line 131 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccpsd_monitor_afteropticsource_nx = mcippixels;
#line 132 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccpsd_monitor_afteropticsource_ny = mcippixels;
#line 128 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  if("sourceafteroptic_psd.dat") strncpy(mccpsd_monitor_afteropticsource_filename, "sourceafteroptic_psd.dat" ? "sourceafteroptic_psd.dat" : "", 16384); else mccpsd_monitor_afteropticsource_filename[0]='\0';
#line 50 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccpsd_monitor_afteropticsource_xmin = -0.05;
#line 50 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccpsd_monitor_afteropticsource_xmax = 0.05;
#line 50 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccpsd_monitor_afteropticsource_ymin = -0.05;
#line 50 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccpsd_monitor_afteropticsource_ymax = 0.05;
#line 129 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccpsd_monitor_afteropticsource_xwidth = mcipdet_width;
#line 130 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccpsd_monitor_afteropticsource_yheight = mcipdet_width;
#line 133 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccpsd_monitor_afteropticsource_restore_neutron = 1;
#line 10654 "./ess_extraction.c"

  SIG_MESSAGE("psd_monitor_afteropticsource (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 10661 "./ess_extraction.c"
  rot_mul(mctr1, mcrotaafter_optic_source, mcrotapsd_monitor_afteropticsource);
  rot_transpose(mcrotaafter_optic_source, mctr1);
  rot_mul(mcrotapsd_monitor_afteropticsource, mctr1, mcrotrpsd_monitor_afteropticsource);
  mctc1 = coords_set(
#line 134 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
    0,
#line 134 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
    0,
#line 134 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
    0);
#line 10672 "./ess_extraction.c"
  rot_transpose(mcrotaafter_optic_source, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposapsd_monitor_afteropticsource = coords_add(mcposaafter_optic_source, mctc2);
  mctc1 = coords_sub(mcposaafter_optic_source, mcposapsd_monitor_afteropticsource);
  mcposrpsd_monitor_afteropticsource = rot_apply(mcrotapsd_monitor_afteropticsource, mctc1);
  mcDEBUG_COMPONENT("psd_monitor_afteropticsource", mcposapsd_monitor_afteropticsource, mcrotapsd_monitor_afteropticsource)
  mccomp_posa[8] = mcposapsd_monitor_afteropticsource;
  mccomp_posr[8] = mcposrpsd_monitor_afteropticsource;
  mcNCounter[8]  = mcPCounter[8] = mcP2Counter[8] = 0;
  mcAbsorbProp[8]= 0;
    /* Component guide_gravity_1. */
  /* Setting parameters for component guide_gravity_1. */
  SIG_MESSAGE("guide_gravity_1 (Init:SetPar)");
#line 139 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccguide_gravity_1_w1 = mcipguide_width;
#line 140 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccguide_gravity_1_h1 = mcipguide_width;
#line 113 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccguide_gravity_1_w2 = 0;
#line 113 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccguide_gravity_1_h2 = 0;
#line 143 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccguide_gravity_1_l = mcipguide_length;
#line 114 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccguide_gravity_1_R0 = 0.995;
#line 114 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccguide_gravity_1_Qc = 0.0218;
#line 114 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccguide_gravity_1_alpha = 4.38;
#line 142 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccguide_gravity_1_m = 2;
#line 114 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccguide_gravity_1_W = 0.003;
#line 114 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccguide_gravity_1_nslit = 1;
#line 114 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccguide_gravity_1_d = 0.0005;
#line 115 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccguide_gravity_1_mleft = -1;
#line 115 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccguide_gravity_1_mright = -1;
#line 115 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccguide_gravity_1_mtop = -1;
#line 115 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccguide_gravity_1_mbottom = -1;
#line 115 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccguide_gravity_1_nhslit = 1;
#line 141 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccguide_gravity_1_G = mcipg;
#line 116 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccguide_gravity_1_aleft = -1;
#line 116 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccguide_gravity_1_aright = -1;
#line 116 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccguide_gravity_1_atop = -1;
#line 116 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccguide_gravity_1_abottom = -1;
#line 117 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccguide_gravity_1_wavy = 0;
#line 117 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccguide_gravity_1_wavy_z = 0;
#line 117 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccguide_gravity_1_wavy_tb = 0;
#line 117 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccguide_gravity_1_wavy_lr = 0;
#line 118 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccguide_gravity_1_chamfers = 0;
#line 118 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccguide_gravity_1_chamfers_z = 0;
#line 118 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccguide_gravity_1_chamfers_lr = 0;
#line 118 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccguide_gravity_1_chamfers_tb = 0;
#line 118 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccguide_gravity_1_nelements = 1;
#line 119 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccguide_gravity_1_nu = 0;
#line 119 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccguide_gravity_1_phase = 0;
#line 119 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  if("NULL") strncpy(mccguide_gravity_1_reflect, "NULL" ? "NULL" : "", 16384); else mccguide_gravity_1_reflect[0]='\0';
#line 10754 "./ess_extraction.c"

  SIG_MESSAGE("guide_gravity_1 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 10761 "./ess_extraction.c"
  rot_mul(mctr1, mcrotaafter_optic_source, mcrotaguide_gravity_1);
  rot_transpose(mcrotapsd_monitor_afteropticsource, mctr1);
  rot_mul(mcrotaguide_gravity_1, mctr1, mcrotrguide_gravity_1);
  mctc1 = coords_set(
#line 144 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
    0,
#line 144 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
    0,
#line 144 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
    0);
#line 10772 "./ess_extraction.c"
  rot_transpose(mcrotaafter_optic_source, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaguide_gravity_1 = coords_add(mcposaafter_optic_source, mctc2);
  mctc1 = coords_sub(mcposapsd_monitor_afteropticsource, mcposaguide_gravity_1);
  mcposrguide_gravity_1 = rot_apply(mcrotaguide_gravity_1, mctc1);
  mcDEBUG_COMPONENT("guide_gravity_1", mcposaguide_gravity_1, mcrotaguide_gravity_1)
  mccomp_posa[9] = mcposaguide_gravity_1;
  mccomp_posr[9] = mcposrguide_gravity_1;
  mcNCounter[9]  = mcPCounter[9] = mcP2Counter[9] = 0;
  mcAbsorbProp[9]= 0;
    /* Component monitor_1. */
  /* Setting parameters for component monitor_1. */
  SIG_MESSAGE("monitor_1 (Init:SetPar)");

  SIG_MESSAGE("monitor_1 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 10792 "./ess_extraction.c"
  rot_mul(mctr1, mcrotaafter_optic_source, mcrotamonitor_1);
  rot_transpose(mcrotaguide_gravity_1, mctr1);
  rot_mul(mcrotamonitor_1, mctr1, mcrotrmonitor_1);
  mctc1 = coords_set(
#line 147 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
    0,
#line 147 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
    0,
#line 147 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
    mcipguide_length + 0.0001);
#line 10803 "./ess_extraction.c"
  rot_transpose(mcrotaafter_optic_source, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposamonitor_1 = coords_add(mcposaafter_optic_source, mctc2);
  mctc1 = coords_sub(mcposaguide_gravity_1, mcposamonitor_1);
  mcposrmonitor_1 = rot_apply(mcrotamonitor_1, mctc1);
  mcDEBUG_COMPONENT("monitor_1", mcposamonitor_1, mcrotamonitor_1)
  mccomp_posa[10] = mcposamonitor_1;
  mccomp_posr[10] = mcposrmonitor_1;
  mcNCounter[10]  = mcPCounter[10] = mcP2Counter[10] = 0;
  mcAbsorbProp[10]= 0;
    /* Component divpos_monitor. */
  /* Setting parameters for component divpos_monitor. */
  SIG_MESSAGE("divpos_monitor (Init:SetPar)");
#line 152 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  if("g1_divpos.dat") strncpy(mccdivpos_monitor_filename, "g1_divpos.dat" ? "g1_divpos.dat" : "", 16384); else mccdivpos_monitor_filename[0]='\0';
#line 55 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccdivpos_monitor_xmin = -0.05;
#line 55 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccdivpos_monitor_xmax = 0.05;
#line 55 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccdivpos_monitor_ymin = -0.05;
#line 55 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccdivpos_monitor_ymax = 0.05;
#line 153 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccdivpos_monitor_xwidth = mcipdet_width;
#line 154 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccdivpos_monitor_yheight = mcipdet_width;
#line 155 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccdivpos_monitor_maxdiv_h = mcipdivergence_max * 2;
#line 156 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccdivpos_monitor_restore_neutron = 1;
#line 56 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccdivpos_monitor_nx = 0;
#line 56 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccdivpos_monitor_ny = 0;
#line 56 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccdivpos_monitor_nz = 1;
#line 56 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccdivpos_monitor_nowritefile = 0;
#line 10843 "./ess_extraction.c"

  SIG_MESSAGE("divpos_monitor (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 158 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
    (0)*DEG2RAD,
#line 158 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
    (0)*DEG2RAD,
#line 158 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
    (90)*DEG2RAD);
#line 10853 "./ess_extraction.c"
  rot_mul(mctr1, mcrotamonitor_1, mcrotadivpos_monitor);
  rot_transpose(mcrotamonitor_1, mctr1);
  rot_mul(mcrotadivpos_monitor, mctr1, mcrotrdivpos_monitor);
  mctc1 = coords_set(
#line 157 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
    0,
#line 157 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
    0,
#line 157 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
    0);
#line 10864 "./ess_extraction.c"
  rot_transpose(mcrotamonitor_1, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposadivpos_monitor = coords_add(mcposamonitor_1, mctc2);
  mctc1 = coords_sub(mcposamonitor_1, mcposadivpos_monitor);
  mcposrdivpos_monitor = rot_apply(mcrotadivpos_monitor, mctc1);
  mcDEBUG_COMPONENT("divpos_monitor", mcposadivpos_monitor, mcrotadivpos_monitor)
  mccomp_posa[11] = mcposadivpos_monitor;
  mccomp_posr[11] = mcposrdivpos_monitor;
  mcNCounter[11]  = mcPCounter[11] = mcP2Counter[11] = 0;
  mcAbsorbProp[11]= 0;
    /* Component psd_monitor_g1. */
  /* Setting parameters for component psd_monitor_g1. */
  SIG_MESSAGE("psd_monitor_g1 (Init:SetPar)");
#line 221 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccpsd_monitor_g1_nx = mcippixels;
#line 222 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccpsd_monitor_g1_ny = mcippixels;
#line 218 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  if("g1_psd.dat") strncpy(mccpsd_monitor_g1_filename, "g1_psd.dat" ? "g1_psd.dat" : "", 16384); else mccpsd_monitor_g1_filename[0]='\0';
#line 50 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccpsd_monitor_g1_xmin = -0.05;
#line 50 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccpsd_monitor_g1_xmax = 0.05;
#line 50 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccpsd_monitor_g1_ymin = -0.05;
#line 50 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccpsd_monitor_g1_ymax = 0.05;
#line 219 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccpsd_monitor_g1_xwidth = mcipdet_width;
#line 220 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccpsd_monitor_g1_yheight = mcipdet_width;
#line 223 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccpsd_monitor_g1_restore_neutron = 1;
#line 10898 "./ess_extraction.c"

  SIG_MESSAGE("psd_monitor_g1 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 10905 "./ess_extraction.c"
  rot_mul(mctr1, mcrotamonitor_1, mcrotapsd_monitor_g1);
  rot_transpose(mcrotadivpos_monitor, mctr1);
  rot_mul(mcrotapsd_monitor_g1, mctr1, mcrotrpsd_monitor_g1);
  mctc1 = coords_set(
#line 224 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
    0,
#line 224 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
    0,
#line 224 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
    0);
#line 10916 "./ess_extraction.c"
  rot_transpose(mcrotamonitor_1, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposapsd_monitor_g1 = coords_add(mcposamonitor_1, mctc2);
  mctc1 = coords_sub(mcposadivpos_monitor, mcposapsd_monitor_g1);
  mcposrpsd_monitor_g1 = rot_apply(mcrotapsd_monitor_g1, mctc1);
  mcDEBUG_COMPONENT("psd_monitor_g1", mcposapsd_monitor_g1, mcrotapsd_monitor_g1)
  mccomp_posa[12] = mcposapsd_monitor_g1;
  mccomp_posr[12] = mcposrpsd_monitor_g1;
  mcNCounter[12]  = mcPCounter[12] = mcP2Counter[12] = 0;
  mcAbsorbProp[12]= 0;
    /* Component divhlambda_monitor_g1. */
  /* Setting parameters for component divhlambda_monitor_g1. */
  SIG_MESSAGE("divhlambda_monitor_g1 (Init:SetPar)");
#line 229 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  if("g1_divv_lambda.dat") strncpy(mccdivhlambda_monitor_g1_filename, "g1_divv_lambda.dat" ? "g1_divv_lambda.dat" : "", 16384); else mccdivhlambda_monitor_g1_filename[0]='\0';
#line 57 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccdivhlambda_monitor_g1_xmin = -0.05;
#line 57 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccdivhlambda_monitor_g1_xmax = 0.05;
#line 57 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccdivhlambda_monitor_g1_ymin = -0.05;
#line 57 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccdivhlambda_monitor_g1_ymax = 0.05;
#line 230 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccdivhlambda_monitor_g1_xwidth = mcipdet_width;
#line 231 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccdivhlambda_monitor_g1_yheight = mcipdet_width;
#line 232 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccdivhlambda_monitor_g1_maxdiv_h = mcipdivergence_max * 2;
#line 233 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccdivhlambda_monitor_g1_Lmin = mcipL_min;
#line 234 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccdivhlambda_monitor_g1_Lmax = mcipL_max;
#line 235 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccdivhlambda_monitor_g1_restore_neutron = 1;
#line 57 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccdivhlambda_monitor_g1_nx = 0;
#line 57 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccdivhlambda_monitor_g1_ny = 0;
#line 57 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccdivhlambda_monitor_g1_nz = 1;
#line 57 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccdivhlambda_monitor_g1_nowritefile = 0;
#line 10960 "./ess_extraction.c"

  SIG_MESSAGE("divhlambda_monitor_g1 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 237 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
    (0)*DEG2RAD,
#line 237 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
    (0)*DEG2RAD,
#line 237 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
    (90)*DEG2RAD);
#line 10970 "./ess_extraction.c"
  rot_mul(mctr1, mcrotamonitor_1, mcrotadivhlambda_monitor_g1);
  rot_transpose(mcrotapsd_monitor_g1, mctr1);
  rot_mul(mcrotadivhlambda_monitor_g1, mctr1, mcrotrdivhlambda_monitor_g1);
  mctc1 = coords_set(
#line 236 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
    0,
#line 236 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
    0,
#line 236 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
    0);
#line 10981 "./ess_extraction.c"
  rot_transpose(mcrotamonitor_1, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposadivhlambda_monitor_g1 = coords_add(mcposamonitor_1, mctc2);
  mctc1 = coords_sub(mcposapsd_monitor_g1, mcposadivhlambda_monitor_g1);
  mcposrdivhlambda_monitor_g1 = rot_apply(mcrotadivhlambda_monitor_g1, mctc1);
  mcDEBUG_COMPONENT("divhlambda_monitor_g1", mcposadivhlambda_monitor_g1, mcrotadivhlambda_monitor_g1)
  mccomp_posa[13] = mcposadivhlambda_monitor_g1;
  mccomp_posr[13] = mcposrdivhlambda_monitor_g1;
  mcNCounter[13]  = mcPCounter[13] = mcP2Counter[13] = 0;
  mcAbsorbProp[13]= 0;
    /* Component parabolic_optic1. */
  /* Setting parameters for component parabolic_optic1. */
  SIG_MESSAGE("parabolic_optic1 (Init:SetPar)");
#line 255 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccparabolic_optic1_sourceDist = - mcipincoming_length;
#line 256 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccparabolic_optic1_LStart = - mcipincoming_length;
#line 257 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccparabolic_optic1_LEnd = mcipfocal_length + 1;
#line 258 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccparabolic_optic1_lStart = 0;
#line 259 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccparabolic_optic1_lEnd = 1;
#line 260 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccparabolic_optic1_r_0 = mcipguide_width / 2 + 0.001 + 0.005;
#line 262 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccparabolic_optic1_nummirror = mcipmirrors;
#line 264 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccparabolic_optic1_mf = 4.1;
#line 265 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccparabolic_optic1_mb = 0;
#line 261 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccparabolic_optic1_mirror_width = 0;
#line 263 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccparabolic_optic1_doubleReflections = 1;
#line 11017 "./ess_extraction.c"

  SIG_MESSAGE("parabolic_optic1 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 269 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
    (0)*DEG2RAD,
#line 269 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
    (0)*DEG2RAD,
#line 269 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
    (90)*DEG2RAD);
#line 11027 "./ess_extraction.c"
  rot_mul(mctr1, mcrotamonitor_1, mcrotaparabolic_optic1);
  rot_transpose(mcrotadivhlambda_monitor_g1, mctr1);
  rot_mul(mcrotaparabolic_optic1, mctr1, mcrotrparabolic_optic1);
  mctc1 = coords_set(
#line 268 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
    0,
#line 268 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
    0,
#line 268 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
    0);
#line 11038 "./ess_extraction.c"
  rot_transpose(mcrotamonitor_1, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaparabolic_optic1 = coords_add(mcposamonitor_1, mctc2);
  mctc1 = coords_sub(mcposadivhlambda_monitor_g1, mcposaparabolic_optic1);
  mcposrparabolic_optic1 = rot_apply(mcrotaparabolic_optic1, mctc1);
  mcDEBUG_COMPONENT("parabolic_optic1", mcposaparabolic_optic1, mcrotaparabolic_optic1)
  mccomp_posa[14] = mcposaparabolic_optic1;
  mccomp_posr[14] = mcposrparabolic_optic1;
  mcNCounter[14]  = mcPCounter[14] = mcP2Counter[14] = 0;
  mcAbsorbProp[14]= 0;
    /* Component parabolic_optic2. */
  /* Setting parameters for component parabolic_optic2. */
  SIG_MESSAGE("parabolic_optic2 (Init:SetPar)");
#line 272 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccparabolic_optic2_sourceDist = - mcipincoming_length -1;
#line 273 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccparabolic_optic2_LStart = - mcipincoming_length -1;
#line 274 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccparabolic_optic2_LEnd = mcipfocal_length + 0;
#line 275 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccparabolic_optic2_lStart = 0;
#line 276 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccparabolic_optic2_lEnd = 1;
#line 277 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccparabolic_optic2_r_0 = mcipguide_width / 2 + 0.001 + 0.005;
#line 279 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccparabolic_optic2_nummirror = mcipmirrors;
#line 281 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccparabolic_optic2_mf = 4.1;
#line 282 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccparabolic_optic2_mb = 0;
#line 278 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccparabolic_optic2_mirror_width = 0;
#line 280 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccparabolic_optic2_doubleReflections = 1;
#line 11074 "./ess_extraction.c"

  SIG_MESSAGE("parabolic_optic2 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 286 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
    (0)*DEG2RAD,
#line 286 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
    (0)*DEG2RAD,
#line 286 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
    (0)*DEG2RAD);
#line 11084 "./ess_extraction.c"
  rot_mul(mctr1, mcrotamonitor_1, mcrotaparabolic_optic2);
  rot_transpose(mcrotaparabolic_optic1, mctr1);
  rot_mul(mcrotaparabolic_optic2, mctr1, mcrotrparabolic_optic2);
  mctc1 = coords_set(
#line 285 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
    0,
#line 285 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
    0,
#line 285 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
    1);
#line 11095 "./ess_extraction.c"
  rot_transpose(mcrotamonitor_1, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaparabolic_optic2 = coords_add(mcposamonitor_1, mctc2);
  mctc1 = coords_sub(mcposaparabolic_optic1, mcposaparabolic_optic2);
  mcposrparabolic_optic2 = rot_apply(mcrotaparabolic_optic2, mctc1);
  mcDEBUG_COMPONENT("parabolic_optic2", mcposaparabolic_optic2, mcrotaparabolic_optic2)
  mccomp_posa[15] = mcposaparabolic_optic2;
  mccomp_posr[15] = mcposrparabolic_optic2;
  mcNCounter[15]  = mcPCounter[15] = mcP2Counter[15] = 0;
  mcAbsorbProp[15]= 0;
    /* Component monitor_2. */
  /* Setting parameters for component monitor_2. */
  SIG_MESSAGE("monitor_2 (Init:SetPar)");

  SIG_MESSAGE("monitor_2 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 293 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
    (0)*DEG2RAD,
#line 293 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
    (0)*DEG2RAD,
#line 293 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
    (0)*DEG2RAD);
#line 11118 "./ess_extraction.c"
  rot_mul(mctr1, mcrotamonitor_1, mcrotamonitor_2);
  rot_transpose(mcrotaparabolic_optic2, mctr1);
  rot_mul(mcrotamonitor_2, mctr1, mcrotrmonitor_2);
  mctc1 = coords_set(
#line 292 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
    0,
#line 292 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
    0,
#line 292 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
    mcipfocal_length + 1);
#line 11129 "./ess_extraction.c"
  rot_transpose(mcrotamonitor_1, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposamonitor_2 = coords_add(mcposamonitor_1, mctc2);
  mctc1 = coords_sub(mcposaparabolic_optic2, mcposamonitor_2);
  mcposrmonitor_2 = rot_apply(mcrotamonitor_2, mctc1);
  mcDEBUG_COMPONENT("monitor_2", mcposamonitor_2, mcrotamonitor_2)
  mccomp_posa[16] = mcposamonitor_2;
  mccomp_posr[16] = mcposrmonitor_2;
  mcNCounter[16]  = mcPCounter[16] = mcP2Counter[16] = 0;
  mcAbsorbProp[16]= 0;
    /* Component psd_monitor_f_zoom. */
  /* Setting parameters for component psd_monitor_f_zoom. */
  SIG_MESSAGE("psd_monitor_f_zoom (Init:SetPar)");
#line 355 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccpsd_monitor_f_zoom_nx = mcippixels;
#line 356 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccpsd_monitor_f_zoom_ny = mcippixels;
#line 352 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  if("f_psd_zoom.dat") strncpy(mccpsd_monitor_f_zoom_filename, "f_psd_zoom.dat" ? "f_psd_zoom.dat" : "", 16384); else mccpsd_monitor_f_zoom_filename[0]='\0';
#line 50 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccpsd_monitor_f_zoom_xmin = -0.05;
#line 50 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccpsd_monitor_f_zoom_xmax = 0.05;
#line 50 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccpsd_monitor_f_zoom_ymin = -0.05;
#line 50 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccpsd_monitor_f_zoom_ymax = 0.05;
#line 353 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccpsd_monitor_f_zoom_xwidth = mcipdet_width_focus;
#line 354 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccpsd_monitor_f_zoom_yheight = mcipdet_width_focus;
#line 357 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccpsd_monitor_f_zoom_restore_neutron = 1;
#line 11163 "./ess_extraction.c"

  SIG_MESSAGE("psd_monitor_f_zoom (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 359 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
    (0)*DEG2RAD,
#line 359 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
    (0)*DEG2RAD,
#line 359 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
    (0)*DEG2RAD);
#line 11173 "./ess_extraction.c"
  rot_mul(mctr1, mcrotamonitor_2, mcrotapsd_monitor_f_zoom);
  rot_transpose(mcrotamonitor_2, mctr1);
  rot_mul(mcrotapsd_monitor_f_zoom, mctr1, mcrotrpsd_monitor_f_zoom);
  mctc1 = coords_set(
#line 358 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
    0,
#line 358 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
    0,
#line 358 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
    0);
#line 11184 "./ess_extraction.c"
  rot_transpose(mcrotamonitor_2, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposapsd_monitor_f_zoom = coords_add(mcposamonitor_2, mctc2);
  mctc1 = coords_sub(mcposamonitor_2, mcposapsd_monitor_f_zoom);
  mcposrpsd_monitor_f_zoom = rot_apply(mcrotapsd_monitor_f_zoom, mctc1);
  mcDEBUG_COMPONENT("psd_monitor_f_zoom", mcposapsd_monitor_f_zoom, mcrotapsd_monitor_f_zoom)
  mccomp_posa[17] = mcposapsd_monitor_f_zoom;
  mccomp_posr[17] = mcposrpsd_monitor_f_zoom;
  mcNCounter[17]  = mcPCounter[17] = mcP2Counter[17] = 0;
  mcAbsorbProp[17]= 0;
    /* Component psd_monitor_f. */
  /* Setting parameters for component psd_monitor_f. */
  SIG_MESSAGE("psd_monitor_f (Init:SetPar)");
#line 365 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccpsd_monitor_f_nx = mcippixels;
#line 366 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccpsd_monitor_f_ny = mcippixels;
#line 362 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  if("f_psd.dat") strncpy(mccpsd_monitor_f_filename, "f_psd.dat" ? "f_psd.dat" : "", 16384); else mccpsd_monitor_f_filename[0]='\0';
#line 50 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccpsd_monitor_f_xmin = -0.05;
#line 50 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccpsd_monitor_f_xmax = 0.05;
#line 50 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccpsd_monitor_f_ymin = -0.05;
#line 50 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccpsd_monitor_f_ymax = 0.05;
#line 363 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccpsd_monitor_f_xwidth = mcipdet_width;
#line 364 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccpsd_monitor_f_yheight = mcipdet_width;
#line 367 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccpsd_monitor_f_restore_neutron = 1;
#line 11218 "./ess_extraction.c"

  SIG_MESSAGE("psd_monitor_f (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 369 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
    (0)*DEG2RAD,
#line 369 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
    (0)*DEG2RAD,
#line 369 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
    (0)*DEG2RAD);
#line 11228 "./ess_extraction.c"
  rot_mul(mctr1, mcrotamonitor_2, mcrotapsd_monitor_f);
  rot_transpose(mcrotapsd_monitor_f_zoom, mctr1);
  rot_mul(mcrotapsd_monitor_f, mctr1, mcrotrpsd_monitor_f);
  mctc1 = coords_set(
#line 368 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
    0,
#line 368 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
    0,
#line 368 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
    0);
#line 11239 "./ess_extraction.c"
  rot_transpose(mcrotamonitor_2, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposapsd_monitor_f = coords_add(mcposamonitor_2, mctc2);
  mctc1 = coords_sub(mcposapsd_monitor_f_zoom, mcposapsd_monitor_f);
  mcposrpsd_monitor_f = rot_apply(mcrotapsd_monitor_f, mctc1);
  mcDEBUG_COMPONENT("psd_monitor_f", mcposapsd_monitor_f, mcrotapsd_monitor_f)
  mccomp_posa[18] = mcposapsd_monitor_f;
  mccomp_posr[18] = mcposrpsd_monitor_f;
  mcNCounter[18]  = mcPCounter[18] = mcP2Counter[18] = 0;
  mcAbsorbProp[18]= 0;
    /* Component f_divpos. */
  /* Setting parameters for component f_divpos. */
  SIG_MESSAGE("f_divpos (Init:SetPar)");
#line 374 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  if("f_divpos.dat") strncpy(mccf_divpos_filename, "f_divpos.dat" ? "f_divpos.dat" : "", 16384); else mccf_divpos_filename[0]='\0';
#line 55 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccf_divpos_xmin = -0.05;
#line 55 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccf_divpos_xmax = 0.05;
#line 55 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccf_divpos_ymin = -0.05;
#line 55 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccf_divpos_ymax = 0.05;
#line 375 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccf_divpos_xwidth = mcipdet_width_focus;
#line 376 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccf_divpos_yheight = mcipdet_width_focus;
#line 377 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccf_divpos_maxdiv_h = mcipdivergence_max * 2;
#line 378 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccf_divpos_restore_neutron = 1;
#line 56 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccf_divpos_nx = 0;
#line 56 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccf_divpos_ny = 0;
#line 56 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccf_divpos_nz = 1;
#line 56 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccf_divpos_nowritefile = 0;
#line 11279 "./ess_extraction.c"

  SIG_MESSAGE("f_divpos (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 380 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
    (0)*DEG2RAD,
#line 380 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
    (0)*DEG2RAD,
#line 380 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
    (90)*DEG2RAD);
#line 11289 "./ess_extraction.c"
  rot_mul(mctr1, mcrotamonitor_2, mcrotaf_divpos);
  rot_transpose(mcrotapsd_monitor_f, mctr1);
  rot_mul(mcrotaf_divpos, mctr1, mcrotrf_divpos);
  mctc1 = coords_set(
#line 379 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
    0,
#line 379 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
    0,
#line 379 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
    0);
#line 11300 "./ess_extraction.c"
  rot_transpose(mcrotamonitor_2, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaf_divpos = coords_add(mcposamonitor_2, mctc2);
  mctc1 = coords_sub(mcposapsd_monitor_f, mcposaf_divpos);
  mcposrf_divpos = rot_apply(mcrotaf_divpos, mctc1);
  mcDEBUG_COMPONENT("f_divpos", mcposaf_divpos, mcrotaf_divpos)
  mccomp_posa[19] = mcposaf_divpos;
  mccomp_posr[19] = mcposrf_divpos;
  mcNCounter[19]  = mcPCounter[19] = mcP2Counter[19] = 0;
  mcAbsorbProp[19]= 0;
    /* Component divhlambda_monitor_f. */
  /* Setting parameters for component divhlambda_monitor_f. */
  SIG_MESSAGE("divhlambda_monitor_f (Init:SetPar)");
#line 385 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  if("f_divv_lambda.dat") strncpy(mccdivhlambda_monitor_f_filename, "f_divv_lambda.dat" ? "f_divv_lambda.dat" : "", 16384); else mccdivhlambda_monitor_f_filename[0]='\0';
#line 57 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccdivhlambda_monitor_f_xmin = -0.05;
#line 57 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccdivhlambda_monitor_f_xmax = 0.05;
#line 57 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccdivhlambda_monitor_f_ymin = -0.05;
#line 57 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccdivhlambda_monitor_f_ymax = 0.05;
#line 386 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccdivhlambda_monitor_f_xwidth = mcipdet_width_focus;
#line 387 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccdivhlambda_monitor_f_yheight = mcipdet_width_focus;
#line 388 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccdivhlambda_monitor_f_maxdiv_h = mcipdivergence_max * 2;
#line 389 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccdivhlambda_monitor_f_Lmin = mcipL_min;
#line 390 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccdivhlambda_monitor_f_Lmax = mcipL_max;
#line 391 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccdivhlambda_monitor_f_restore_neutron = 1;
#line 57 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccdivhlambda_monitor_f_nx = 0;
#line 57 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccdivhlambda_monitor_f_ny = 0;
#line 57 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccdivhlambda_monitor_f_nz = 1;
#line 57 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
  mccdivhlambda_monitor_f_nowritefile = 0;
#line 11344 "./ess_extraction.c"

  SIG_MESSAGE("divhlambda_monitor_f (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 393 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
    (0)*DEG2RAD,
#line 393 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
    (0)*DEG2RAD,
#line 393 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
    (90)*DEG2RAD);
#line 11354 "./ess_extraction.c"
  rot_mul(mctr1, mcrotamonitor_2, mcrotadivhlambda_monitor_f);
  rot_transpose(mcrotaf_divpos, mctr1);
  rot_mul(mcrotadivhlambda_monitor_f, mctr1, mcrotrdivhlambda_monitor_f);
  mctc1 = coords_set(
#line 392 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
    0,
#line 392 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
    0,
#line 392 "/home/cherb/Documents/McStas/ess_moderator_extraction/ess_extraction.instr"
    0);
#line 11365 "./ess_extraction.c"
  rot_transpose(mcrotamonitor_2, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposadivhlambda_monitor_f = coords_add(mcposamonitor_2, mctc2);
  mctc1 = coords_sub(mcposaf_divpos, mcposadivhlambda_monitor_f);
  mcposrdivhlambda_monitor_f = rot_apply(mcrotadivhlambda_monitor_f, mctc1);
  mcDEBUG_COMPONENT("divhlambda_monitor_f", mcposadivhlambda_monitor_f, mcrotadivhlambda_monitor_f)
  mccomp_posa[20] = mcposadivhlambda_monitor_f;
  mccomp_posr[20] = mcposrdivhlambda_monitor_f;
  mcNCounter[20]  = mcPCounter[20] = mcP2Counter[20] = 0;
  mcAbsorbProp[20]= 0;
  /* Component initializations. */
  /* Initializations for component origin. */
  SIG_MESSAGE("origin (Init)");
#define mccompcurname  origin
#define mccompcurtype  Progress_bar
#define mccompcurindex 1
#define IntermediateCnts mccorigin_IntermediateCnts
#define StartTime mccorigin_StartTime
#define EndTime mccorigin_EndTime
#define CurrentTime mccorigin_CurrentTime
#define profile mccorigin_profile
#define percent mccorigin_percent
#define flag_save mccorigin_flag_save
#define minutes mccorigin_minutes
#line 57 "/usr/share/mcstas/2.6.1/tools/Python/mcrun/../mccodelib/../../../misc/Progress_bar.comp"
{
IntermediateCnts=0;
StartTime=0;
EndTime=0;
CurrentTime=0;

fprintf(stdout, "[%s] Initialize\n", mcinstrument_name);
  if (percent*mcget_ncount()/100 < 1e5) {
    percent=1e5*100.0/mcget_ncount();
  }
}
#line 11402 "./ess_extraction.c"
#undef minutes
#undef flag_save
#undef percent
#undef profile
#undef CurrentTime
#undef EndTime
#undef StartTime
#undef IntermediateCnts
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component source. */
  SIG_MESSAGE("source (Init)");

  /* Initializations for component source_div. */
  SIG_MESSAGE("source_div (Init)");
#define mccompcurname  source_div
#define mccompcurtype  Source_div
#define mccompcurindex 3
#define thetah mccsource_div_thetah
#define thetav mccsource_div_thetav
#define sigmah mccsource_div_sigmah
#define sigmav mccsource_div_sigmav
#define tan_h mccsource_div_tan_h
#define tan_v mccsource_div_tan_v
#define p_init mccsource_div_p_init
#define dist mccsource_div_dist
#define focus_xw mccsource_div_focus_xw
#define focus_yh mccsource_div_focus_yh
#define xwidth mccsource_div_xwidth
#define yheight mccsource_div_yheight
#define focus_aw mccsource_div_focus_aw
#define focus_ah mccsource_div_focus_ah
#define E0 mccsource_div_E0
#define dE mccsource_div_dE
#define lambda0 mccsource_div_lambda0
#define dlambda mccsource_div_dlambda
#define gauss mccsource_div_gauss
#define flux mccsource_div_flux
#line 72 "/usr/share/mcstas/2.6.1/tools/Python/mcrun/../mccodelib/../../../sources/Source_div.comp"
{
sigmah = DEG2RAD*focus_aw/(sqrt(8.0*log(2.0)));
  sigmav = DEG2RAD*focus_ah/(sqrt(8.0*log(2.0)));

  if (xwidth < 0 || yheight < 0 || focus_aw < 0 || focus_ah < 0) {
      printf("Source_div: %s: Error in input parameter values!\n"
             "ERROR       Exiting\n",
           NAME_CURRENT_COMP);
      exit(-1);
  }
  if ((!lambda0 && !E0 && !dE && !dlambda)) {
    printf("Source_div: %s: You must specify either a wavelength or energy range!\n ERROR - Exiting\n",
           NAME_CURRENT_COMP);
    exit(-1);
  }
  if ((!lambda0 && !dlambda && (E0 <= 0 || dE < 0 || E0-dE <= 0))
    || (!E0 && !dE && (lambda0 <= 0 || dlambda < 0 || lambda0-dlambda <= 0))) {
    printf("Source_div: %s: Unmeaningful definition of wavelength or energy range!\n ERROR - Exiting\n",
           NAME_CURRENT_COMP);
      exit(-1);
  }
  /* compute distance to next component */
  Coords ToTarget;
  double tx,ty,tz;
  ToTarget = coords_sub(POS_A_COMP_INDEX(INDEX_CURRENT_COMP+1),POS_A_CURRENT_COMP);
  ToTarget = rot_apply(ROT_A_CURRENT_COMP, ToTarget);
  coords_get(ToTarget, &tx, &ty, &tz);
  dist=sqrt(tx*tx+ty*ty+tz*tz);
  /* compute target area */
  if (dist) {
    focus_xw=dist*tan(focus_aw*DEG2RAD);
    focus_yh=dist*tan(focus_ah*DEG2RAD);
  }

  p_init  = flux*1e4*xwidth*yheight/mcget_ncount();
  if (!focus_aw || !focus_ah)
    exit(printf("Source_div: %s: Zero divergence defined. \n"
                "ERROR       Use non zero values for focus_aw and focus_ah.\n",
           NAME_CURRENT_COMP));
  p_init *= 2*fabs(DEG2RAD*focus_aw*sin(DEG2RAD*focus_ah/2));  /* solid angle */
  if (dlambda)
    p_init *= 2*dlambda;
  else if (dE)
    p_init *= 2*dE;
}
#line 11489 "./ess_extraction.c"
#undef flux
#undef gauss
#undef dlambda
#undef lambda0
#undef dE
#undef E0
#undef focus_ah
#undef focus_aw
#undef yheight
#undef xwidth
#undef focus_yh
#undef focus_xw
#undef dist
#undef p_init
#undef tan_v
#undef tan_h
#undef sigmav
#undef sigmah
#undef thetav
#undef thetah
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component psd_monitor_source. */
  SIG_MESSAGE("psd_monitor_source (Init)");
#define mccompcurname  psd_monitor_source
#define mccompcurtype  PSD_monitor
#define mccompcurindex 4
#define PSD_N mccpsd_monitor_source_PSD_N
#define PSD_p mccpsd_monitor_source_PSD_p
#define PSD_p2 mccpsd_monitor_source_PSD_p2
#define nx mccpsd_monitor_source_nx
#define ny mccpsd_monitor_source_ny
#define filename mccpsd_monitor_source_filename
#define xmin mccpsd_monitor_source_xmin
#define xmax mccpsd_monitor_source_xmax
#define ymin mccpsd_monitor_source_ymin
#define ymax mccpsd_monitor_source_ymax
#define xwidth mccpsd_monitor_source_xwidth
#define yheight mccpsd_monitor_source_yheight
#define restore_neutron mccpsd_monitor_source_restore_neutron
#line 68 "/usr/share/mcstas/2.6.1/tools/Python/mcrun/../mccodelib/../../../monitors/PSD_monitor.comp"
{
  if (xwidth  > 0) { xmax = xwidth/2;  xmin = -xmax; }
  if (yheight > 0) { ymax = yheight/2; ymin = -ymax; }

  if ((xmin >= xmax) || (ymin >= ymax)){
    printf("PSD_monitor: %s: Null detection area !\n"
           "ERROR        (xwidth,yheight,xmin,xmax,ymin,ymax). Exiting",
    NAME_CURRENT_COMP);
    exit(-1);
  }

  PSD_N = create_darr2d(nx, ny);
  PSD_p = create_darr2d(nx, ny);
  PSD_p2 = create_darr2d(nx, ny);

  int i, j;
  for (i=0; i<nx; i++){
    for (j=0; j<ny; j++){
      PSD_N[i][j] = 0;
      PSD_p[i][j] = 0;
      PSD_p2[i][j] = 0;
    }
  }
}
#line 11557 "./ess_extraction.c"
#undef restore_neutron
#undef yheight
#undef xwidth
#undef ymax
#undef ymin
#undef xmax
#undef xmin
#undef filename
#undef ny
#undef nx
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component parabolic_optic_before_guide_v. */
  SIG_MESSAGE("parabolic_optic_before_guide_v (Init)");
#define mccompcurname  parabolic_optic_before_guide_v
#define mccompcurtype  FlatEllipse_finite_mirror
#define mccompcurindex 5
#define reflect mccparabolic_optic_before_guide_v_reflect
#define s mccparabolic_optic_before_guide_v_s
#define pTable mccparabolic_optic_before_guide_v_pTable
#define R0 mccparabolic_optic_before_guide_v_R0
#define Qc mccparabolic_optic_before_guide_v_Qc
#define W mccparabolic_optic_before_guide_v_W
#define alpha mccparabolic_optic_before_guide_v_alpha
#define transmit mccparabolic_optic_before_guide_v_transmit
#define sourceDist mccparabolic_optic_before_guide_v_sourceDist
#define LStart mccparabolic_optic_before_guide_v_LStart
#define LEnd mccparabolic_optic_before_guide_v_LEnd
#define lStart mccparabolic_optic_before_guide_v_lStart
#define lEnd mccparabolic_optic_before_guide_v_lEnd
#define r_0 mccparabolic_optic_before_guide_v_r_0
#define nummirror mccparabolic_optic_before_guide_v_nummirror
#define mf mccparabolic_optic_before_guide_v_mf
#define mb mccparabolic_optic_before_guide_v_mb
#define mirror_width mccparabolic_optic_before_guide_v_mirror_width
#define doubleReflections mccparabolic_optic_before_guide_v_doubleReflections
#line 79 "FlatEllipse_finite_mirror.comp"
{
    if (sourceDist == 0){
        sourceDist = sqrt(LStart*LStart);
    }
    pointer_lStart = &lStart;
    //Load Reflectivity Data File
    /*if (reflect && strlen(reflect)) {
        if (Table_Read(&pTable, reflect, 1) <= 0)
            exit(fprintf(stderr, "Can not read file: %s\n", reflect));
    }

    //Custom function for tracing neutrons using table data for reflectivity
    void traceNeutronConicWithTables(Particle* pa, ConicSurf c) {
        double tl = getTimeOfFirstCollisionConic(*pa, c);
        if (tl < 0)
            return;
        else {
            //Move Particle to Surface Edge
            moveParticleT(tl,pa);

            if (c.m==0) {
                absorbParticle(pa);
                return;
            }

            //Handle Reflectivity
            Vec n = getNormConic(getParticlePos(*pa),c);
            double vdotn = dotVec(n,getParticleVel(*pa));

            double q = fabs(2*vdotn*V2Q);
            double B;
            if (reflect && strlen(reflect))
                B=Table_Value(pTable, q, 1);
            else {
                B = R0;
                if (q > Qc) {
                    double arg = (q-c.m*Qc)/W;
                    if(arg < 10)
                        B *= .5*(1-tanh(arg))*(1-alpha*(q-Qc));
                    else
                        B=0;
                }
            }
            if (B < 0)
                B=0;
            else if (B > 1)
                B=1;
            if (!transmit) {
                if (!B) absorbParticle(pa);
                pa->w *= B;
                reflectParticle(n,pa);
            } else {
                if (B == 0 || rand01() >= B) { /*unreflected*/ /*}
                else { reflectParticle(n,pa); }
            }
        }
    }
    */
    //Make new scene
    silicon = (mirror_width==0) ? 0 : -1; //neutron starts in air by default
    s = makeScene();
    rs = get_r_at_z0(nummirror, 0, r_0, sourceDist, LEnd, lStart, lEnd);

    //Set Scene to use custom trace function for conic
    //s.traceNeutronConic = traceNeutronConicWithTables;

    //Add Geometry Here
    Point p1;
    for (int i = 0; i < nummirror; i++) {
		    p1 = makePoint(rs[i], 0, 0);
            addFlatEllipse(LStart, LEnd, p1, lStart, lEnd, -0.1, 0.1, mf, doubleReflections, &s); //inner side of the mirror
		    printf("b[%d] = %f\n", i, rs[i]);
    }
    if (mirror_width > 0){
        for (int i = 0; i < nummirror; i++){
            p1 = makePoint(rs[i]+mirror_width, 0, 0);
            addFlatEllipse(LStart, LEnd, p1, lStart, lEnd, -0.1, 0.1, mb, doubleReflections, &s); //backside of the above mirror
        }
    }
    addDisk(lEnd, 0.0, 2.0, &s); //neutrons will be propagated important if the are in silicon
	//addEllipsoid(-L, L,p1, -l,+l, 40,&s);
}
#line 11682 "./ess_extraction.c"
#undef doubleReflections
#undef mirror_width
#undef mb
#undef mf
#undef nummirror
#undef r_0
#undef lEnd
#undef lStart
#undef LEnd
#undef LStart
#undef sourceDist
#undef transmit
#undef alpha
#undef W
#undef Qc
#undef R0
#undef pTable
#undef s
#undef reflect
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component parabolic_optic_before_guide_h. */
  SIG_MESSAGE("parabolic_optic_before_guide_h (Init)");
#define mccompcurname  parabolic_optic_before_guide_h
#define mccompcurtype  FlatEllipse_finite_mirror
#define mccompcurindex 6
#define reflect mccparabolic_optic_before_guide_h_reflect
#define s mccparabolic_optic_before_guide_h_s
#define pTable mccparabolic_optic_before_guide_h_pTable
#define R0 mccparabolic_optic_before_guide_h_R0
#define Qc mccparabolic_optic_before_guide_h_Qc
#define W mccparabolic_optic_before_guide_h_W
#define alpha mccparabolic_optic_before_guide_h_alpha
#define transmit mccparabolic_optic_before_guide_h_transmit
#define sourceDist mccparabolic_optic_before_guide_h_sourceDist
#define LStart mccparabolic_optic_before_guide_h_LStart
#define LEnd mccparabolic_optic_before_guide_h_LEnd
#define lStart mccparabolic_optic_before_guide_h_lStart
#define lEnd mccparabolic_optic_before_guide_h_lEnd
#define r_0 mccparabolic_optic_before_guide_h_r_0
#define nummirror mccparabolic_optic_before_guide_h_nummirror
#define mf mccparabolic_optic_before_guide_h_mf
#define mb mccparabolic_optic_before_guide_h_mb
#define mirror_width mccparabolic_optic_before_guide_h_mirror_width
#define doubleReflections mccparabolic_optic_before_guide_h_doubleReflections
#line 79 "FlatEllipse_finite_mirror.comp"
{
    if (sourceDist == 0){
        sourceDist = sqrt(LStart*LStart);
    }
    pointer_lStart = &lStart;
    //Load Reflectivity Data File
    /*if (reflect && strlen(reflect)) {
        if (Table_Read(&pTable, reflect, 1) <= 0)
            exit(fprintf(stderr, "Can not read file: %s\n", reflect));
    }

    //Custom function for tracing neutrons using table data for reflectivity
    void traceNeutronConicWithTables(Particle* pa, ConicSurf c) {
        double tl = getTimeOfFirstCollisionConic(*pa, c);
        if (tl < 0)
            return;
        else {
            //Move Particle to Surface Edge
            moveParticleT(tl,pa);

            if (c.m==0) {
                absorbParticle(pa);
                return;
            }

            //Handle Reflectivity
            Vec n = getNormConic(getParticlePos(*pa),c);
            double vdotn = dotVec(n,getParticleVel(*pa));

            double q = fabs(2*vdotn*V2Q);
            double B;
            if (reflect && strlen(reflect))
                B=Table_Value(pTable, q, 1);
            else {
                B = R0;
                if (q > Qc) {
                    double arg = (q-c.m*Qc)/W;
                    if(arg < 10)
                        B *= .5*(1-tanh(arg))*(1-alpha*(q-Qc));
                    else
                        B=0;
                }
            }
            if (B < 0)
                B=0;
            else if (B > 1)
                B=1;
            if (!transmit) {
                if (!B) absorbParticle(pa);
                pa->w *= B;
                reflectParticle(n,pa);
            } else {
                if (B == 0 || rand01() >= B) { /*unreflected*/ /*}
                else { reflectParticle(n,pa); }
            }
        }
    }
    */
    //Make new scene
    silicon = (mirror_width==0) ? 0 : -1; //neutron starts in air by default
    s = makeScene();
    rs = get_r_at_z0(nummirror, 0, r_0, sourceDist, LEnd, lStart, lEnd);

    //Set Scene to use custom trace function for conic
    //s.traceNeutronConic = traceNeutronConicWithTables;

    //Add Geometry Here
    Point p1;
    for (int i = 0; i < nummirror; i++) {
		    p1 = makePoint(rs[i], 0, 0);
            addFlatEllipse(LStart, LEnd, p1, lStart, lEnd, -0.1, 0.1, mf, doubleReflections, &s); //inner side of the mirror
		    printf("b[%d] = %f\n", i, rs[i]);
    }
    if (mirror_width > 0){
        for (int i = 0; i < nummirror; i++){
            p1 = makePoint(rs[i]+mirror_width, 0, 0);
            addFlatEllipse(LStart, LEnd, p1, lStart, lEnd, -0.1, 0.1, mb, doubleReflections, &s); //backside of the above mirror
        }
    }
    addDisk(lEnd, 0.0, 2.0, &s); //neutrons will be propagated important if the are in silicon
	//addEllipsoid(-L, L,p1, -l,+l, 40,&s);
}
#line 11813 "./ess_extraction.c"
#undef doubleReflections
#undef mirror_width
#undef mb
#undef mf
#undef nummirror
#undef r_0
#undef lEnd
#undef lStart
#undef LEnd
#undef LStart
#undef sourceDist
#undef transmit
#undef alpha
#undef W
#undef Qc
#undef R0
#undef pTable
#undef s
#undef reflect
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component after_optic_source. */
  SIG_MESSAGE("after_optic_source (Init)");

  /* Initializations for component psd_monitor_afteropticsource. */
  SIG_MESSAGE("psd_monitor_afteropticsource (Init)");
#define mccompcurname  psd_monitor_afteropticsource
#define mccompcurtype  PSD_monitor
#define mccompcurindex 8
#define PSD_N mccpsd_monitor_afteropticsource_PSD_N
#define PSD_p mccpsd_monitor_afteropticsource_PSD_p
#define PSD_p2 mccpsd_monitor_afteropticsource_PSD_p2
#define nx mccpsd_monitor_afteropticsource_nx
#define ny mccpsd_monitor_afteropticsource_ny
#define filename mccpsd_monitor_afteropticsource_filename
#define xmin mccpsd_monitor_afteropticsource_xmin
#define xmax mccpsd_monitor_afteropticsource_xmax
#define ymin mccpsd_monitor_afteropticsource_ymin
#define ymax mccpsd_monitor_afteropticsource_ymax
#define xwidth mccpsd_monitor_afteropticsource_xwidth
#define yheight mccpsd_monitor_afteropticsource_yheight
#define restore_neutron mccpsd_monitor_afteropticsource_restore_neutron
#line 68 "/usr/share/mcstas/2.6.1/tools/Python/mcrun/../mccodelib/../../../monitors/PSD_monitor.comp"
{
  if (xwidth  > 0) { xmax = xwidth/2;  xmin = -xmax; }
  if (yheight > 0) { ymax = yheight/2; ymin = -ymax; }

  if ((xmin >= xmax) || (ymin >= ymax)){
    printf("PSD_monitor: %s: Null detection area !\n"
           "ERROR        (xwidth,yheight,xmin,xmax,ymin,ymax). Exiting",
    NAME_CURRENT_COMP);
    exit(-1);
  }

  PSD_N = create_darr2d(nx, ny);
  PSD_p = create_darr2d(nx, ny);
  PSD_p2 = create_darr2d(nx, ny);

  int i, j;
  for (i=0; i<nx; i++){
    for (j=0; j<ny; j++){
      PSD_N[i][j] = 0;
      PSD_p[i][j] = 0;
      PSD_p2[i][j] = 0;
    }
  }
}
#line 11883 "./ess_extraction.c"
#undef restore_neutron
#undef yheight
#undef xwidth
#undef ymax
#undef ymin
#undef xmax
#undef xmin
#undef filename
#undef ny
#undef nx
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component guide_gravity_1. */
  SIG_MESSAGE("guide_gravity_1 (Init)");
#define mccompcurname  guide_gravity_1
#define mccompcurtype  Guide_gravity
#define mccompcurindex 9
#define GVars mccguide_gravity_1_GVars
#define pTable mccguide_gravity_1_pTable
#define w1 mccguide_gravity_1_w1
#define h1 mccguide_gravity_1_h1
#define w2 mccguide_gravity_1_w2
#define h2 mccguide_gravity_1_h2
#define l mccguide_gravity_1_l
#define R0 mccguide_gravity_1_R0
#define Qc mccguide_gravity_1_Qc
#define alpha mccguide_gravity_1_alpha
#define m mccguide_gravity_1_m
#define W mccguide_gravity_1_W
#define nslit mccguide_gravity_1_nslit
#define d mccguide_gravity_1_d
#define mleft mccguide_gravity_1_mleft
#define mright mccguide_gravity_1_mright
#define mtop mccguide_gravity_1_mtop
#define mbottom mccguide_gravity_1_mbottom
#define nhslit mccguide_gravity_1_nhslit
#define G mccguide_gravity_1_G
#define aleft mccguide_gravity_1_aleft
#define aright mccguide_gravity_1_aright
#define atop mccguide_gravity_1_atop
#define abottom mccguide_gravity_1_abottom
#define wavy mccguide_gravity_1_wavy
#define wavy_z mccguide_gravity_1_wavy_z
#define wavy_tb mccguide_gravity_1_wavy_tb
#define wavy_lr mccguide_gravity_1_wavy_lr
#define chamfers mccguide_gravity_1_chamfers
#define chamfers_z mccguide_gravity_1_chamfers_z
#define chamfers_lr mccguide_gravity_1_chamfers_lr
#define chamfers_tb mccguide_gravity_1_chamfers_tb
#define nelements mccguide_gravity_1_nelements
#define nu mccguide_gravity_1_nu
#define phase mccguide_gravity_1_phase
#define reflect mccguide_gravity_1_reflect
#line 339 "/usr/share/mcstas/2.6.1/tools/Python/mcrun/../mccodelib/../../../optics/Guide_gravity.comp"
{
  double Gx=0, Gy=-GRAVITY, Gz=0;
  Coords mcLocG;
  int i;

  if (reflect && strlen(reflect) && strcmp(reflect,"NULL") && strcmp(reflect,"0")) {
    if (Table_Read(&pTable, reflect, 1) <= 0) /* read 1st block data from file into pTable */
      exit(fprintf(stderr,"Guide_gravity: %s: can not read file %s\n", NAME_CURRENT_COMP, reflect));
  } else {
    if (W < 0 || R0 < 0 || Qc < 0)
    { fprintf(stderr,"Guide_gravity: %s: W R0 Qc must be >0.\n", NAME_CURRENT_COMP);
      exit(-1); }
  }

  if (nslit <= 0 || nhslit <= 0)
  { fprintf(stderr,"Guide_gravity: %s: nslit nhslit must be >0.\n", NAME_CURRENT_COMP);
    exit(-1); }

  if (!w1 || !h1)
  { fprintf(stderr,"Guide_gravity: %s: input window is closed (w1=h1=0).\n", NAME_CURRENT_COMP);
    exit(-1); }

  if (d*nslit > w1) exit(fprintf(stderr, "Guide_gravity: %s: absorbing walls fill input window. No space left for transmission (d*nslit > w1).\n", NAME_CURRENT_COMP));

  if (!w2) w2=w1;
  if (!h2) h2=h1;

  if (mcgravitation) G=-GRAVITY;
  mcLocG = rot_apply(ROT_A_CURRENT_COMP, coords_set(0,G,0));
  coords_get(mcLocG, &Gx, &Gy, &Gz);

  strcpy(GVars.compcurname, NAME_CURRENT_COMP);

  if (l > 0 && nelements > 0) {

    Gravity_guide_Init(&GVars,
      w1, h1, w2, h2, l, R0,
      Qc, alpha, m, W, nslit, d,
      Gx, Gy, Gz, mleft, mright, mtop,
      mbottom, nhslit, wavy_lr, wavy_tb, wavy_z, wavy,
      chamfers_z, chamfers_lr, chamfers_tb, chamfers,nu,phase,aleft,aright,atop,abottom);
    if (!G) for (i=0; i<5; GVars.A[i++] = 0);
    if (GVars.fc_freq != 0 || GVars.fc_phase != 0) {
      if (w1 != w2 || h1 != h2)
      exit(fprintf(stderr,"Guide_gravity: %s: rotating slit pack must be straight (w1=w2 and h1=h2).\n", NAME_CURRENT_COMP));
      printf("Guide_gravity: %s: Fermi Chopper mode: frequency=%g [Hz] phase=%g [deg]\n",
        NAME_CURRENT_COMP, GVars.fc_freq, GVars.fc_phase);
    }
  } else printf("Guide_gravity: %s: unactivated (l=0 or nelements=0)\n", NAME_CURRENT_COMP);

}
#line 11994 "./ess_extraction.c"
#undef reflect
#undef phase
#undef nu
#undef nelements
#undef chamfers_tb
#undef chamfers_lr
#undef chamfers_z
#undef chamfers
#undef wavy_lr
#undef wavy_tb
#undef wavy_z
#undef wavy
#undef abottom
#undef atop
#undef aright
#undef aleft
#undef G
#undef nhslit
#undef mbottom
#undef mtop
#undef mright
#undef mleft
#undef d
#undef nslit
#undef W
#undef m
#undef alpha
#undef Qc
#undef R0
#undef l
#undef h2
#undef w2
#undef h1
#undef w1
#undef pTable
#undef GVars
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component monitor_1. */
  SIG_MESSAGE("monitor_1 (Init)");

  /* Initializations for component divpos_monitor. */
  SIG_MESSAGE("divpos_monitor (Init)");
#define mccompcurname  divpos_monitor
#define mccompcurtype  DivPos_monitor
#define mccompcurindex 11
#define nh mccdivpos_monitor_nh
#define ndiv mccdivpos_monitor_ndiv
#define Div_N mccdivpos_monitor_Div_N
#define Div_p mccdivpos_monitor_Div_p
#define Div_p2 mccdivpos_monitor_Div_p2
#define filename mccdivpos_monitor_filename
#define xmin mccdivpos_monitor_xmin
#define xmax mccdivpos_monitor_xmax
#define ymin mccdivpos_monitor_ymin
#define ymax mccdivpos_monitor_ymax
#define xwidth mccdivpos_monitor_xwidth
#define yheight mccdivpos_monitor_yheight
#define maxdiv_h mccdivpos_monitor_maxdiv_h
#define restore_neutron mccdivpos_monitor_restore_neutron
#define nx mccdivpos_monitor_nx
#define ny mccdivpos_monitor_ny
#define nz mccdivpos_monitor_nz
#define nowritefile mccdivpos_monitor_nowritefile
#line 68 "/usr/share/mcstas/2.6.1/tools/Python/mcrun/../mccodelib/../../../monitors/DivPos_monitor.comp"
{
int i,j;

if (xwidth  > 0) { xmax = xwidth/2;  xmin = -xmax; }
    if (yheight > 0) { ymax = yheight/2; ymin = -ymax; }

    if ((xmin >= xmax) || (ymin >= ymax)) {
            printf("DivPos_monitor: %s: Null detection area !\n"
                   "ERROR           (xwidth,yheight,xmin,xmax,ymin,ymax). Exiting",
           NAME_CURRENT_COMP);
      exit(-1);
    }

    for (i=0; i<nh; i++)
     for (j=0; j<ndiv; j++)
     {
      Div_N[i][j] = 0;
      Div_p[i][j] = 0;
      Div_p2[i][j] = 0;
     }
    NORM(nx,ny,nz);
}
#line 12084 "./ess_extraction.c"
#undef nowritefile
#undef nz
#undef ny
#undef nx
#undef restore_neutron
#undef maxdiv_h
#undef yheight
#undef xwidth
#undef ymax
#undef ymin
#undef xmax
#undef xmin
#undef filename
#undef Div_p2
#undef Div_p
#undef Div_N
#undef ndiv
#undef nh
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component psd_monitor_g1. */
  SIG_MESSAGE("psd_monitor_g1 (Init)");
#define mccompcurname  psd_monitor_g1
#define mccompcurtype  PSD_monitor
#define mccompcurindex 12
#define PSD_N mccpsd_monitor_g1_PSD_N
#define PSD_p mccpsd_monitor_g1_PSD_p
#define PSD_p2 mccpsd_monitor_g1_PSD_p2
#define nx mccpsd_monitor_g1_nx
#define ny mccpsd_monitor_g1_ny
#define filename mccpsd_monitor_g1_filename
#define xmin mccpsd_monitor_g1_xmin
#define xmax mccpsd_monitor_g1_xmax
#define ymin mccpsd_monitor_g1_ymin
#define ymax mccpsd_monitor_g1_ymax
#define xwidth mccpsd_monitor_g1_xwidth
#define yheight mccpsd_monitor_g1_yheight
#define restore_neutron mccpsd_monitor_g1_restore_neutron
#line 68 "/usr/share/mcstas/2.6.1/tools/Python/mcrun/../mccodelib/../../../monitors/PSD_monitor.comp"
{
  if (xwidth  > 0) { xmax = xwidth/2;  xmin = -xmax; }
  if (yheight > 0) { ymax = yheight/2; ymin = -ymax; }

  if ((xmin >= xmax) || (ymin >= ymax)){
    printf("PSD_monitor: %s: Null detection area !\n"
           "ERROR        (xwidth,yheight,xmin,xmax,ymin,ymax). Exiting",
    NAME_CURRENT_COMP);
    exit(-1);
  }

  PSD_N = create_darr2d(nx, ny);
  PSD_p = create_darr2d(nx, ny);
  PSD_p2 = create_darr2d(nx, ny);

  int i, j;
  for (i=0; i<nx; i++){
    for (j=0; j<ny; j++){
      PSD_N[i][j] = 0;
      PSD_p[i][j] = 0;
      PSD_p2[i][j] = 0;
    }
  }
}
#line 12150 "./ess_extraction.c"
#undef restore_neutron
#undef yheight
#undef xwidth
#undef ymax
#undef ymin
#undef xmax
#undef xmin
#undef filename
#undef ny
#undef nx
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component divhlambda_monitor_g1. */
  SIG_MESSAGE("divhlambda_monitor_g1 (Init)");
#define mccompcurname  divhlambda_monitor_g1
#define mccompcurtype  DivLambda_monitor
#define mccompcurindex 13
#define nL mccdivhlambda_monitor_g1_nL
#define nh mccdivhlambda_monitor_g1_nh
#define Div_N mccdivhlambda_monitor_g1_Div_N
#define Div_p mccdivhlambda_monitor_g1_Div_p
#define Div_p2 mccdivhlambda_monitor_g1_Div_p2
#define filename mccdivhlambda_monitor_g1_filename
#define xmin mccdivhlambda_monitor_g1_xmin
#define xmax mccdivhlambda_monitor_g1_xmax
#define ymin mccdivhlambda_monitor_g1_ymin
#define ymax mccdivhlambda_monitor_g1_ymax
#define xwidth mccdivhlambda_monitor_g1_xwidth
#define yheight mccdivhlambda_monitor_g1_yheight
#define maxdiv_h mccdivhlambda_monitor_g1_maxdiv_h
#define Lmin mccdivhlambda_monitor_g1_Lmin
#define Lmax mccdivhlambda_monitor_g1_Lmax
#define restore_neutron mccdivhlambda_monitor_g1_restore_neutron
#define nx mccdivhlambda_monitor_g1_nx
#define ny mccdivhlambda_monitor_g1_ny
#define nz mccdivhlambda_monitor_g1_nz
#define nowritefile mccdivhlambda_monitor_g1_nowritefile
#line 67 "/usr/share/mcstas/2.6.1/tools/Python/mcrun/../mccodelib/../../../monitors/DivLambda_monitor.comp"
{
  int i,j;

  if (xwidth  > 0) { xmax = xwidth/2;  xmin = -xmax; }
  if (yheight > 0) { ymax = yheight/2; ymin = -ymax; }
  if ((xmin >= xmax) || (ymin >= ymax)) {
    printf("ERROR: (%s): Null detection area! Aborting.\n",
        NAME_CURRENT_COMP);
    exit(-1);
  }

  for (i=0; i<nL; i++)
    for (j=0; j<nh; j++)
    {
      Div_N[i][j] = 0;
      Div_p[i][j] = 0;
      Div_p2[i][j] = 0;
    }
  NORM(nx,ny,nz);

}
#line 12215 "./ess_extraction.c"
#undef nowritefile
#undef nz
#undef ny
#undef nx
#undef restore_neutron
#undef Lmax
#undef Lmin
#undef maxdiv_h
#undef yheight
#undef xwidth
#undef ymax
#undef ymin
#undef xmax
#undef xmin
#undef filename
#undef Div_p2
#undef Div_p
#undef Div_N
#undef nh
#undef nL
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component parabolic_optic1. */
  SIG_MESSAGE("parabolic_optic1 (Init)");
#define mccompcurname  parabolic_optic1
#define mccompcurtype  FlatEllipse_finite_mirror
#define mccompcurindex 14
#define reflect mccparabolic_optic1_reflect
#define s mccparabolic_optic1_s
#define pTable mccparabolic_optic1_pTable
#define R0 mccparabolic_optic1_R0
#define Qc mccparabolic_optic1_Qc
#define W mccparabolic_optic1_W
#define alpha mccparabolic_optic1_alpha
#define transmit mccparabolic_optic1_transmit
#define sourceDist mccparabolic_optic1_sourceDist
#define LStart mccparabolic_optic1_LStart
#define LEnd mccparabolic_optic1_LEnd
#define lStart mccparabolic_optic1_lStart
#define lEnd mccparabolic_optic1_lEnd
#define r_0 mccparabolic_optic1_r_0
#define nummirror mccparabolic_optic1_nummirror
#define mf mccparabolic_optic1_mf
#define mb mccparabolic_optic1_mb
#define mirror_width mccparabolic_optic1_mirror_width
#define doubleReflections mccparabolic_optic1_doubleReflections
#line 79 "FlatEllipse_finite_mirror.comp"
{
    if (sourceDist == 0){
        sourceDist = sqrt(LStart*LStart);
    }
    pointer_lStart = &lStart;
    //Load Reflectivity Data File
    /*if (reflect && strlen(reflect)) {
        if (Table_Read(&pTable, reflect, 1) <= 0)
            exit(fprintf(stderr, "Can not read file: %s\n", reflect));
    }

    //Custom function for tracing neutrons using table data for reflectivity
    void traceNeutronConicWithTables(Particle* pa, ConicSurf c) {
        double tl = getTimeOfFirstCollisionConic(*pa, c);
        if (tl < 0)
            return;
        else {
            //Move Particle to Surface Edge
            moveParticleT(tl,pa);

            if (c.m==0) {
                absorbParticle(pa);
                return;
            }

            //Handle Reflectivity
            Vec n = getNormConic(getParticlePos(*pa),c);
            double vdotn = dotVec(n,getParticleVel(*pa));

            double q = fabs(2*vdotn*V2Q);
            double B;
            if (reflect && strlen(reflect))
                B=Table_Value(pTable, q, 1);
            else {
                B = R0;
                if (q > Qc) {
                    double arg = (q-c.m*Qc)/W;
                    if(arg < 10)
                        B *= .5*(1-tanh(arg))*(1-alpha*(q-Qc));
                    else
                        B=0;
                }
            }
            if (B < 0)
                B=0;
            else if (B > 1)
                B=1;
            if (!transmit) {
                if (!B) absorbParticle(pa);
                pa->w *= B;
                reflectParticle(n,pa);
            } else {
                if (B == 0 || rand01() >= B) { /*unreflected*/ /*}
                else { reflectParticle(n,pa); }
            }
        }
    }
    */
    //Make new scene
    silicon = (mirror_width==0) ? 0 : -1; //neutron starts in air by default
    s = makeScene();
    rs = get_r_at_z0(nummirror, 0, r_0, sourceDist, LEnd, lStart, lEnd);

    //Set Scene to use custom trace function for conic
    //s.traceNeutronConic = traceNeutronConicWithTables;

    //Add Geometry Here
    Point p1;
    for (int i = 0; i < nummirror; i++) {
		    p1 = makePoint(rs[i], 0, 0);
            addFlatEllipse(LStart, LEnd, p1, lStart, lEnd, -0.1, 0.1, mf, doubleReflections, &s); //inner side of the mirror
		    printf("b[%d] = %f\n", i, rs[i]);
    }
    if (mirror_width > 0){
        for (int i = 0; i < nummirror; i++){
            p1 = makePoint(rs[i]+mirror_width, 0, 0);
            addFlatEllipse(LStart, LEnd, p1, lStart, lEnd, -0.1, 0.1, mb, doubleReflections, &s); //backside of the above mirror
        }
    }
    addDisk(lEnd, 0.0, 2.0, &s); //neutrons will be propagated important if the are in silicon
	//addEllipsoid(-L, L,p1, -l,+l, 40,&s);
}
#line 12347 "./ess_extraction.c"
#undef doubleReflections
#undef mirror_width
#undef mb
#undef mf
#undef nummirror
#undef r_0
#undef lEnd
#undef lStart
#undef LEnd
#undef LStart
#undef sourceDist
#undef transmit
#undef alpha
#undef W
#undef Qc
#undef R0
#undef pTable
#undef s
#undef reflect
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component parabolic_optic2. */
  SIG_MESSAGE("parabolic_optic2 (Init)");
#define mccompcurname  parabolic_optic2
#define mccompcurtype  FlatEllipse_finite_mirror
#define mccompcurindex 15
#define reflect mccparabolic_optic2_reflect
#define s mccparabolic_optic2_s
#define pTable mccparabolic_optic2_pTable
#define R0 mccparabolic_optic2_R0
#define Qc mccparabolic_optic2_Qc
#define W mccparabolic_optic2_W
#define alpha mccparabolic_optic2_alpha
#define transmit mccparabolic_optic2_transmit
#define sourceDist mccparabolic_optic2_sourceDist
#define LStart mccparabolic_optic2_LStart
#define LEnd mccparabolic_optic2_LEnd
#define lStart mccparabolic_optic2_lStart
#define lEnd mccparabolic_optic2_lEnd
#define r_0 mccparabolic_optic2_r_0
#define nummirror mccparabolic_optic2_nummirror
#define mf mccparabolic_optic2_mf
#define mb mccparabolic_optic2_mb
#define mirror_width mccparabolic_optic2_mirror_width
#define doubleReflections mccparabolic_optic2_doubleReflections
#line 79 "FlatEllipse_finite_mirror.comp"
{
    if (sourceDist == 0){
        sourceDist = sqrt(LStart*LStart);
    }
    pointer_lStart = &lStart;
    //Load Reflectivity Data File
    /*if (reflect && strlen(reflect)) {
        if (Table_Read(&pTable, reflect, 1) <= 0)
            exit(fprintf(stderr, "Can not read file: %s\n", reflect));
    }

    //Custom function for tracing neutrons using table data for reflectivity
    void traceNeutronConicWithTables(Particle* pa, ConicSurf c) {
        double tl = getTimeOfFirstCollisionConic(*pa, c);
        if (tl < 0)
            return;
        else {
            //Move Particle to Surface Edge
            moveParticleT(tl,pa);

            if (c.m==0) {
                absorbParticle(pa);
                return;
            }

            //Handle Reflectivity
            Vec n = getNormConic(getParticlePos(*pa),c);
            double vdotn = dotVec(n,getParticleVel(*pa));

            double q = fabs(2*vdotn*V2Q);
            double B;
            if (reflect && strlen(reflect))
                B=Table_Value(pTable, q, 1);
            else {
                B = R0;
                if (q > Qc) {
                    double arg = (q-c.m*Qc)/W;
                    if(arg < 10)
                        B *= .5*(1-tanh(arg))*(1-alpha*(q-Qc));
                    else
                        B=0;
                }
            }
            if (B < 0)
                B=0;
            else if (B > 1)
                B=1;
            if (!transmit) {
                if (!B) absorbParticle(pa);
                pa->w *= B;
                reflectParticle(n,pa);
            } else {
                if (B == 0 || rand01() >= B) { /*unreflected*/ /*}
                else { reflectParticle(n,pa); }
            }
        }
    }
    */
    //Make new scene
    silicon = (mirror_width==0) ? 0 : -1; //neutron starts in air by default
    s = makeScene();
    rs = get_r_at_z0(nummirror, 0, r_0, sourceDist, LEnd, lStart, lEnd);

    //Set Scene to use custom trace function for conic
    //s.traceNeutronConic = traceNeutronConicWithTables;

    //Add Geometry Here
    Point p1;
    for (int i = 0; i < nummirror; i++) {
		    p1 = makePoint(rs[i], 0, 0);
            addFlatEllipse(LStart, LEnd, p1, lStart, lEnd, -0.1, 0.1, mf, doubleReflections, &s); //inner side of the mirror
		    printf("b[%d] = %f\n", i, rs[i]);
    }
    if (mirror_width > 0){
        for (int i = 0; i < nummirror; i++){
            p1 = makePoint(rs[i]+mirror_width, 0, 0);
            addFlatEllipse(LStart, LEnd, p1, lStart, lEnd, -0.1, 0.1, mb, doubleReflections, &s); //backside of the above mirror
        }
    }
    addDisk(lEnd, 0.0, 2.0, &s); //neutrons will be propagated important if the are in silicon
	//addEllipsoid(-L, L,p1, -l,+l, 40,&s);
}
#line 12478 "./ess_extraction.c"
#undef doubleReflections
#undef mirror_width
#undef mb
#undef mf
#undef nummirror
#undef r_0
#undef lEnd
#undef lStart
#undef LEnd
#undef LStart
#undef sourceDist
#undef transmit
#undef alpha
#undef W
#undef Qc
#undef R0
#undef pTable
#undef s
#undef reflect
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component monitor_2. */
  SIG_MESSAGE("monitor_2 (Init)");

  /* Initializations for component psd_monitor_f_zoom. */
  SIG_MESSAGE("psd_monitor_f_zoom (Init)");
#define mccompcurname  psd_monitor_f_zoom
#define mccompcurtype  PSD_monitor
#define mccompcurindex 17
#define PSD_N mccpsd_monitor_f_zoom_PSD_N
#define PSD_p mccpsd_monitor_f_zoom_PSD_p
#define PSD_p2 mccpsd_monitor_f_zoom_PSD_p2
#define nx mccpsd_monitor_f_zoom_nx
#define ny mccpsd_monitor_f_zoom_ny
#define filename mccpsd_monitor_f_zoom_filename
#define xmin mccpsd_monitor_f_zoom_xmin
#define xmax mccpsd_monitor_f_zoom_xmax
#define ymin mccpsd_monitor_f_zoom_ymin
#define ymax mccpsd_monitor_f_zoom_ymax
#define xwidth mccpsd_monitor_f_zoom_xwidth
#define yheight mccpsd_monitor_f_zoom_yheight
#define restore_neutron mccpsd_monitor_f_zoom_restore_neutron
#line 68 "/usr/share/mcstas/2.6.1/tools/Python/mcrun/../mccodelib/../../../monitors/PSD_monitor.comp"
{
  if (xwidth  > 0) { xmax = xwidth/2;  xmin = -xmax; }
  if (yheight > 0) { ymax = yheight/2; ymin = -ymax; }

  if ((xmin >= xmax) || (ymin >= ymax)){
    printf("PSD_monitor: %s: Null detection area !\n"
           "ERROR        (xwidth,yheight,xmin,xmax,ymin,ymax). Exiting",
    NAME_CURRENT_COMP);
    exit(-1);
  }

  PSD_N = create_darr2d(nx, ny);
  PSD_p = create_darr2d(nx, ny);
  PSD_p2 = create_darr2d(nx, ny);

  int i, j;
  for (i=0; i<nx; i++){
    for (j=0; j<ny; j++){
      PSD_N[i][j] = 0;
      PSD_p[i][j] = 0;
      PSD_p2[i][j] = 0;
    }
  }
}
#line 12548 "./ess_extraction.c"
#undef restore_neutron
#undef yheight
#undef xwidth
#undef ymax
#undef ymin
#undef xmax
#undef xmin
#undef filename
#undef ny
#undef nx
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component psd_monitor_f. */
  SIG_MESSAGE("psd_monitor_f (Init)");
#define mccompcurname  psd_monitor_f
#define mccompcurtype  PSD_monitor
#define mccompcurindex 18
#define PSD_N mccpsd_monitor_f_PSD_N
#define PSD_p mccpsd_monitor_f_PSD_p
#define PSD_p2 mccpsd_monitor_f_PSD_p2
#define nx mccpsd_monitor_f_nx
#define ny mccpsd_monitor_f_ny
#define filename mccpsd_monitor_f_filename
#define xmin mccpsd_monitor_f_xmin
#define xmax mccpsd_monitor_f_xmax
#define ymin mccpsd_monitor_f_ymin
#define ymax mccpsd_monitor_f_ymax
#define xwidth mccpsd_monitor_f_xwidth
#define yheight mccpsd_monitor_f_yheight
#define restore_neutron mccpsd_monitor_f_restore_neutron
#line 68 "/usr/share/mcstas/2.6.1/tools/Python/mcrun/../mccodelib/../../../monitors/PSD_monitor.comp"
{
  if (xwidth  > 0) { xmax = xwidth/2;  xmin = -xmax; }
  if (yheight > 0) { ymax = yheight/2; ymin = -ymax; }

  if ((xmin >= xmax) || (ymin >= ymax)){
    printf("PSD_monitor: %s: Null detection area !\n"
           "ERROR        (xwidth,yheight,xmin,xmax,ymin,ymax). Exiting",
    NAME_CURRENT_COMP);
    exit(-1);
  }

  PSD_N = create_darr2d(nx, ny);
  PSD_p = create_darr2d(nx, ny);
  PSD_p2 = create_darr2d(nx, ny);

  int i, j;
  for (i=0; i<nx; i++){
    for (j=0; j<ny; j++){
      PSD_N[i][j] = 0;
      PSD_p[i][j] = 0;
      PSD_p2[i][j] = 0;
    }
  }
}
#line 12609 "./ess_extraction.c"
#undef restore_neutron
#undef yheight
#undef xwidth
#undef ymax
#undef ymin
#undef xmax
#undef xmin
#undef filename
#undef ny
#undef nx
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component f_divpos. */
  SIG_MESSAGE("f_divpos (Init)");
#define mccompcurname  f_divpos
#define mccompcurtype  DivPos_monitor
#define mccompcurindex 19
#define nh mccf_divpos_nh
#define ndiv mccf_divpos_ndiv
#define Div_N mccf_divpos_Div_N
#define Div_p mccf_divpos_Div_p
#define Div_p2 mccf_divpos_Div_p2
#define filename mccf_divpos_filename
#define xmin mccf_divpos_xmin
#define xmax mccf_divpos_xmax
#define ymin mccf_divpos_ymin
#define ymax mccf_divpos_ymax
#define xwidth mccf_divpos_xwidth
#define yheight mccf_divpos_yheight
#define maxdiv_h mccf_divpos_maxdiv_h
#define restore_neutron mccf_divpos_restore_neutron
#define nx mccf_divpos_nx
#define ny mccf_divpos_ny
#define nz mccf_divpos_nz
#define nowritefile mccf_divpos_nowritefile
#line 68 "/usr/share/mcstas/2.6.1/tools/Python/mcrun/../mccodelib/../../../monitors/DivPos_monitor.comp"
{
int i,j;

if (xwidth  > 0) { xmax = xwidth/2;  xmin = -xmax; }
    if (yheight > 0) { ymax = yheight/2; ymin = -ymax; }

    if ((xmin >= xmax) || (ymin >= ymax)) {
            printf("DivPos_monitor: %s: Null detection area !\n"
                   "ERROR           (xwidth,yheight,xmin,xmax,ymin,ymax). Exiting",
           NAME_CURRENT_COMP);
      exit(-1);
    }

    for (i=0; i<nh; i++)
     for (j=0; j<ndiv; j++)
     {
      Div_N[i][j] = 0;
      Div_p[i][j] = 0;
      Div_p2[i][j] = 0;
     }
    NORM(nx,ny,nz);
}
#line 12673 "./ess_extraction.c"
#undef nowritefile
#undef nz
#undef ny
#undef nx
#undef restore_neutron
#undef maxdiv_h
#undef yheight
#undef xwidth
#undef ymax
#undef ymin
#undef xmax
#undef xmin
#undef filename
#undef Div_p2
#undef Div_p
#undef Div_N
#undef ndiv
#undef nh
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component divhlambda_monitor_f. */
  SIG_MESSAGE("divhlambda_monitor_f (Init)");
#define mccompcurname  divhlambda_monitor_f
#define mccompcurtype  DivLambda_monitor
#define mccompcurindex 20
#define nL mccdivhlambda_monitor_f_nL
#define nh mccdivhlambda_monitor_f_nh
#define Div_N mccdivhlambda_monitor_f_Div_N
#define Div_p mccdivhlambda_monitor_f_Div_p
#define Div_p2 mccdivhlambda_monitor_f_Div_p2
#define filename mccdivhlambda_monitor_f_filename
#define xmin mccdivhlambda_monitor_f_xmin
#define xmax mccdivhlambda_monitor_f_xmax
#define ymin mccdivhlambda_monitor_f_ymin
#define ymax mccdivhlambda_monitor_f_ymax
#define xwidth mccdivhlambda_monitor_f_xwidth
#define yheight mccdivhlambda_monitor_f_yheight
#define maxdiv_h mccdivhlambda_monitor_f_maxdiv_h
#define Lmin mccdivhlambda_monitor_f_Lmin
#define Lmax mccdivhlambda_monitor_f_Lmax
#define restore_neutron mccdivhlambda_monitor_f_restore_neutron
#define nx mccdivhlambda_monitor_f_nx
#define ny mccdivhlambda_monitor_f_ny
#define nz mccdivhlambda_monitor_f_nz
#define nowritefile mccdivhlambda_monitor_f_nowritefile
#line 67 "/usr/share/mcstas/2.6.1/tools/Python/mcrun/../mccodelib/../../../monitors/DivLambda_monitor.comp"
{
  int i,j;

  if (xwidth  > 0) { xmax = xwidth/2;  xmin = -xmax; }
  if (yheight > 0) { ymax = yheight/2; ymin = -ymax; }
  if ((xmin >= xmax) || (ymin >= ymax)) {
    printf("ERROR: (%s): Null detection area! Aborting.\n",
        NAME_CURRENT_COMP);
    exit(-1);
  }

  for (i=0; i<nL; i++)
    for (j=0; j<nh; j++)
    {
      Div_N[i][j] = 0;
      Div_p[i][j] = 0;
      Div_p2[i][j] = 0;
    }
  NORM(nx,ny,nz);

}
#line 12743 "./ess_extraction.c"
#undef nowritefile
#undef nz
#undef ny
#undef nx
#undef restore_neutron
#undef Lmax
#undef Lmin
#undef maxdiv_h
#undef yheight
#undef xwidth
#undef ymax
#undef ymin
#undef xmax
#undef xmin
#undef filename
#undef Div_p2
#undef Div_p
#undef Div_N
#undef nh
#undef nL
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

    if(mcdotrace) mcdisplay();
    mcDEBUG_INSTR_END()
  }

} /* end init */

void mcraytrace(void) {
  /* Neutronics-specific defines */
#ifdef NEUTRONICS
extern double mcnx, mcny, mcnz, mcnvx, mcnvy, mcnvz;
extern double mcnt, mcnsx, mcnsy, mcnsz, mcnp;
#endif
  /* End of Neutronics-specific defines */
  /* Copy neutron state to local variables. */
  MCNUM mcnlx = mcnx;
  MCNUM mcnly = mcny;
  MCNUM mcnlz = mcnz;
  MCNUM mcnlvx = mcnvx;
  MCNUM mcnlvy = mcnvy;
  MCNUM mcnlvz = mcnvz;
  MCNUM mcnlt = mcnt;
  MCNUM mcnlsx = mcnsx;
  MCNUM mcnlsy = mcnsy;
  MCNUM mcnlsz = mcnsz;
  MCNUM mcnlp = mcnp;

  mcDEBUG_ENTER()
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define mcabsorb mcabsorbAll
  /* TRACE Component origin [1] */
  mccoordschange(mcposrorigin, mcrotrorigin,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component origin (without coords transformations) */
  mcJumpTrace_origin:
  SIG_MESSAGE("origin (Trace)");
  mcDEBUG_COMP("origin")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp

#define mcabsorbComp mcabsorbComporigin
  STORE_NEUTRON(1,
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp);
  mcScattered=0;
  mcRestore=0;
  mcNCounter[1]++;
  mcPCounter[1] += p;
  mcP2Counter[1] += p*p;
#define mccompcurname  origin
#define mccompcurtype  Progress_bar
#define mccompcurindex 1
#define IntermediateCnts mccorigin_IntermediateCnts
#define StartTime mccorigin_StartTime
#define EndTime mccorigin_EndTime
#define CurrentTime mccorigin_CurrentTime
{   /* Declarations of origin=Progress_bar() SETTING parameters. */
char* profile = mccorigin_profile;
MCNUM percent = mccorigin_percent;
MCNUM flag_save = mccorigin_flag_save;
MCNUM minutes = mccorigin_minutes;
#line 70 "/usr/share/mcstas/2.6.1/tools/Python/mcrun/../mccodelib/../../../misc/Progress_bar.comp"
{
  double ncount;
  ncount = mcget_run_num();
  if (!StartTime) {
    time(&StartTime); /* compute starting time */
    IntermediateCnts = 1e3;
  }
  time_t NowTime;
  time(&NowTime);
  /* compute initial estimate of computation duration */
  if (!EndTime && ncount >= IntermediateCnts) {
    CurrentTime = NowTime;
    if (difftime(NowTime,StartTime) > 10 && ncount) { /* wait 10 sec before writing ETA */
      EndTime = StartTime + (time_t)(difftime(NowTime,StartTime)
				     *(double)mcget_ncount()/ncount);
      IntermediateCnts = 0;
      fprintf(stdout, "\nTrace ETA ");
      if (difftime(EndTime,StartTime) < 60.0)
        fprintf(stdout, "%g [s] %% ", difftime(EndTime,StartTime));
      else if (difftime(EndTime,StartTime) > 3600.0)
        fprintf(stdout, "%g [h] %% ", difftime(EndTime,StartTime)/3600.0);
      else
        fprintf(stdout, "%g [min] %% ", difftime(EndTime,StartTime)/60.0);
    } else IntermediateCnts += 1e3;
    fflush(stdout);
  }

  /* display percentage when percent or minutes have reached step */
  if (EndTime && mcget_ncount() &&
    (    (minutes && difftime(NowTime,CurrentTime) > minutes*60)
      || (percent && !minutes && ncount >= IntermediateCnts))   )
  {
    fprintf(stdout, "%d ", (int)(ncount*100.0/mcget_ncount())); fflush(stdout);
    CurrentTime = NowTime;

    IntermediateCnts = ncount + percent*mcget_ncount()/100;
    /* check that next intermediate ncount check is a multiple of the desired percentage */
    IntermediateCnts = floor(IntermediateCnts*100/percent/mcget_ncount())*percent*mcget_ncount()/100;
    /* raise flag to indicate that we did something */
    SCATTER;
    if (flag_save) mcsave(NULL);
  }
}
#line 12921 "./ess_extraction.c"
}   /* End of origin=Progress_bar() SETTING parameter declarations. */
#undef CurrentTime
#undef EndTime
#undef StartTime
#undef IntermediateCnts
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbComporigin:
  if (RESTORE) /* restore if needed */
  { RESTORE_NEUTRON(1,
      mcnlx,
      mcnly,
      mcnlz,
      mcnlvx,
      mcnlvy,
      mcnlvz,
      mcnlt,
      mcnlsx,
      mcnlsy,
      mcnlsz,
      mcnlp); }
#undef mcabsorbComp
#undef p
#undef sz
#undef sy
#undef sx
#undef t
#undef vz
#undef vy
#undef vx
#undef z
#undef y
#undef x
  mcDEBUG_STATE(
mcnlx,
mcnly,
mcnlz,
mcnlvx,
mcnlvy,
mcnlvz,
mcnlt,
mcnlsx,
mcnlsy,
mcnlsz,
mcnlp)

  /* TRACE Component source [2] */
  mccoordschange(mcposrsource, mcrotrsource,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component source (without coords transformations) */
  mcJumpTrace_source:
  SIG_MESSAGE("source (Trace)");
  mcDEBUG_COMP("source")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp

#define mcabsorbComp mcabsorbCompsource
  STORE_NEUTRON(2,
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp);
  mcScattered=0;
  mcRestore=0;
  mcNCounter[2]++;
  mcPCounter[2] += p;
  mcP2Counter[2] += p*p;
#define mccompcurname  source
#define mccompcurtype  Arm
#define mccompcurindex 2
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompsource:
  if (RESTORE) /* restore if needed */
  { RESTORE_NEUTRON(2,
      mcnlx,
      mcnly,
      mcnlz,
      mcnlvx,
      mcnlvy,
      mcnlvz,
      mcnlt,
      mcnlsx,
      mcnlsy,
      mcnlsz,
      mcnlp); }
#undef mcabsorbComp
#undef p
#undef sz
#undef sy
#undef sx
#undef t
#undef vz
#undef vy
#undef vx
#undef z
#undef y
#undef x
  mcDEBUG_STATE(
mcnlx,
mcnly,
mcnlz,
mcnlvx,
mcnlvy,
mcnlvz,
mcnlt,
mcnlsx,
mcnlsy,
mcnlsz,
mcnlp)

  /* TRACE Component source_div [3] */
  mccoordschange(mcposrsource_div, mcrotrsource_div,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component source_div (without coords transformations) */
  mcJumpTrace_source_div:
  SIG_MESSAGE("source_div (Trace)");
  mcDEBUG_COMP("source_div")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp

#define mcabsorbComp mcabsorbCompsource_div
  STORE_NEUTRON(3,
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp);
  mcScattered=0;
  mcRestore=0;
  mcNCounter[3]++;
  mcPCounter[3] += p;
  mcP2Counter[3] += p*p;
#define mccompcurname  source_div
#define mccompcurtype  Source_div
#define mccompcurindex 3
#define thetah mccsource_div_thetah
#define thetav mccsource_div_thetav
#define sigmah mccsource_div_sigmah
#define sigmav mccsource_div_sigmav
#define tan_h mccsource_div_tan_h
#define tan_v mccsource_div_tan_v
#define p_init mccsource_div_p_init
#define dist mccsource_div_dist
#define focus_xw mccsource_div_focus_xw
#define focus_yh mccsource_div_focus_yh
{   /* Declarations of source_div=Source_div() SETTING parameters. */
MCNUM xwidth = mccsource_div_xwidth;
MCNUM yheight = mccsource_div_yheight;
MCNUM focus_aw = mccsource_div_focus_aw;
MCNUM focus_ah = mccsource_div_focus_ah;
MCNUM E0 = mccsource_div_E0;
MCNUM dE = mccsource_div_dE;
MCNUM lambda0 = mccsource_div_lambda0;
MCNUM dlambda = mccsource_div_dlambda;
MCNUM gauss = mccsource_div_gauss;
MCNUM flux = mccsource_div_flux;
#line 118 "/usr/share/mcstas/2.6.1/tools/Python/mcrun/../mccodelib/../../../sources/Source_div.comp"
{
  double E,lambda,v;

  p=p_init;
  z=0;
  t=0;

  x=randpm1()*xwidth/2.0;
  y=randpm1()*yheight/2.0;
  if(lambda0==0) {
    if (!gauss) {
      E=E0+dE*randpm1();              /*  Choose from uniform distribution */
    } else {
      E=E0+randnorm()*dE;
    }
    v=sqrt(E)*SE2V;
  } else {
    if (!gauss) {
      lambda=lambda0+dlambda*randpm1();
    } else {
      lambda=lambda0+randnorm()*dlambda;
    }
    v = K2V*(2*PI/lambda);
  }

  if (gauss==1) {
      thetah = randnorm()*sigmah;
      thetav = randnorm()*sigmav;
  } else {
      /*find limits of uniform sampling scheme for vertical divergence.
        thetav should be acos(1-2*U) for U\in[0,1]. for theta measured from vertical axis
        we only use a sub-interval for U and measure from horizontal plane.*/
      double sample_lim1,u2;
      sample_lim1=(1-cos(M_PI_2 - focus_ah/2.0*DEG2RAD))*0.5;
      u2=randpm1()*(sample_lim1-0.5) + 0.5;
      thetav = acos(1-2*u2) - M_PI_2;
      thetah = randpm1()*focus_aw*DEG2RAD/2;
  }

  tan_h = tan(thetah);
  tan_v = tan(thetav);

  /* Perform the correct treatment - no small angle approx. here! */
  vz = v / sqrt(1 + tan_v*tan_v + tan_h*tan_h);
  vy = tan_v * vz;
  vx = tan_h * vz;
}
#line 13202 "./ess_extraction.c"
}   /* End of source_div=Source_div() SETTING parameter declarations. */
#undef focus_yh
#undef focus_xw
#undef dist
#undef p_init
#undef tan_v
#undef tan_h
#undef sigmav
#undef sigmah
#undef thetav
#undef thetah
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompsource_div:
  if (RESTORE) /* restore if needed */
  { RESTORE_NEUTRON(3,
      mcnlx,
      mcnly,
      mcnlz,
      mcnlvx,
      mcnlvy,
      mcnlvz,
      mcnlt,
      mcnlsx,
      mcnlsy,
      mcnlsz,
      mcnlp); }
#undef mcabsorbComp
#undef p
#undef sz
#undef sy
#undef sx
#undef t
#undef vz
#undef vy
#undef vx
#undef z
#undef y
#undef x
  mcDEBUG_STATE(
mcnlx,
mcnly,
mcnlz,
mcnlvx,
mcnlvy,
mcnlvz,
mcnlt,
mcnlsx,
mcnlsy,
mcnlsz,
mcnlp)

  /* TRACE Component psd_monitor_source [4] */
  mccoordschange(mcposrpsd_monitor_source, mcrotrpsd_monitor_source,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component psd_monitor_source (without coords transformations) */
  mcJumpTrace_psd_monitor_source:
  SIG_MESSAGE("psd_monitor_source (Trace)");
  mcDEBUG_COMP("psd_monitor_source")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp

#define mcabsorbComp mcabsorbComppsd_monitor_source
  STORE_NEUTRON(4,
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp);
  mcScattered=0;
  mcRestore=0;
  mcNCounter[4]++;
  mcPCounter[4] += p;
  mcP2Counter[4] += p*p;
#define mccompcurname  psd_monitor_source
#define mccompcurtype  PSD_monitor
#define mccompcurindex 4
#define PSD_N mccpsd_monitor_source_PSD_N
#define PSD_p mccpsd_monitor_source_PSD_p
#define PSD_p2 mccpsd_monitor_source_PSD_p2
{   /* Declarations of psd_monitor_source=PSD_monitor() SETTING parameters. */
int nx = mccpsd_monitor_source_nx;
int ny = mccpsd_monitor_source_ny;
char* filename = mccpsd_monitor_source_filename;
MCNUM xmin = mccpsd_monitor_source_xmin;
MCNUM xmax = mccpsd_monitor_source_xmax;
MCNUM ymin = mccpsd_monitor_source_ymin;
MCNUM ymax = mccpsd_monitor_source_ymax;
MCNUM xwidth = mccpsd_monitor_source_xwidth;
MCNUM yheight = mccpsd_monitor_source_yheight;
MCNUM restore_neutron = mccpsd_monitor_source_restore_neutron;
#line 94 "/usr/share/mcstas/2.6.1/tools/Python/mcrun/../mccodelib/../../../monitors/PSD_monitor.comp"
{
  PROP_Z0;
  if (x>xmin && x<xmax && y>ymin && y<ymax){
    int i = floor((x - xmin)*nx/(xmax - xmin));
    int j = floor((y - ymin)*ny/(ymax - ymin));
    PSD_N[i][j]++;
    PSD_p[i][j] += p;
    PSD_p2[i][j] += p*p;
    SCATTER;
  }
  if (restore_neutron) {
    RESTORE_NEUTRON(INDEX_CURRENT_COMP, x, y, z, vx, vy, vz, t, sx, sy, sz, p);
  }
}
#line 13346 "./ess_extraction.c"
}   /* End of psd_monitor_source=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbComppsd_monitor_source:
  if (RESTORE) /* restore if needed */
  { RESTORE_NEUTRON(4,
      mcnlx,
      mcnly,
      mcnlz,
      mcnlvx,
      mcnlvy,
      mcnlvz,
      mcnlt,
      mcnlsx,
      mcnlsy,
      mcnlsz,
      mcnlp); }
#undef mcabsorbComp
#undef p
#undef sz
#undef sy
#undef sx
#undef t
#undef vz
#undef vy
#undef vx
#undef z
#undef y
#undef x
  mcDEBUG_STATE(
mcnlx,
mcnly,
mcnlz,
mcnlvx,
mcnlvy,
mcnlvz,
mcnlt,
mcnlsx,
mcnlsy,
mcnlsz,
mcnlp)

  /* TRACE Component parabolic_optic_before_guide_v [5] */
  mccoordschange(mcposrparabolic_optic_before_guide_v, mcrotrparabolic_optic_before_guide_v,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component parabolic_optic_before_guide_v (without coords transformations) */
  mcJumpTrace_parabolic_optic_before_guide_v:
  SIG_MESSAGE("parabolic_optic_before_guide_v (Trace)");
  mcDEBUG_COMP("parabolic_optic_before_guide_v")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp

#define mcabsorbComp mcabsorbCompparabolic_optic_before_guide_v
  STORE_NEUTRON(5,
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp);
  mcScattered=0;
  mcRestore=0;
  mcNCounter[5]++;
  mcPCounter[5] += p;
  mcP2Counter[5] += p*p;
#define mccompcurname  parabolic_optic_before_guide_v
#define mccompcurtype  FlatEllipse_finite_mirror
#define mccompcurindex 5
#define reflect mccparabolic_optic_before_guide_v_reflect
#define s mccparabolic_optic_before_guide_v_s
#define pTable mccparabolic_optic_before_guide_v_pTable
#define R0 mccparabolic_optic_before_guide_v_R0
#define Qc mccparabolic_optic_before_guide_v_Qc
#define W mccparabolic_optic_before_guide_v_W
#define alpha mccparabolic_optic_before_guide_v_alpha
#define transmit mccparabolic_optic_before_guide_v_transmit
{   /* Declarations of parabolic_optic_before_guide_v=FlatEllipse_finite_mirror() SETTING parameters. */
MCNUM sourceDist = mccparabolic_optic_before_guide_v_sourceDist;
MCNUM LStart = mccparabolic_optic_before_guide_v_LStart;
MCNUM LEnd = mccparabolic_optic_before_guide_v_LEnd;
MCNUM lStart = mccparabolic_optic_before_guide_v_lStart;
MCNUM lEnd = mccparabolic_optic_before_guide_v_lEnd;
MCNUM r_0 = mccparabolic_optic_before_guide_v_r_0;
MCNUM nummirror = mccparabolic_optic_before_guide_v_nummirror;
MCNUM mf = mccparabolic_optic_before_guide_v_mf;
MCNUM mb = mccparabolic_optic_before_guide_v_mb;
MCNUM mirror_width = mccparabolic_optic_before_guide_v_mirror_width;
MCNUM doubleReflections = mccparabolic_optic_before_guide_v_doubleReflections;
/* 'parabolic_optic_before_guide_v=FlatEllipse_finite_mirror()' component instance has conditional execution */
if (( mcipplaceholder > 0 ))

#line 163 "FlatEllipse_finite_mirror.comp"
{
    dt = (-z + *pointer_lStart)/vz; // first propagate neutron to the entrance window
    if (dt < 0) {
        printf("negative time\n");
    }
    PROP_DT(dt); //propagate neutron to the entrance window of the NMO
    Particle pa = makeParticle(x, y, z, vx, vy, vz, t, sx, sy, sz, silicon, p); //Assume the particle is not arriving in silicon
    if (mirror_width>0){ // if the width of the mirrors is finite neutrons have to know whether they are in silicon or not
        for (int i = 0; i < nummirror; i++){
            dt = sqrt(rs[i]*rs[i]); //make sure the mirror distance to check against is positive
            if (dt +mirror_width >= fabs(x)){ //backside of the mirror further out than neutron
                if (dt <= sqrt(x*x)) { // mirror itself closer to the optical axis than the mirro, i.e., we arrive in silicon
                    //ABSORB;
                    pa = makeParticle(x,y,z,vx,vy,vz, t, sx, sy, sz, -silicon, p); //create a particle knowing it is in silicon
                    break;
                    }
                }
            else{   // we do not arrive in Silicon
                //printf("arrived in Air\n");
                    break;
                    }

        }
    }

    traceSingleNeutron(&pa,s);//trace the neutron through the mirror assembly



    //Communicate particle state to McStas
    x = pa._x;
    y = pa._y;
    z = pa._z;
    vx = pa._vx;
    vy = pa._vy;
    vz = pa._vz;
    t = pa._t;
    sx = pa._sx;
    sy = pa._sy;
    sz = pa._sz;
    p = pa.w;

    if (pa.absorb)
        ABSORB;

    SCATTER;
}
#line 13524 "./ess_extraction.c"
}   /* End of parabolic_optic_before_guide_v=FlatEllipse_finite_mirror() SETTING parameter declarations. */
#undef transmit
#undef alpha
#undef W
#undef Qc
#undef R0
#undef pTable
#undef s
#undef reflect
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompparabolic_optic_before_guide_v:
  if (RESTORE) /* restore if needed */
  { RESTORE_NEUTRON(5,
      mcnlx,
      mcnly,
      mcnlz,
      mcnlvx,
      mcnlvy,
      mcnlvz,
      mcnlt,
      mcnlsx,
      mcnlsy,
      mcnlsz,
      mcnlp); }
#undef mcabsorbComp
#undef p
#undef sz
#undef sy
#undef sx
#undef t
#undef vz
#undef vy
#undef vx
#undef z
#undef y
#undef x
  mcDEBUG_STATE(
mcnlx,
mcnly,
mcnlz,
mcnlvx,
mcnlvy,
mcnlvz,
mcnlt,
mcnlsx,
mcnlsy,
mcnlsz,
mcnlp)

  /* TRACE Component parabolic_optic_before_guide_h [6] */
  mccoordschange(mcposrparabolic_optic_before_guide_h, mcrotrparabolic_optic_before_guide_h,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component parabolic_optic_before_guide_h (without coords transformations) */
  mcJumpTrace_parabolic_optic_before_guide_h:
  SIG_MESSAGE("parabolic_optic_before_guide_h (Trace)");
  mcDEBUG_COMP("parabolic_optic_before_guide_h")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp

#define mcabsorbComp mcabsorbCompparabolic_optic_before_guide_h
  STORE_NEUTRON(6,
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp);
  mcScattered=0;
  mcRestore=0;
  mcNCounter[6]++;
  mcPCounter[6] += p;
  mcP2Counter[6] += p*p;
#define mccompcurname  parabolic_optic_before_guide_h
#define mccompcurtype  FlatEllipse_finite_mirror
#define mccompcurindex 6
#define reflect mccparabolic_optic_before_guide_h_reflect
#define s mccparabolic_optic_before_guide_h_s
#define pTable mccparabolic_optic_before_guide_h_pTable
#define R0 mccparabolic_optic_before_guide_h_R0
#define Qc mccparabolic_optic_before_guide_h_Qc
#define W mccparabolic_optic_before_guide_h_W
#define alpha mccparabolic_optic_before_guide_h_alpha
#define transmit mccparabolic_optic_before_guide_h_transmit
{   /* Declarations of parabolic_optic_before_guide_h=FlatEllipse_finite_mirror() SETTING parameters. */
MCNUM sourceDist = mccparabolic_optic_before_guide_h_sourceDist;
MCNUM LStart = mccparabolic_optic_before_guide_h_LStart;
MCNUM LEnd = mccparabolic_optic_before_guide_h_LEnd;
MCNUM lStart = mccparabolic_optic_before_guide_h_lStart;
MCNUM lEnd = mccparabolic_optic_before_guide_h_lEnd;
MCNUM r_0 = mccparabolic_optic_before_guide_h_r_0;
MCNUM nummirror = mccparabolic_optic_before_guide_h_nummirror;
MCNUM mf = mccparabolic_optic_before_guide_h_mf;
MCNUM mb = mccparabolic_optic_before_guide_h_mb;
MCNUM mirror_width = mccparabolic_optic_before_guide_h_mirror_width;
MCNUM doubleReflections = mccparabolic_optic_before_guide_h_doubleReflections;
/* 'parabolic_optic_before_guide_h=FlatEllipse_finite_mirror()' component instance has conditional execution */
if (( mcipplaceholder > 0 ))

#line 163 "FlatEllipse_finite_mirror.comp"
{
    dt = (-z + *pointer_lStart)/vz; // first propagate neutron to the entrance window
    if (dt < 0) {
        printf("negative time\n");
    }
    PROP_DT(dt); //propagate neutron to the entrance window of the NMO
    Particle pa = makeParticle(x, y, z, vx, vy, vz, t, sx, sy, sz, silicon, p); //Assume the particle is not arriving in silicon
    if (mirror_width>0){ // if the width of the mirrors is finite neutrons have to know whether they are in silicon or not
        for (int i = 0; i < nummirror; i++){
            dt = sqrt(rs[i]*rs[i]); //make sure the mirror distance to check against is positive
            if (dt +mirror_width >= fabs(x)){ //backside of the mirror further out than neutron
                if (dt <= sqrt(x*x)) { // mirror itself closer to the optical axis than the mirro, i.e., we arrive in silicon
                    //ABSORB;
                    pa = makeParticle(x,y,z,vx,vy,vz, t, sx, sy, sz, -silicon, p); //create a particle knowing it is in silicon
                    break;
                    }
                }
            else{   // we do not arrive in Silicon
                //printf("arrived in Air\n");
                    break;
                    }

        }
    }

    traceSingleNeutron(&pa,s);//trace the neutron through the mirror assembly



    //Communicate particle state to McStas
    x = pa._x;
    y = pa._y;
    z = pa._z;
    vx = pa._vx;
    vy = pa._vy;
    vz = pa._vz;
    t = pa._t;
    sx = pa._sx;
    sy = pa._sy;
    sz = pa._sz;
    p = pa.w;

    if (pa.absorb)
        ABSORB;

    SCATTER;
}
#line 13707 "./ess_extraction.c"
}   /* End of parabolic_optic_before_guide_h=FlatEllipse_finite_mirror() SETTING parameter declarations. */
#undef transmit
#undef alpha
#undef W
#undef Qc
#undef R0
#undef pTable
#undef s
#undef reflect
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompparabolic_optic_before_guide_h:
  if (RESTORE) /* restore if needed */
  { RESTORE_NEUTRON(6,
      mcnlx,
      mcnly,
      mcnlz,
      mcnlvx,
      mcnlvy,
      mcnlvz,
      mcnlt,
      mcnlsx,
      mcnlsy,
      mcnlsz,
      mcnlp); }
#undef mcabsorbComp
#undef p
#undef sz
#undef sy
#undef sx
#undef t
#undef vz
#undef vy
#undef vx
#undef z
#undef y
#undef x
  mcDEBUG_STATE(
mcnlx,
mcnly,
mcnlz,
mcnlvx,
mcnlvy,
mcnlvz,
mcnlt,
mcnlsx,
mcnlsy,
mcnlsz,
mcnlp)

  /* TRACE Component after_optic_source [7] */
  mccoordschange(mcposrafter_optic_source, mcrotrafter_optic_source,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component after_optic_source (without coords transformations) */
  mcJumpTrace_after_optic_source:
  SIG_MESSAGE("after_optic_source (Trace)");
  mcDEBUG_COMP("after_optic_source")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp

#define mcabsorbComp mcabsorbCompafter_optic_source
  STORE_NEUTRON(7,
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp);
  mcScattered=0;
  mcRestore=0;
  mcNCounter[7]++;
  mcPCounter[7] += p;
  mcP2Counter[7] += p*p;
#define mccompcurname  after_optic_source
#define mccompcurtype  Arm
#define mccompcurindex 7
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompafter_optic_source:
  if (RESTORE) /* restore if needed */
  { RESTORE_NEUTRON(7,
      mcnlx,
      mcnly,
      mcnlz,
      mcnlvx,
      mcnlvy,
      mcnlvz,
      mcnlt,
      mcnlsx,
      mcnlsy,
      mcnlsz,
      mcnlp); }
#undef mcabsorbComp
#undef p
#undef sz
#undef sy
#undef sx
#undef t
#undef vz
#undef vy
#undef vx
#undef z
#undef y
#undef x
  mcDEBUG_STATE(
mcnlx,
mcnly,
mcnlz,
mcnlvx,
mcnlvy,
mcnlvz,
mcnlt,
mcnlsx,
mcnlsy,
mcnlsz,
mcnlp)

  /* TRACE Component psd_monitor_afteropticsource [8] */
  mccoordschange(mcposrpsd_monitor_afteropticsource, mcrotrpsd_monitor_afteropticsource,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component psd_monitor_afteropticsource (without coords transformations) */
  mcJumpTrace_psd_monitor_afteropticsource:
  SIG_MESSAGE("psd_monitor_afteropticsource (Trace)");
  mcDEBUG_COMP("psd_monitor_afteropticsource")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp

#define mcabsorbComp mcabsorbComppsd_monitor_afteropticsource
  STORE_NEUTRON(8,
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp);
  mcScattered=0;
  mcRestore=0;
  mcNCounter[8]++;
  mcPCounter[8] += p;
  mcP2Counter[8] += p*p;
#define mccompcurname  psd_monitor_afteropticsource
#define mccompcurtype  PSD_monitor
#define mccompcurindex 8
#define PSD_N mccpsd_monitor_afteropticsource_PSD_N
#define PSD_p mccpsd_monitor_afteropticsource_PSD_p
#define PSD_p2 mccpsd_monitor_afteropticsource_PSD_p2
{   /* Declarations of psd_monitor_afteropticsource=PSD_monitor() SETTING parameters. */
int nx = mccpsd_monitor_afteropticsource_nx;
int ny = mccpsd_monitor_afteropticsource_ny;
char* filename = mccpsd_monitor_afteropticsource_filename;
MCNUM xmin = mccpsd_monitor_afteropticsource_xmin;
MCNUM xmax = mccpsd_monitor_afteropticsource_xmax;
MCNUM ymin = mccpsd_monitor_afteropticsource_ymin;
MCNUM ymax = mccpsd_monitor_afteropticsource_ymax;
MCNUM xwidth = mccpsd_monitor_afteropticsource_xwidth;
MCNUM yheight = mccpsd_monitor_afteropticsource_yheight;
MCNUM restore_neutron = mccpsd_monitor_afteropticsource_restore_neutron;
#line 94 "/usr/share/mcstas/2.6.1/tools/Python/mcrun/../mccodelib/../../../monitors/PSD_monitor.comp"
{
  PROP_Z0;
  if (x>xmin && x<xmax && y>ymin && y<ymax){
    int i = floor((x - xmin)*nx/(xmax - xmin));
    int j = floor((y - ymin)*ny/(ymax - ymin));
    PSD_N[i][j]++;
    PSD_p[i][j] += p;
    PSD_p2[i][j] += p*p;
    SCATTER;
  }
  if (restore_neutron) {
    RESTORE_NEUTRON(INDEX_CURRENT_COMP, x, y, z, vx, vy, vz, t, sx, sy, sz, p);
  }
}
#line 13952 "./ess_extraction.c"
}   /* End of psd_monitor_afteropticsource=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbComppsd_monitor_afteropticsource:
  if (RESTORE) /* restore if needed */
  { RESTORE_NEUTRON(8,
      mcnlx,
      mcnly,
      mcnlz,
      mcnlvx,
      mcnlvy,
      mcnlvz,
      mcnlt,
      mcnlsx,
      mcnlsy,
      mcnlsz,
      mcnlp); }
#undef mcabsorbComp
#undef p
#undef sz
#undef sy
#undef sx
#undef t
#undef vz
#undef vy
#undef vx
#undef z
#undef y
#undef x
  mcDEBUG_STATE(
mcnlx,
mcnly,
mcnlz,
mcnlvx,
mcnlvy,
mcnlvz,
mcnlt,
mcnlsx,
mcnlsy,
mcnlsz,
mcnlp)

  /* TRACE Component guide_gravity_1 [9] */
  mccoordschange(mcposrguide_gravity_1, mcrotrguide_gravity_1,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component guide_gravity_1 (without coords transformations) */
  mcJumpTrace_guide_gravity_1:
  SIG_MESSAGE("guide_gravity_1 (Trace)");
  mcDEBUG_COMP("guide_gravity_1")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp

#define mcabsorbComp mcabsorbCompguide_gravity_1
  STORE_NEUTRON(9,
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp);
  mcScattered=0;
  mcRestore=0;
  mcNCounter[9]++;
  mcPCounter[9] += p;
  mcP2Counter[9] += p*p;
#define mccompcurname  guide_gravity_1
#define mccompcurtype  Guide_gravity
#define mccompcurindex 9
#define GVars mccguide_gravity_1_GVars
#define pTable mccguide_gravity_1_pTable
{   /* Declarations of guide_gravity_1=Guide_gravity() SETTING parameters. */
MCNUM w1 = mccguide_gravity_1_w1;
MCNUM h1 = mccguide_gravity_1_h1;
MCNUM w2 = mccguide_gravity_1_w2;
MCNUM h2 = mccguide_gravity_1_h2;
MCNUM l = mccguide_gravity_1_l;
MCNUM R0 = mccguide_gravity_1_R0;
MCNUM Qc = mccguide_gravity_1_Qc;
MCNUM alpha = mccguide_gravity_1_alpha;
MCNUM m = mccguide_gravity_1_m;
MCNUM W = mccguide_gravity_1_W;
MCNUM nslit = mccguide_gravity_1_nslit;
MCNUM d = mccguide_gravity_1_d;
MCNUM mleft = mccguide_gravity_1_mleft;
MCNUM mright = mccguide_gravity_1_mright;
MCNUM mtop = mccguide_gravity_1_mtop;
MCNUM mbottom = mccguide_gravity_1_mbottom;
MCNUM nhslit = mccguide_gravity_1_nhslit;
MCNUM G = mccguide_gravity_1_G;
MCNUM aleft = mccguide_gravity_1_aleft;
MCNUM aright = mccguide_gravity_1_aright;
MCNUM atop = mccguide_gravity_1_atop;
MCNUM abottom = mccguide_gravity_1_abottom;
MCNUM wavy = mccguide_gravity_1_wavy;
MCNUM wavy_z = mccguide_gravity_1_wavy_z;
MCNUM wavy_tb = mccguide_gravity_1_wavy_tb;
MCNUM wavy_lr = mccguide_gravity_1_wavy_lr;
MCNUM chamfers = mccguide_gravity_1_chamfers;
MCNUM chamfers_z = mccguide_gravity_1_chamfers_z;
MCNUM chamfers_lr = mccguide_gravity_1_chamfers_lr;
MCNUM chamfers_tb = mccguide_gravity_1_chamfers_tb;
MCNUM nelements = mccguide_gravity_1_nelements;
MCNUM nu = mccguide_gravity_1_nu;
MCNUM phase = mccguide_gravity_1_phase;
char* reflect = mccguide_gravity_1_reflect;
#line 392 "/usr/share/mcstas/2.6.1/tools/Python/mcrun/../mccodelib/../../../optics/Guide_gravity.comp"
{
  if (l > 0 && nelements > 0) {
    double B, C, dt;
    int    ret, bounces = 0, i=0;
    double this_width, this_height;
    double angle=0;

    if (GVars.fc_freq != 0 || GVars.fc_phase != 0) { /* rotate neutron w/r to guide element */
      /* approximation of rotating straight Fermi Chopper */
      Coords   X = coords_set(x,y,z-l/2);  /* current coordinates of neutron in centered static frame */
      Rotation R;
      double dt=(-z+l/2)/vz; /* time shift to each center of slit package */
      angle=fmod(360*GVars.fc_freq*(t+dt)+GVars.fc_phase, 360); /* in deg */
      /* modify angle so that Z0 guide side is always in front of incoming neutron */
      if (angle > 90 && angle < 270) { angle -= 180; }
      angle *= DEG2RAD;
      rot_set_rotation(R, 0, -angle, 0); /* will rotate neutron instead of comp: negative side */
      /* apply rotation to centered coordinates */
      Coords   RX = rot_apply(R, X);
      coords_get(RX, &x, &y, &z);
      z = z+l/2;
      /* rotate speed */
      X  = coords_set(vx,vy,vz);
      RX = rot_apply(R, X);
      coords_get(RX, &vx, &vy, &vz);
    }

    for (i=0; i<7; GVars.N_reflection[i++] = 0);

    /* propagate to box input (with gravitation) in comp local coords */
    /* A = 0.5 n.g; B = n.v; C = n.(r-W); */
    /* 0=Z0 side: n=(0, 0, -l) ; W = (0, 0, 0) (at z=0, guide input)*/
    B = -l*vz; C = -l*z;

    ret = solve_2nd_order(&dt, NULL, GVars.A[0], B, C);
    if (ret==0) ABSORB;

    if (dt>0.0) PROP_GRAV_DT(dt, GVars.gx, GVars.gy, GVars.gz); else if (angle) ABSORB;
    GVars.N_reflection[6]++;

    this_width  = w1;
    this_height = h1;

  /* check if we are in the box input, else absorb */
    if (fabs(x) > this_width/2 || fabs(y) > this_height/2)
      ABSORB;
    else
    {
      double w_edge, w_adj; /* Channel displacement on X */
      double h_edge, h_adj; /* Channel displacement on Y */
      double w_chnum,h_chnum; /* channel indexes */

      SCATTER;

      /* X: Shift origin to center of channel hit (absorb if hit dividing walls) */
      x += w1/2.0;
      w_chnum = floor(x/(GVars.w1c+d));  /* 0= right side, nslit+1=left side  */
      w_edge  = w_chnum*(GVars.w1c+d);
      if(x - w_edge > GVars.w1c)
      {
        x -= w1/2.0; /* Re-adjust origin */
        ABSORB;
      }
      w_adj = w_edge + (GVars.w1c)/2.0;
      x -= w_adj; w_adj -=  w1/2.0;

      /* Y: Shift origin to center of channel hit (absorb if hit dividing walls) */
      y += h1/2.0;
      h_chnum = floor(y/(GVars.h1c+d));  /* 0= lower side, nslit+1=upper side  */
      h_edge  = h_chnum*(GVars.h1c+d);
      if(y - h_edge > GVars.h1c)
      {
        y -= h1/2.0; /* Re-adjust origin */
        ABSORB;
      }
      h_adj = h_edge + (GVars.h1c)/2.0;
      y -= h_adj; h_adj -=  h1/2.0;

      /* neutron is now in the input window of the guide */
      /* do loops on reflections in the box */
      for(;;)
      {
        /* get intersections for all box sides */
        double q, nx,ny,nz;
        double this_length;
        int side=0;

        bounces++;
        /* now look for intersection with guide sides and exit */
        side = Gravity_guide_Trace(&dt, &GVars, x, y, z,
            vx, vy, vz, w_chnum, nslit, h_chnum, nhslit,
            &nx, &ny, &nz);

        /* only positive dt are valid */
        /* exit reflection loops if no intersection (neutron is after box) */
        if (side == 0 || dt <= 0)
          { if (GVars.warnings < 100)
              fprintf(stderr,"%s: warning: neutron has entered guide, but can not exit !\n", GVars.compcurname);
            GVars.warnings++;
            x += w_adj; y += h_adj; ABSORB; } /* should never occur */

        /* propagate to dt */
        PROP_GRAV_DT(dt, GVars.gx, GVars.gy, GVars.gz);

        /* do reflection on speed for l/r/u/d sides */
        if (side == 5) /* neutron reaches end of guide: end loop and exit comp */
          { GVars.N_reflection[side]++; x += w_adj; y += h_adj; SCATTER; x -= w_adj; y -= h_adj; break; }
        /* else reflection on a guide wall */
        if(GVars.M[side] == 0 || Qc == 0 || R0 == 0)  /* walls are absorbing */
          { x += w_adj; y += h_adj; ABSORB; }
        /* handle chamfers */
        this_width = w1+(w2-w1)*z/l;
        this_height= h1+(h2-h1)*z/l;
        this_length= fmod(z, l/nelements);
        /* absorb on input/output of element parts */
        if (GVars.chamfer_z && (this_length<GVars.chamfer_z || this_length>l/nelements-GVars.chamfer_z))
        { x += w_adj; y += h_adj; ABSORB; }
        /* absorb on l/r/t/b sides */
        if (GVars.chamfer_lr && (side==1 || side==2) && (fabs(y+h_adj)>this_height/2-GVars.chamfer_lr))
        { x += w_adj; y += h_adj; ABSORB; }
        if (GVars.chamfer_tb && (side==3 || side==4) && (fabs(x+w_adj)>this_width/2- GVars.chamfer_tb))
        { x += w_adj; y += h_adj; ABSORB; }
        /* change/mirror velocity: h_f = v - n.2*n.v/|n|^2 */
        GVars.N_reflection[side]++; /* GVars.norm_n2 > 0 was checked at INIT */
        /* compute n.v using current values */
        B = scalar_prod(vx,vy,vz,nx,ny,nz);
        dt = 2*B/GVars.norm_n2[side]; /* 2*n.v/|n|^2 */
        vx -= nx*dt;
        vy -= ny*dt;
        vz -= nz*dt;

        /* compute q and modify neutron weight */
        /* scattering q=|n_i-n_f| = V2Q*|vf - v| = V2Q*2*n.v/|n| */
        q = 2*V2Q*fabs(B)/GVars.norm_n[side];

        if (reflect && strlen(reflect) && strcmp(reflect,"NULL") && strcmp(reflect,"0"))
          TableReflecFunc(q, &pTable, &B);
        else {
          double par[] = {R0, Qc, GVars.Alpha[side], GVars.M[side], W};
          StdReflecFunc(q, par, &B);
        }
        if (B <= 0) { x += w_adj; y += h_adj; ABSORB; }
        else p *= B;
        x += w_adj; y += h_adj; SCATTER; x -= w_adj; y -= h_adj;
        GVars.N_reflection[0]++;
        /* go to the next reflection */
        if (bounces > 1000) ABSORB;
      } /* end for */
      x += w_adj; y += h_adj; /* Re-adjust origin after SCATTER */
    }

    if (GVars.fc_freq != 0 || GVars.fc_phase != 0) { /* rotate back neutron w/r to guide element */
      /* approximation of rotating straight Fermi Chopper */
      Coords   X = coords_set(x,y,z-l/2);  /* current coordinates of neutron in centered static frame */
      Rotation R;
      rot_set_rotation(R, 0, angle, 0); /* will rotate back neutron: positive side */
      /* apply rotation to centered coordinates */
      Coords   RX = rot_apply(R, X);
      coords_get(RX, &x, &y, &z);
      z = z+l/2;
      /* rotate speed */
      X  = coords_set(vx,vy,vz);
      RX = rot_apply(R, X);
      coords_get(RX, &vx, &vy, &vz);
    }

  } /* if l */
}
#line 14266 "./ess_extraction.c"
}   /* End of guide_gravity_1=Guide_gravity() SETTING parameter declarations. */
#undef pTable
#undef GVars
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompguide_gravity_1:
  if (RESTORE) /* restore if needed */
  { RESTORE_NEUTRON(9,
      mcnlx,
      mcnly,
      mcnlz,
      mcnlvx,
      mcnlvy,
      mcnlvz,
      mcnlt,
      mcnlsx,
      mcnlsy,
      mcnlsz,
      mcnlp); }
#undef mcabsorbComp
#undef p
#undef sz
#undef sy
#undef sx
#undef t
#undef vz
#undef vy
#undef vx
#undef z
#undef y
#undef x
  mcDEBUG_STATE(
mcnlx,
mcnly,
mcnlz,
mcnlvx,
mcnlvy,
mcnlvz,
mcnlt,
mcnlsx,
mcnlsy,
mcnlsz,
mcnlp)

  /* TRACE Component monitor_1 [10] */
  mccoordschange(mcposrmonitor_1, mcrotrmonitor_1,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component monitor_1 (without coords transformations) */
  mcJumpTrace_monitor_1:
  SIG_MESSAGE("monitor_1 (Trace)");
  mcDEBUG_COMP("monitor_1")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp

#define mcabsorbComp mcabsorbCompmonitor_1
  STORE_NEUTRON(10,
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp);
  mcScattered=0;
  mcRestore=0;
  mcNCounter[10]++;
  mcPCounter[10] += p;
  mcP2Counter[10] += p*p;
#define mccompcurname  monitor_1
#define mccompcurtype  Arm
#define mccompcurindex 10
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompmonitor_1:
  if (RESTORE) /* restore if needed */
  { RESTORE_NEUTRON(10,
      mcnlx,
      mcnly,
      mcnlz,
      mcnlvx,
      mcnlvy,
      mcnlvz,
      mcnlt,
      mcnlsx,
      mcnlsy,
      mcnlsz,
      mcnlp); }
#undef mcabsorbComp
#undef p
#undef sz
#undef sy
#undef sx
#undef t
#undef vz
#undef vy
#undef vx
#undef z
#undef y
#undef x
  mcDEBUG_STATE(
mcnlx,
mcnly,
mcnlz,
mcnlvx,
mcnlvy,
mcnlvz,
mcnlt,
mcnlsx,
mcnlsy,
mcnlsz,
mcnlp)

  /* TRACE Component divpos_monitor [11] */
  mccoordschange(mcposrdivpos_monitor, mcrotrdivpos_monitor,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component divpos_monitor (without coords transformations) */
  mcJumpTrace_divpos_monitor:
  SIG_MESSAGE("divpos_monitor (Trace)");
  mcDEBUG_COMP("divpos_monitor")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp

#define mcabsorbComp mcabsorbCompdivpos_monitor
  STORE_NEUTRON(11,
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp);
  mcScattered=0;
  mcRestore=0;
  mcNCounter[11]++;
  mcPCounter[11] += p;
  mcP2Counter[11] += p*p;
#define mccompcurname  divpos_monitor
#define mccompcurtype  DivPos_monitor
#define mccompcurindex 11
#define nh mccdivpos_monitor_nh
#define ndiv mccdivpos_monitor_ndiv
#define Div_N mccdivpos_monitor_Div_N
#define Div_p mccdivpos_monitor_Div_p
#define Div_p2 mccdivpos_monitor_Div_p2
{   /* Declarations of divpos_monitor=DivPos_monitor() SETTING parameters. */
char* filename = mccdivpos_monitor_filename;
MCNUM xmin = mccdivpos_monitor_xmin;
MCNUM xmax = mccdivpos_monitor_xmax;
MCNUM ymin = mccdivpos_monitor_ymin;
MCNUM ymax = mccdivpos_monitor_ymax;
MCNUM xwidth = mccdivpos_monitor_xwidth;
MCNUM yheight = mccdivpos_monitor_yheight;
MCNUM maxdiv_h = mccdivpos_monitor_maxdiv_h;
MCNUM restore_neutron = mccdivpos_monitor_restore_neutron;
MCNUM nx = mccdivpos_monitor_nx;
MCNUM ny = mccdivpos_monitor_ny;
MCNUM nz = mccdivpos_monitor_nz;
int nowritefile = mccdivpos_monitor_nowritefile;
#line 92 "/usr/share/mcstas/2.6.1/tools/Python/mcrun/../mccodelib/../../../monitors/DivPos_monitor.comp"
{
    int i,j;
    double div;
    double v, vn;

    PROP_Z0;
    if (x>xmin && x<xmax && y>ymin && y<ymax)
    {
      /* Find length of projection onto the [nx ny nz] axis */
      vn = scalar_prod(vx, vy, vz, nx, ny, nz);
      div = RAD2DEG*atan2(vx,vn);

      if (div < maxdiv_h && div > -maxdiv_h)
      {
        i = floor((x - xmin)*nh/(xmax - xmin));
        j = floor((div + maxdiv_h)*ndiv/(2.0*maxdiv_h));
        Div_N[i][j]++;
        Div_p[i][j] += p;
        Div_p2[i][j] += p*p;
        SCATTER;
      }
    }
    if (restore_neutron) {
      RESTORE_NEUTRON(INDEX_CURRENT_COMP, x, y, z, vx, vy, vz, t, sx, sy, sz, p);
    }
}
#line 14522 "./ess_extraction.c"
}   /* End of divpos_monitor=DivPos_monitor() SETTING parameter declarations. */
#undef Div_p2
#undef Div_p
#undef Div_N
#undef ndiv
#undef nh
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompdivpos_monitor:
  if (RESTORE) /* restore if needed */
  { RESTORE_NEUTRON(11,
      mcnlx,
      mcnly,
      mcnlz,
      mcnlvx,
      mcnlvy,
      mcnlvz,
      mcnlt,
      mcnlsx,
      mcnlsy,
      mcnlsz,
      mcnlp); }
#undef mcabsorbComp
#undef p
#undef sz
#undef sy
#undef sx
#undef t
#undef vz
#undef vy
#undef vx
#undef z
#undef y
#undef x
  mcDEBUG_STATE(
mcnlx,
mcnly,
mcnlz,
mcnlvx,
mcnlvy,
mcnlvz,
mcnlt,
mcnlsx,
mcnlsy,
mcnlsz,
mcnlp)

  /* TRACE Component psd_monitor_g1 [12] */
  mccoordschange(mcposrpsd_monitor_g1, mcrotrpsd_monitor_g1,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component psd_monitor_g1 (without coords transformations) */
  mcJumpTrace_psd_monitor_g1:
  SIG_MESSAGE("psd_monitor_g1 (Trace)");
  mcDEBUG_COMP("psd_monitor_g1")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp

#define mcabsorbComp mcabsorbComppsd_monitor_g1
  STORE_NEUTRON(12,
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp);
  mcScattered=0;
  mcRestore=0;
  mcNCounter[12]++;
  mcPCounter[12] += p;
  mcP2Counter[12] += p*p;
#define mccompcurname  psd_monitor_g1
#define mccompcurtype  PSD_monitor
#define mccompcurindex 12
#define PSD_N mccpsd_monitor_g1_PSD_N
#define PSD_p mccpsd_monitor_g1_PSD_p
#define PSD_p2 mccpsd_monitor_g1_PSD_p2
{   /* Declarations of psd_monitor_g1=PSD_monitor() SETTING parameters. */
int nx = mccpsd_monitor_g1_nx;
int ny = mccpsd_monitor_g1_ny;
char* filename = mccpsd_monitor_g1_filename;
MCNUM xmin = mccpsd_monitor_g1_xmin;
MCNUM xmax = mccpsd_monitor_g1_xmax;
MCNUM ymin = mccpsd_monitor_g1_ymin;
MCNUM ymax = mccpsd_monitor_g1_ymax;
MCNUM xwidth = mccpsd_monitor_g1_xwidth;
MCNUM yheight = mccpsd_monitor_g1_yheight;
MCNUM restore_neutron = mccpsd_monitor_g1_restore_neutron;
#line 94 "/usr/share/mcstas/2.6.1/tools/Python/mcrun/../mccodelib/../../../monitors/PSD_monitor.comp"
{
  PROP_Z0;
  if (x>xmin && x<xmax && y>ymin && y<ymax){
    int i = floor((x - xmin)*nx/(xmax - xmin));
    int j = floor((y - ymin)*ny/(ymax - ymin));
    PSD_N[i][j]++;
    PSD_p[i][j] += p;
    PSD_p2[i][j] += p*p;
    SCATTER;
  }
  if (restore_neutron) {
    RESTORE_NEUTRON(INDEX_CURRENT_COMP, x, y, z, vx, vy, vz, t, sx, sy, sz, p);
  }
}
#line 14661 "./ess_extraction.c"
}   /* End of psd_monitor_g1=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbComppsd_monitor_g1:
  if (RESTORE) /* restore if needed */
  { RESTORE_NEUTRON(12,
      mcnlx,
      mcnly,
      mcnlz,
      mcnlvx,
      mcnlvy,
      mcnlvz,
      mcnlt,
      mcnlsx,
      mcnlsy,
      mcnlsz,
      mcnlp); }
#undef mcabsorbComp
#undef p
#undef sz
#undef sy
#undef sx
#undef t
#undef vz
#undef vy
#undef vx
#undef z
#undef y
#undef x
  mcDEBUG_STATE(
mcnlx,
mcnly,
mcnlz,
mcnlvx,
mcnlvy,
mcnlvz,
mcnlt,
mcnlsx,
mcnlsy,
mcnlsz,
mcnlp)

  /* TRACE Component divhlambda_monitor_g1 [13] */
  mccoordschange(mcposrdivhlambda_monitor_g1, mcrotrdivhlambda_monitor_g1,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component divhlambda_monitor_g1 (without coords transformations) */
  mcJumpTrace_divhlambda_monitor_g1:
  SIG_MESSAGE("divhlambda_monitor_g1 (Trace)");
  mcDEBUG_COMP("divhlambda_monitor_g1")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp

#define mcabsorbComp mcabsorbCompdivhlambda_monitor_g1
  STORE_NEUTRON(13,
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp);
  mcScattered=0;
  mcRestore=0;
  mcNCounter[13]++;
  mcPCounter[13] += p;
  mcP2Counter[13] += p*p;
#define mccompcurname  divhlambda_monitor_g1
#define mccompcurtype  DivLambda_monitor
#define mccompcurindex 13
#define nL mccdivhlambda_monitor_g1_nL
#define nh mccdivhlambda_monitor_g1_nh
#define Div_N mccdivhlambda_monitor_g1_Div_N
#define Div_p mccdivhlambda_monitor_g1_Div_p
#define Div_p2 mccdivhlambda_monitor_g1_Div_p2
{   /* Declarations of divhlambda_monitor_g1=DivLambda_monitor() SETTING parameters. */
char* filename = mccdivhlambda_monitor_g1_filename;
MCNUM xmin = mccdivhlambda_monitor_g1_xmin;
MCNUM xmax = mccdivhlambda_monitor_g1_xmax;
MCNUM ymin = mccdivhlambda_monitor_g1_ymin;
MCNUM ymax = mccdivhlambda_monitor_g1_ymax;
MCNUM xwidth = mccdivhlambda_monitor_g1_xwidth;
MCNUM yheight = mccdivhlambda_monitor_g1_yheight;
MCNUM maxdiv_h = mccdivhlambda_monitor_g1_maxdiv_h;
MCNUM Lmin = mccdivhlambda_monitor_g1_Lmin;
MCNUM Lmax = mccdivhlambda_monitor_g1_Lmax;
MCNUM restore_neutron = mccdivhlambda_monitor_g1_restore_neutron;
MCNUM nx = mccdivhlambda_monitor_g1_nx;
MCNUM ny = mccdivhlambda_monitor_g1_ny;
MCNUM nz = mccdivhlambda_monitor_g1_nz;
int nowritefile = mccdivhlambda_monitor_g1_nowritefile;
#line 89 "/usr/share/mcstas/2.6.1/tools/Python/mcrun/../mccodelib/../../../monitors/DivLambda_monitor.comp"
{
    int i,j;
    double div;
    double lambda;
    double v, vn;

    PROP_Z0;
    lambda = (2*PI/V2K)/sqrt(vx*vx + vy*vy + vz*vz);
    if (x>xmin && x<xmax && y>ymin && y<ymax &&
        lambda > Lmin && lambda < Lmax)
    {
      /* Find length of projection onto the [nx ny nz] axis */
      vn = scalar_prod(vx, vy, vz, nx, ny, nz);
      div = RAD2DEG*atan2(vx,vn);

      if (div < maxdiv_h && div > -maxdiv_h)
      {
        i = floor((lambda - Lmin)*nL/(Lmax - Lmin));
        j = floor((div + maxdiv_h)*nh/(2.0*maxdiv_h));
        Div_N[i][j]++;
        Div_p[i][j] += p;
        Div_p2[i][j] += p*p;
        SCATTER;
      }
    }
    if (restore_neutron) {
      RESTORE_NEUTRON(INDEX_CURRENT_COMP, x, y, z, vx, vy, vz, t, sx, sy, sz, p);
    }
}
#line 14820 "./ess_extraction.c"
}   /* End of divhlambda_monitor_g1=DivLambda_monitor() SETTING parameter declarations. */
#undef Div_p2
#undef Div_p
#undef Div_N
#undef nh
#undef nL
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompdivhlambda_monitor_g1:
  if (RESTORE) /* restore if needed */
  { RESTORE_NEUTRON(13,
      mcnlx,
      mcnly,
      mcnlz,
      mcnlvx,
      mcnlvy,
      mcnlvz,
      mcnlt,
      mcnlsx,
      mcnlsy,
      mcnlsz,
      mcnlp); }
#undef mcabsorbComp
#undef p
#undef sz
#undef sy
#undef sx
#undef t
#undef vz
#undef vy
#undef vx
#undef z
#undef y
#undef x
  mcDEBUG_STATE(
mcnlx,
mcnly,
mcnlz,
mcnlvx,
mcnlvy,
mcnlvz,
mcnlt,
mcnlsx,
mcnlsy,
mcnlsz,
mcnlp)

  /* TRACE Component parabolic_optic1 [14] */
  mccoordschange(mcposrparabolic_optic1, mcrotrparabolic_optic1,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component parabolic_optic1 (without coords transformations) */
  mcJumpTrace_parabolic_optic1:
  SIG_MESSAGE("parabolic_optic1 (Trace)");
  mcDEBUG_COMP("parabolic_optic1")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp

#define mcabsorbComp mcabsorbCompparabolic_optic1
  STORE_NEUTRON(14,
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp);
  mcScattered=0;
  mcRestore=0;
  mcNCounter[14]++;
  mcPCounter[14] += p;
  mcP2Counter[14] += p*p;
#define mccompcurname  parabolic_optic1
#define mccompcurtype  FlatEllipse_finite_mirror
#define mccompcurindex 14
#define reflect mccparabolic_optic1_reflect
#define s mccparabolic_optic1_s
#define pTable mccparabolic_optic1_pTable
#define R0 mccparabolic_optic1_R0
#define Qc mccparabolic_optic1_Qc
#define W mccparabolic_optic1_W
#define alpha mccparabolic_optic1_alpha
#define transmit mccparabolic_optic1_transmit
{   /* Declarations of parabolic_optic1=FlatEllipse_finite_mirror() SETTING parameters. */
MCNUM sourceDist = mccparabolic_optic1_sourceDist;
MCNUM LStart = mccparabolic_optic1_LStart;
MCNUM LEnd = mccparabolic_optic1_LEnd;
MCNUM lStart = mccparabolic_optic1_lStart;
MCNUM lEnd = mccparabolic_optic1_lEnd;
MCNUM r_0 = mccparabolic_optic1_r_0;
MCNUM nummirror = mccparabolic_optic1_nummirror;
MCNUM mf = mccparabolic_optic1_mf;
MCNUM mb = mccparabolic_optic1_mb;
MCNUM mirror_width = mccparabolic_optic1_mirror_width;
MCNUM doubleReflections = mccparabolic_optic1_doubleReflections;
/* 'parabolic_optic1=FlatEllipse_finite_mirror()' component instance has conditional execution */
if (( mcipplaceholder > 0 ))

#line 163 "FlatEllipse_finite_mirror.comp"
{
    dt = (-z + *pointer_lStart)/vz; // first propagate neutron to the entrance window
    if (dt < 0) {
        printf("negative time\n");
    }
    PROP_DT(dt); //propagate neutron to the entrance window of the NMO
    Particle pa = makeParticle(x, y, z, vx, vy, vz, t, sx, sy, sz, silicon, p); //Assume the particle is not arriving in silicon
    if (mirror_width>0){ // if the width of the mirrors is finite neutrons have to know whether they are in silicon or not
        for (int i = 0; i < nummirror; i++){
            dt = sqrt(rs[i]*rs[i]); //make sure the mirror distance to check against is positive
            if (dt +mirror_width >= fabs(x)){ //backside of the mirror further out than neutron
                if (dt <= sqrt(x*x)) { // mirror itself closer to the optical axis than the mirro, i.e., we arrive in silicon
                    //ABSORB;
                    pa = makeParticle(x,y,z,vx,vy,vz, t, sx, sy, sz, -silicon, p); //create a particle knowing it is in silicon
                    break;
                    }
                }
            else{   // we do not arrive in Silicon
                //printf("arrived in Air\n");
                    break;
                    }

        }
    }

    traceSingleNeutron(&pa,s);//trace the neutron through the mirror assembly



    //Communicate particle state to McStas
    x = pa._x;
    y = pa._y;
    z = pa._z;
    vx = pa._vx;
    vy = pa._vy;
    vz = pa._vz;
    t = pa._t;
    sx = pa._sx;
    sy = pa._sy;
    sz = pa._sz;
    p = pa.w;

    if (pa.absorb)
        ABSORB;

    SCATTER;
}
#line 15000 "./ess_extraction.c"
}   /* End of parabolic_optic1=FlatEllipse_finite_mirror() SETTING parameter declarations. */
#undef transmit
#undef alpha
#undef W
#undef Qc
#undef R0
#undef pTable
#undef s
#undef reflect
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompparabolic_optic1:
  if (RESTORE) /* restore if needed */
  { RESTORE_NEUTRON(14,
      mcnlx,
      mcnly,
      mcnlz,
      mcnlvx,
      mcnlvy,
      mcnlvz,
      mcnlt,
      mcnlsx,
      mcnlsy,
      mcnlsz,
      mcnlp); }
#undef mcabsorbComp
#undef p
#undef sz
#undef sy
#undef sx
#undef t
#undef vz
#undef vy
#undef vx
#undef z
#undef y
#undef x
  mcDEBUG_STATE(
mcnlx,
mcnly,
mcnlz,
mcnlvx,
mcnlvy,
mcnlvz,
mcnlt,
mcnlsx,
mcnlsy,
mcnlsz,
mcnlp)

  /* TRACE Component parabolic_optic2 [15] */
  mccoordschange(mcposrparabolic_optic2, mcrotrparabolic_optic2,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component parabolic_optic2 (without coords transformations) */
  mcJumpTrace_parabolic_optic2:
  SIG_MESSAGE("parabolic_optic2 (Trace)");
  mcDEBUG_COMP("parabolic_optic2")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp

#define mcabsorbComp mcabsorbCompparabolic_optic2
  STORE_NEUTRON(15,
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp);
  mcScattered=0;
  mcRestore=0;
  mcNCounter[15]++;
  mcPCounter[15] += p;
  mcP2Counter[15] += p*p;
#define mccompcurname  parabolic_optic2
#define mccompcurtype  FlatEllipse_finite_mirror
#define mccompcurindex 15
#define reflect mccparabolic_optic2_reflect
#define s mccparabolic_optic2_s
#define pTable mccparabolic_optic2_pTable
#define R0 mccparabolic_optic2_R0
#define Qc mccparabolic_optic2_Qc
#define W mccparabolic_optic2_W
#define alpha mccparabolic_optic2_alpha
#define transmit mccparabolic_optic2_transmit
{   /* Declarations of parabolic_optic2=FlatEllipse_finite_mirror() SETTING parameters. */
MCNUM sourceDist = mccparabolic_optic2_sourceDist;
MCNUM LStart = mccparabolic_optic2_LStart;
MCNUM LEnd = mccparabolic_optic2_LEnd;
MCNUM lStart = mccparabolic_optic2_lStart;
MCNUM lEnd = mccparabolic_optic2_lEnd;
MCNUM r_0 = mccparabolic_optic2_r_0;
MCNUM nummirror = mccparabolic_optic2_nummirror;
MCNUM mf = mccparabolic_optic2_mf;
MCNUM mb = mccparabolic_optic2_mb;
MCNUM mirror_width = mccparabolic_optic2_mirror_width;
MCNUM doubleReflections = mccparabolic_optic2_doubleReflections;
/* 'parabolic_optic2=FlatEllipse_finite_mirror()' component instance has conditional execution */
if (( mcipplaceholder > 0 ))

#line 163 "FlatEllipse_finite_mirror.comp"
{
    dt = (-z + *pointer_lStart)/vz; // first propagate neutron to the entrance window
    if (dt < 0) {
        printf("negative time\n");
    }
    PROP_DT(dt); //propagate neutron to the entrance window of the NMO
    Particle pa = makeParticle(x, y, z, vx, vy, vz, t, sx, sy, sz, silicon, p); //Assume the particle is not arriving in silicon
    if (mirror_width>0){ // if the width of the mirrors is finite neutrons have to know whether they are in silicon or not
        for (int i = 0; i < nummirror; i++){
            dt = sqrt(rs[i]*rs[i]); //make sure the mirror distance to check against is positive
            if (dt +mirror_width >= fabs(x)){ //backside of the mirror further out than neutron
                if (dt <= sqrt(x*x)) { // mirror itself closer to the optical axis than the mirro, i.e., we arrive in silicon
                    //ABSORB;
                    pa = makeParticle(x,y,z,vx,vy,vz, t, sx, sy, sz, -silicon, p); //create a particle knowing it is in silicon
                    break;
                    }
                }
            else{   // we do not arrive in Silicon
                //printf("arrived in Air\n");
                    break;
                    }

        }
    }

    traceSingleNeutron(&pa,s);//trace the neutron through the mirror assembly



    //Communicate particle state to McStas
    x = pa._x;
    y = pa._y;
    z = pa._z;
    vx = pa._vx;
    vy = pa._vy;
    vz = pa._vz;
    t = pa._t;
    sx = pa._sx;
    sy = pa._sy;
    sz = pa._sz;
    p = pa.w;

    if (pa.absorb)
        ABSORB;

    SCATTER;
}
#line 15183 "./ess_extraction.c"
}   /* End of parabolic_optic2=FlatEllipse_finite_mirror() SETTING parameter declarations. */
#undef transmit
#undef alpha
#undef W
#undef Qc
#undef R0
#undef pTable
#undef s
#undef reflect
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompparabolic_optic2:
  if (RESTORE) /* restore if needed */
  { RESTORE_NEUTRON(15,
      mcnlx,
      mcnly,
      mcnlz,
      mcnlvx,
      mcnlvy,
      mcnlvz,
      mcnlt,
      mcnlsx,
      mcnlsy,
      mcnlsz,
      mcnlp); }
#undef mcabsorbComp
#undef p
#undef sz
#undef sy
#undef sx
#undef t
#undef vz
#undef vy
#undef vx
#undef z
#undef y
#undef x
  mcDEBUG_STATE(
mcnlx,
mcnly,
mcnlz,
mcnlvx,
mcnlvy,
mcnlvz,
mcnlt,
mcnlsx,
mcnlsy,
mcnlsz,
mcnlp)

  /* TRACE Component monitor_2 [16] */
  mccoordschange(mcposrmonitor_2, mcrotrmonitor_2,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component monitor_2 (without coords transformations) */
  mcJumpTrace_monitor_2:
  SIG_MESSAGE("monitor_2 (Trace)");
  mcDEBUG_COMP("monitor_2")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp

#define mcabsorbComp mcabsorbCompmonitor_2
  STORE_NEUTRON(16,
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp);
  mcScattered=0;
  mcRestore=0;
  mcNCounter[16]++;
  mcPCounter[16] += p;
  mcP2Counter[16] += p*p;
#define mccompcurname  monitor_2
#define mccompcurtype  Arm
#define mccompcurindex 16
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompmonitor_2:
  if (RESTORE) /* restore if needed */
  { RESTORE_NEUTRON(16,
      mcnlx,
      mcnly,
      mcnlz,
      mcnlvx,
      mcnlvy,
      mcnlvz,
      mcnlt,
      mcnlsx,
      mcnlsy,
      mcnlsz,
      mcnlp); }
#undef mcabsorbComp
#undef p
#undef sz
#undef sy
#undef sx
#undef t
#undef vz
#undef vy
#undef vx
#undef z
#undef y
#undef x
  mcDEBUG_STATE(
mcnlx,
mcnly,
mcnlz,
mcnlvx,
mcnlvy,
mcnlvz,
mcnlt,
mcnlsx,
mcnlsy,
mcnlsz,
mcnlp)

  /* TRACE Component psd_monitor_f_zoom [17] */
  mccoordschange(mcposrpsd_monitor_f_zoom, mcrotrpsd_monitor_f_zoom,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component psd_monitor_f_zoom (without coords transformations) */
  mcJumpTrace_psd_monitor_f_zoom:
  SIG_MESSAGE("psd_monitor_f_zoom (Trace)");
  mcDEBUG_COMP("psd_monitor_f_zoom")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp

#define mcabsorbComp mcabsorbComppsd_monitor_f_zoom
  STORE_NEUTRON(17,
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp);
  mcScattered=0;
  mcRestore=0;
  mcNCounter[17]++;
  mcPCounter[17] += p;
  mcP2Counter[17] += p*p;
#define mccompcurname  psd_monitor_f_zoom
#define mccompcurtype  PSD_monitor
#define mccompcurindex 17
#define PSD_N mccpsd_monitor_f_zoom_PSD_N
#define PSD_p mccpsd_monitor_f_zoom_PSD_p
#define PSD_p2 mccpsd_monitor_f_zoom_PSD_p2
{   /* Declarations of psd_monitor_f_zoom=PSD_monitor() SETTING parameters. */
int nx = mccpsd_monitor_f_zoom_nx;
int ny = mccpsd_monitor_f_zoom_ny;
char* filename = mccpsd_monitor_f_zoom_filename;
MCNUM xmin = mccpsd_monitor_f_zoom_xmin;
MCNUM xmax = mccpsd_monitor_f_zoom_xmax;
MCNUM ymin = mccpsd_monitor_f_zoom_ymin;
MCNUM ymax = mccpsd_monitor_f_zoom_ymax;
MCNUM xwidth = mccpsd_monitor_f_zoom_xwidth;
MCNUM yheight = mccpsd_monitor_f_zoom_yheight;
MCNUM restore_neutron = mccpsd_monitor_f_zoom_restore_neutron;
#line 94 "/usr/share/mcstas/2.6.1/tools/Python/mcrun/../mccodelib/../../../monitors/PSD_monitor.comp"
{
  PROP_Z0;
  if (x>xmin && x<xmax && y>ymin && y<ymax){
    int i = floor((x - xmin)*nx/(xmax - xmin));
    int j = floor((y - ymin)*ny/(ymax - ymin));
    PSD_N[i][j]++;
    PSD_p[i][j] += p;
    PSD_p2[i][j] += p*p;
    SCATTER;
  }
  if (restore_neutron) {
    RESTORE_NEUTRON(INDEX_CURRENT_COMP, x, y, z, vx, vy, vz, t, sx, sy, sz, p);
  }
}
#line 15428 "./ess_extraction.c"
}   /* End of psd_monitor_f_zoom=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbComppsd_monitor_f_zoom:
  if (RESTORE) /* restore if needed */
  { RESTORE_NEUTRON(17,
      mcnlx,
      mcnly,
      mcnlz,
      mcnlvx,
      mcnlvy,
      mcnlvz,
      mcnlt,
      mcnlsx,
      mcnlsy,
      mcnlsz,
      mcnlp); }
#undef mcabsorbComp
#undef p
#undef sz
#undef sy
#undef sx
#undef t
#undef vz
#undef vy
#undef vx
#undef z
#undef y
#undef x
  mcDEBUG_STATE(
mcnlx,
mcnly,
mcnlz,
mcnlvx,
mcnlvy,
mcnlvz,
mcnlt,
mcnlsx,
mcnlsy,
mcnlsz,
mcnlp)

  /* TRACE Component psd_monitor_f [18] */
  mccoordschange(mcposrpsd_monitor_f, mcrotrpsd_monitor_f,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component psd_monitor_f (without coords transformations) */
  mcJumpTrace_psd_monitor_f:
  SIG_MESSAGE("psd_monitor_f (Trace)");
  mcDEBUG_COMP("psd_monitor_f")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp

#define mcabsorbComp mcabsorbComppsd_monitor_f
  STORE_NEUTRON(18,
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp);
  mcScattered=0;
  mcRestore=0;
  mcNCounter[18]++;
  mcPCounter[18] += p;
  mcP2Counter[18] += p*p;
#define mccompcurname  psd_monitor_f
#define mccompcurtype  PSD_monitor
#define mccompcurindex 18
#define PSD_N mccpsd_monitor_f_PSD_N
#define PSD_p mccpsd_monitor_f_PSD_p
#define PSD_p2 mccpsd_monitor_f_PSD_p2
{   /* Declarations of psd_monitor_f=PSD_monitor() SETTING parameters. */
int nx = mccpsd_monitor_f_nx;
int ny = mccpsd_monitor_f_ny;
char* filename = mccpsd_monitor_f_filename;
MCNUM xmin = mccpsd_monitor_f_xmin;
MCNUM xmax = mccpsd_monitor_f_xmax;
MCNUM ymin = mccpsd_monitor_f_ymin;
MCNUM ymax = mccpsd_monitor_f_ymax;
MCNUM xwidth = mccpsd_monitor_f_xwidth;
MCNUM yheight = mccpsd_monitor_f_yheight;
MCNUM restore_neutron = mccpsd_monitor_f_restore_neutron;
#line 94 "/usr/share/mcstas/2.6.1/tools/Python/mcrun/../mccodelib/../../../monitors/PSD_monitor.comp"
{
  PROP_Z0;
  if (x>xmin && x<xmax && y>ymin && y<ymax){
    int i = floor((x - xmin)*nx/(xmax - xmin));
    int j = floor((y - ymin)*ny/(ymax - ymin));
    PSD_N[i][j]++;
    PSD_p[i][j] += p;
    PSD_p2[i][j] += p*p;
    SCATTER;
  }
  if (restore_neutron) {
    RESTORE_NEUTRON(INDEX_CURRENT_COMP, x, y, z, vx, vy, vz, t, sx, sy, sz, p);
  }
}
#line 15565 "./ess_extraction.c"
}   /* End of psd_monitor_f=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbComppsd_monitor_f:
  if (RESTORE) /* restore if needed */
  { RESTORE_NEUTRON(18,
      mcnlx,
      mcnly,
      mcnlz,
      mcnlvx,
      mcnlvy,
      mcnlvz,
      mcnlt,
      mcnlsx,
      mcnlsy,
      mcnlsz,
      mcnlp); }
#undef mcabsorbComp
#undef p
#undef sz
#undef sy
#undef sx
#undef t
#undef vz
#undef vy
#undef vx
#undef z
#undef y
#undef x
  mcDEBUG_STATE(
mcnlx,
mcnly,
mcnlz,
mcnlvx,
mcnlvy,
mcnlvz,
mcnlt,
mcnlsx,
mcnlsy,
mcnlsz,
mcnlp)

  /* TRACE Component f_divpos [19] */
  mccoordschange(mcposrf_divpos, mcrotrf_divpos,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component f_divpos (without coords transformations) */
  mcJumpTrace_f_divpos:
  SIG_MESSAGE("f_divpos (Trace)");
  mcDEBUG_COMP("f_divpos")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp

#define mcabsorbComp mcabsorbCompf_divpos
  STORE_NEUTRON(19,
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp);
  mcScattered=0;
  mcRestore=0;
  mcNCounter[19]++;
  mcPCounter[19] += p;
  mcP2Counter[19] += p*p;
#define mccompcurname  f_divpos
#define mccompcurtype  DivPos_monitor
#define mccompcurindex 19
#define nh mccf_divpos_nh
#define ndiv mccf_divpos_ndiv
#define Div_N mccf_divpos_Div_N
#define Div_p mccf_divpos_Div_p
#define Div_p2 mccf_divpos_Div_p2
{   /* Declarations of f_divpos=DivPos_monitor() SETTING parameters. */
char* filename = mccf_divpos_filename;
MCNUM xmin = mccf_divpos_xmin;
MCNUM xmax = mccf_divpos_xmax;
MCNUM ymin = mccf_divpos_ymin;
MCNUM ymax = mccf_divpos_ymax;
MCNUM xwidth = mccf_divpos_xwidth;
MCNUM yheight = mccf_divpos_yheight;
MCNUM maxdiv_h = mccf_divpos_maxdiv_h;
MCNUM restore_neutron = mccf_divpos_restore_neutron;
MCNUM nx = mccf_divpos_nx;
MCNUM ny = mccf_divpos_ny;
MCNUM nz = mccf_divpos_nz;
int nowritefile = mccf_divpos_nowritefile;
#line 92 "/usr/share/mcstas/2.6.1/tools/Python/mcrun/../mccodelib/../../../monitors/DivPos_monitor.comp"
{
    int i,j;
    double div;
    double v, vn;

    PROP_Z0;
    if (x>xmin && x<xmax && y>ymin && y<ymax)
    {
      /* Find length of projection onto the [nx ny nz] axis */
      vn = scalar_prod(vx, vy, vz, nx, ny, nz);
      div = RAD2DEG*atan2(vx,vn);

      if (div < maxdiv_h && div > -maxdiv_h)
      {
        i = floor((x - xmin)*nh/(xmax - xmin));
        j = floor((div + maxdiv_h)*ndiv/(2.0*maxdiv_h));
        Div_N[i][j]++;
        Div_p[i][j] += p;
        Div_p2[i][j] += p*p;
        SCATTER;
      }
    }
    if (restore_neutron) {
      RESTORE_NEUTRON(INDEX_CURRENT_COMP, x, y, z, vx, vy, vz, t, sx, sy, sz, p);
    }
}
#line 15719 "./ess_extraction.c"
}   /* End of f_divpos=DivPos_monitor() SETTING parameter declarations. */
#undef Div_p2
#undef Div_p
#undef Div_N
#undef ndiv
#undef nh
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompf_divpos:
  if (RESTORE) /* restore if needed */
  { RESTORE_NEUTRON(19,
      mcnlx,
      mcnly,
      mcnlz,
      mcnlvx,
      mcnlvy,
      mcnlvz,
      mcnlt,
      mcnlsx,
      mcnlsy,
      mcnlsz,
      mcnlp); }
#undef mcabsorbComp
#undef p
#undef sz
#undef sy
#undef sx
#undef t
#undef vz
#undef vy
#undef vx
#undef z
#undef y
#undef x
  mcDEBUG_STATE(
mcnlx,
mcnly,
mcnlz,
mcnlvx,
mcnlvy,
mcnlvz,
mcnlt,
mcnlsx,
mcnlsy,
mcnlsz,
mcnlp)

  /* TRACE Component divhlambda_monitor_f [20] */
  mccoordschange(mcposrdivhlambda_monitor_f, mcrotrdivhlambda_monitor_f,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component divhlambda_monitor_f (without coords transformations) */
  mcJumpTrace_divhlambda_monitor_f:
  SIG_MESSAGE("divhlambda_monitor_f (Trace)");
  mcDEBUG_COMP("divhlambda_monitor_f")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp

#define mcabsorbComp mcabsorbCompdivhlambda_monitor_f
  STORE_NEUTRON(20,
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp);
  mcScattered=0;
  mcRestore=0;
  mcNCounter[20]++;
  mcPCounter[20] += p;
  mcP2Counter[20] += p*p;
#define mccompcurname  divhlambda_monitor_f
#define mccompcurtype  DivLambda_monitor
#define mccompcurindex 20
#define nL mccdivhlambda_monitor_f_nL
#define nh mccdivhlambda_monitor_f_nh
#define Div_N mccdivhlambda_monitor_f_Div_N
#define Div_p mccdivhlambda_monitor_f_Div_p
#define Div_p2 mccdivhlambda_monitor_f_Div_p2
{   /* Declarations of divhlambda_monitor_f=DivLambda_monitor() SETTING parameters. */
char* filename = mccdivhlambda_monitor_f_filename;
MCNUM xmin = mccdivhlambda_monitor_f_xmin;
MCNUM xmax = mccdivhlambda_monitor_f_xmax;
MCNUM ymin = mccdivhlambda_monitor_f_ymin;
MCNUM ymax = mccdivhlambda_monitor_f_ymax;
MCNUM xwidth = mccdivhlambda_monitor_f_xwidth;
MCNUM yheight = mccdivhlambda_monitor_f_yheight;
MCNUM maxdiv_h = mccdivhlambda_monitor_f_maxdiv_h;
MCNUM Lmin = mccdivhlambda_monitor_f_Lmin;
MCNUM Lmax = mccdivhlambda_monitor_f_Lmax;
MCNUM restore_neutron = mccdivhlambda_monitor_f_restore_neutron;
MCNUM nx = mccdivhlambda_monitor_f_nx;
MCNUM ny = mccdivhlambda_monitor_f_ny;
MCNUM nz = mccdivhlambda_monitor_f_nz;
int nowritefile = mccdivhlambda_monitor_f_nowritefile;
#line 89 "/usr/share/mcstas/2.6.1/tools/Python/mcrun/../mccodelib/../../../monitors/DivLambda_monitor.comp"
{
    int i,j;
    double div;
    double lambda;
    double v, vn;

    PROP_Z0;
    lambda = (2*PI/V2K)/sqrt(vx*vx + vy*vy + vz*vz);
    if (x>xmin && x<xmax && y>ymin && y<ymax &&
        lambda > Lmin && lambda < Lmax)
    {
      /* Find length of projection onto the [nx ny nz] axis */
      vn = scalar_prod(vx, vy, vz, nx, ny, nz);
      div = RAD2DEG*atan2(vx,vn);

      if (div < maxdiv_h && div > -maxdiv_h)
      {
        i = floor((lambda - Lmin)*nL/(Lmax - Lmin));
        j = floor((div + maxdiv_h)*nh/(2.0*maxdiv_h));
        Div_N[i][j]++;
        Div_p[i][j] += p;
        Div_p2[i][j] += p*p;
        SCATTER;
      }
    }
    if (restore_neutron) {
      RESTORE_NEUTRON(INDEX_CURRENT_COMP, x, y, z, vx, vy, vz, t, sx, sy, sz, p);
    }
}
#line 15880 "./ess_extraction.c"
}   /* End of divhlambda_monitor_f=DivLambda_monitor() SETTING parameter declarations. */
#undef Div_p2
#undef Div_p
#undef Div_N
#undef nh
#undef nL
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompdivhlambda_monitor_f:
  if (RESTORE) /* restore if needed */
  { RESTORE_NEUTRON(20,
      mcnlx,
      mcnly,
      mcnlz,
      mcnlvx,
      mcnlvy,
      mcnlvz,
      mcnlt,
      mcnlsx,
      mcnlsy,
      mcnlsz,
      mcnlp); }
#undef mcabsorbComp
#undef p
#undef sz
#undef sy
#undef sx
#undef t
#undef vz
#undef vy
#undef vx
#undef z
#undef y
#undef x
  mcDEBUG_STATE(
mcnlx,
mcnly,
mcnlz,
mcnlvx,
mcnlvy,
mcnlvz,
mcnlt,
mcnlsx,
mcnlsy,
mcnlsz,
mcnlp)

  mcabsorbAll:
  mcDEBUG_LEAVE()
  mcDEBUG_STATE(
mcnlx,
mcnly,
mcnlz,
mcnlvx,
mcnlvy,
mcnlvz,
mcnlt,
mcnlsx,
mcnlsy,
mcnlsz,
mcnlp)
  /* Copy neutron state to global variables. */
  mcnx = mcnlx;
  mcny = mcnly;
  mcnz = mcnlz;
  mcnvx = mcnlvx;
  mcnvy = mcnlvy;
  mcnvz = mcnlvz;
  mcnt = mcnlt;
  mcnsx = mcnlsx;
  mcnsy = mcnlsy;
  mcnsz = mcnlsz;
  mcnp = mcnlp;

} /* end trace */

void mcsave(FILE *handle) {
  if (!handle) mcsiminfo_init(NULL);
  /* User component SAVE code. */

  /* User SAVE code for component 'origin'. */
  SIG_MESSAGE("origin (Save)");
#define mccompcurname  origin
#define mccompcurtype  Progress_bar
#define mccompcurindex 1
#define IntermediateCnts mccorigin_IntermediateCnts
#define StartTime mccorigin_StartTime
#define EndTime mccorigin_EndTime
#define CurrentTime mccorigin_CurrentTime
{   /* Declarations of origin=Progress_bar() SETTING parameters. */
char* profile = mccorigin_profile;
MCNUM percent = mccorigin_percent;
MCNUM flag_save = mccorigin_flag_save;
MCNUM minutes = mccorigin_minutes;
#line 115 "/usr/share/mcstas/2.6.1/tools/Python/mcrun/../mccodelib/../../../misc/Progress_bar.comp"
{
  MPI_MASTER(fprintf(stdout, "\nSave [%s]\n", mcinstrument_name););
  if (profile && strlen(profile) && strcmp(profile,"NULL") && strcmp(profile,"0")) {
    char filename[256];
    if (!strlen(profile) || !strcmp(profile,"NULL") || !strcmp(profile,"0")) strcpy(filename, mcinstrument_name);
    else strcpy(filename, profile);
    DETECTOR_OUT_1D(
        "Intensity profiler",
        "Component index [1]",
        "Intensity",
        "prof", 1, mcNUMCOMP, mcNUMCOMP-1,
        &mcNCounter[1],&mcPCounter[1],&mcP2Counter[1],
        filename);

  }
}
#line 15994 "./ess_extraction.c"
}   /* End of origin=Progress_bar() SETTING parameter declarations. */
#undef CurrentTime
#undef EndTime
#undef StartTime
#undef IntermediateCnts
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* User SAVE code for component 'psd_monitor_source'. */
  SIG_MESSAGE("psd_monitor_source (Save)");
#define mccompcurname  psd_monitor_source
#define mccompcurtype  PSD_monitor
#define mccompcurindex 4
#define PSD_N mccpsd_monitor_source_PSD_N
#define PSD_p mccpsd_monitor_source_PSD_p
#define PSD_p2 mccpsd_monitor_source_PSD_p2
{   /* Declarations of psd_monitor_source=PSD_monitor() SETTING parameters. */
int nx = mccpsd_monitor_source_nx;
int ny = mccpsd_monitor_source_ny;
char* filename = mccpsd_monitor_source_filename;
MCNUM xmin = mccpsd_monitor_source_xmin;
MCNUM xmax = mccpsd_monitor_source_xmax;
MCNUM ymin = mccpsd_monitor_source_ymin;
MCNUM ymax = mccpsd_monitor_source_ymax;
MCNUM xwidth = mccpsd_monitor_source_xwidth;
MCNUM yheight = mccpsd_monitor_source_yheight;
MCNUM restore_neutron = mccpsd_monitor_source_restore_neutron;
#line 110 "/usr/share/mcstas/2.6.1/tools/Python/mcrun/../mccodelib/../../../monitors/PSD_monitor.comp"
{
  DETECTOR_OUT_2D(
    "PSD monitor",
    "X position [cm]",
    "Y position [cm]",
    xmin*100.0, xmax*100.0, ymin*100.0, ymax*100.0,
    nx, ny,
    &PSD_N[0][0],&PSD_p[0][0],&PSD_p2[0][0],
    filename);
}
#line 16034 "./ess_extraction.c"
}   /* End of psd_monitor_source=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* User SAVE code for component 'psd_monitor_afteropticsource'. */
  SIG_MESSAGE("psd_monitor_afteropticsource (Save)");
#define mccompcurname  psd_monitor_afteropticsource
#define mccompcurtype  PSD_monitor
#define mccompcurindex 8
#define PSD_N mccpsd_monitor_afteropticsource_PSD_N
#define PSD_p mccpsd_monitor_afteropticsource_PSD_p
#define PSD_p2 mccpsd_monitor_afteropticsource_PSD_p2
{   /* Declarations of psd_monitor_afteropticsource=PSD_monitor() SETTING parameters. */
int nx = mccpsd_monitor_afteropticsource_nx;
int ny = mccpsd_monitor_afteropticsource_ny;
char* filename = mccpsd_monitor_afteropticsource_filename;
MCNUM xmin = mccpsd_monitor_afteropticsource_xmin;
MCNUM xmax = mccpsd_monitor_afteropticsource_xmax;
MCNUM ymin = mccpsd_monitor_afteropticsource_ymin;
MCNUM ymax = mccpsd_monitor_afteropticsource_ymax;
MCNUM xwidth = mccpsd_monitor_afteropticsource_xwidth;
MCNUM yheight = mccpsd_monitor_afteropticsource_yheight;
MCNUM restore_neutron = mccpsd_monitor_afteropticsource_restore_neutron;
#line 110 "/usr/share/mcstas/2.6.1/tools/Python/mcrun/../mccodelib/../../../monitors/PSD_monitor.comp"
{
  DETECTOR_OUT_2D(
    "PSD monitor",
    "X position [cm]",
    "Y position [cm]",
    xmin*100.0, xmax*100.0, ymin*100.0, ymax*100.0,
    nx, ny,
    &PSD_N[0][0],&PSD_p[0][0],&PSD_p2[0][0],
    filename);
}
#line 16073 "./ess_extraction.c"
}   /* End of psd_monitor_afteropticsource=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* User SAVE code for component 'divpos_monitor'. */
  SIG_MESSAGE("divpos_monitor (Save)");
#define mccompcurname  divpos_monitor
#define mccompcurtype  DivPos_monitor
#define mccompcurindex 11
#define nh mccdivpos_monitor_nh
#define ndiv mccdivpos_monitor_ndiv
#define Div_N mccdivpos_monitor_Div_N
#define Div_p mccdivpos_monitor_Div_p
#define Div_p2 mccdivpos_monitor_Div_p2
{   /* Declarations of divpos_monitor=DivPos_monitor() SETTING parameters. */
char* filename = mccdivpos_monitor_filename;
MCNUM xmin = mccdivpos_monitor_xmin;
MCNUM xmax = mccdivpos_monitor_xmax;
MCNUM ymin = mccdivpos_monitor_ymin;
MCNUM ymax = mccdivpos_monitor_ymax;
MCNUM xwidth = mccdivpos_monitor_xwidth;
MCNUM yheight = mccdivpos_monitor_yheight;
MCNUM maxdiv_h = mccdivpos_monitor_maxdiv_h;
MCNUM restore_neutron = mccdivpos_monitor_restore_neutron;
MCNUM nx = mccdivpos_monitor_nx;
MCNUM ny = mccdivpos_monitor_ny;
MCNUM nz = mccdivpos_monitor_nz;
int nowritefile = mccdivpos_monitor_nowritefile;
#line 120 "/usr/share/mcstas/2.6.1/tools/Python/mcrun/../mccodelib/../../../monitors/DivPos_monitor.comp"
{
    if (!nowritefile) {
    DETECTOR_OUT_2D(
        "Position-divergence monitor",
        "pos [m]",
        "divergence [deg]",
        xmin, xmax, -maxdiv_h, maxdiv_h,
        nh, ndiv,
        &Div_N[0][0],&Div_p[0][0],&Div_p2[0][0],
        filename);
    }
}
#line 16119 "./ess_extraction.c"
}   /* End of divpos_monitor=DivPos_monitor() SETTING parameter declarations. */
#undef Div_p2
#undef Div_p
#undef Div_N
#undef ndiv
#undef nh
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* User SAVE code for component 'psd_monitor_g1'. */
  SIG_MESSAGE("psd_monitor_g1 (Save)");
#define mccompcurname  psd_monitor_g1
#define mccompcurtype  PSD_monitor
#define mccompcurindex 12
#define PSD_N mccpsd_monitor_g1_PSD_N
#define PSD_p mccpsd_monitor_g1_PSD_p
#define PSD_p2 mccpsd_monitor_g1_PSD_p2
{   /* Declarations of psd_monitor_g1=PSD_monitor() SETTING parameters. */
int nx = mccpsd_monitor_g1_nx;
int ny = mccpsd_monitor_g1_ny;
char* filename = mccpsd_monitor_g1_filename;
MCNUM xmin = mccpsd_monitor_g1_xmin;
MCNUM xmax = mccpsd_monitor_g1_xmax;
MCNUM ymin = mccpsd_monitor_g1_ymin;
MCNUM ymax = mccpsd_monitor_g1_ymax;
MCNUM xwidth = mccpsd_monitor_g1_xwidth;
MCNUM yheight = mccpsd_monitor_g1_yheight;
MCNUM restore_neutron = mccpsd_monitor_g1_restore_neutron;
#line 110 "/usr/share/mcstas/2.6.1/tools/Python/mcrun/../mccodelib/../../../monitors/PSD_monitor.comp"
{
  DETECTOR_OUT_2D(
    "PSD monitor",
    "X position [cm]",
    "Y position [cm]",
    xmin*100.0, xmax*100.0, ymin*100.0, ymax*100.0,
    nx, ny,
    &PSD_N[0][0],&PSD_p[0][0],&PSD_p2[0][0],
    filename);
}
#line 16160 "./ess_extraction.c"
}   /* End of psd_monitor_g1=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* User SAVE code for component 'divhlambda_monitor_g1'. */
  SIG_MESSAGE("divhlambda_monitor_g1 (Save)");
#define mccompcurname  divhlambda_monitor_g1
#define mccompcurtype  DivLambda_monitor
#define mccompcurindex 13
#define nL mccdivhlambda_monitor_g1_nL
#define nh mccdivhlambda_monitor_g1_nh
#define Div_N mccdivhlambda_monitor_g1_Div_N
#define Div_p mccdivhlambda_monitor_g1_Div_p
#define Div_p2 mccdivhlambda_monitor_g1_Div_p2
{   /* Declarations of divhlambda_monitor_g1=DivLambda_monitor() SETTING parameters. */
char* filename = mccdivhlambda_monitor_g1_filename;
MCNUM xmin = mccdivhlambda_monitor_g1_xmin;
MCNUM xmax = mccdivhlambda_monitor_g1_xmax;
MCNUM ymin = mccdivhlambda_monitor_g1_ymin;
MCNUM ymax = mccdivhlambda_monitor_g1_ymax;
MCNUM xwidth = mccdivhlambda_monitor_g1_xwidth;
MCNUM yheight = mccdivhlambda_monitor_g1_yheight;
MCNUM maxdiv_h = mccdivhlambda_monitor_g1_maxdiv_h;
MCNUM Lmin = mccdivhlambda_monitor_g1_Lmin;
MCNUM Lmax = mccdivhlambda_monitor_g1_Lmax;
MCNUM restore_neutron = mccdivhlambda_monitor_g1_restore_neutron;
MCNUM nx = mccdivhlambda_monitor_g1_nx;
MCNUM ny = mccdivhlambda_monitor_g1_ny;
MCNUM nz = mccdivhlambda_monitor_g1_nz;
int nowritefile = mccdivhlambda_monitor_g1_nowritefile;
#line 119 "/usr/share/mcstas/2.6.1/tools/Python/mcrun/../mccodelib/../../../monitors/DivLambda_monitor.comp"
{
    if (!nowritefile) {
    DETECTOR_OUT_2D(
        "Wavelength-divergence monitor",
        "Wavelength [AA]",
        "divergence [deg]",
        Lmin, Lmax, -maxdiv_h, maxdiv_h,
        nL, nh,
        &Div_N[0][0],&Div_p[0][0],&Div_p2[0][0],
        filename);
    }
}
#line 16208 "./ess_extraction.c"
}   /* End of divhlambda_monitor_g1=DivLambda_monitor() SETTING parameter declarations. */
#undef Div_p2
#undef Div_p
#undef Div_N
#undef nh
#undef nL
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* User SAVE code for component 'psd_monitor_f_zoom'. */
  SIG_MESSAGE("psd_monitor_f_zoom (Save)");
#define mccompcurname  psd_monitor_f_zoom
#define mccompcurtype  PSD_monitor
#define mccompcurindex 17
#define PSD_N mccpsd_monitor_f_zoom_PSD_N
#define PSD_p mccpsd_monitor_f_zoom_PSD_p
#define PSD_p2 mccpsd_monitor_f_zoom_PSD_p2
{   /* Declarations of psd_monitor_f_zoom=PSD_monitor() SETTING parameters. */
int nx = mccpsd_monitor_f_zoom_nx;
int ny = mccpsd_monitor_f_zoom_ny;
char* filename = mccpsd_monitor_f_zoom_filename;
MCNUM xmin = mccpsd_monitor_f_zoom_xmin;
MCNUM xmax = mccpsd_monitor_f_zoom_xmax;
MCNUM ymin = mccpsd_monitor_f_zoom_ymin;
MCNUM ymax = mccpsd_monitor_f_zoom_ymax;
MCNUM xwidth = mccpsd_monitor_f_zoom_xwidth;
MCNUM yheight = mccpsd_monitor_f_zoom_yheight;
MCNUM restore_neutron = mccpsd_monitor_f_zoom_restore_neutron;
#line 110 "/usr/share/mcstas/2.6.1/tools/Python/mcrun/../mccodelib/../../../monitors/PSD_monitor.comp"
{
  DETECTOR_OUT_2D(
    "PSD monitor",
    "X position [cm]",
    "Y position [cm]",
    xmin*100.0, xmax*100.0, ymin*100.0, ymax*100.0,
    nx, ny,
    &PSD_N[0][0],&PSD_p[0][0],&PSD_p2[0][0],
    filename);
}
#line 16249 "./ess_extraction.c"
}   /* End of psd_monitor_f_zoom=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* User SAVE code for component 'psd_monitor_f'. */
  SIG_MESSAGE("psd_monitor_f (Save)");
#define mccompcurname  psd_monitor_f
#define mccompcurtype  PSD_monitor
#define mccompcurindex 18
#define PSD_N mccpsd_monitor_f_PSD_N
#define PSD_p mccpsd_monitor_f_PSD_p
#define PSD_p2 mccpsd_monitor_f_PSD_p2
{   /* Declarations of psd_monitor_f=PSD_monitor() SETTING parameters. */
int nx = mccpsd_monitor_f_nx;
int ny = mccpsd_monitor_f_ny;
char* filename = mccpsd_monitor_f_filename;
MCNUM xmin = mccpsd_monitor_f_xmin;
MCNUM xmax = mccpsd_monitor_f_xmax;
MCNUM ymin = mccpsd_monitor_f_ymin;
MCNUM ymax = mccpsd_monitor_f_ymax;
MCNUM xwidth = mccpsd_monitor_f_xwidth;
MCNUM yheight = mccpsd_monitor_f_yheight;
MCNUM restore_neutron = mccpsd_monitor_f_restore_neutron;
#line 110 "/usr/share/mcstas/2.6.1/tools/Python/mcrun/../mccodelib/../../../monitors/PSD_monitor.comp"
{
  DETECTOR_OUT_2D(
    "PSD monitor",
    "X position [cm]",
    "Y position [cm]",
    xmin*100.0, xmax*100.0, ymin*100.0, ymax*100.0,
    nx, ny,
    &PSD_N[0][0],&PSD_p[0][0],&PSD_p2[0][0],
    filename);
}
#line 16288 "./ess_extraction.c"
}   /* End of psd_monitor_f=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* User SAVE code for component 'f_divpos'. */
  SIG_MESSAGE("f_divpos (Save)");
#define mccompcurname  f_divpos
#define mccompcurtype  DivPos_monitor
#define mccompcurindex 19
#define nh mccf_divpos_nh
#define ndiv mccf_divpos_ndiv
#define Div_N mccf_divpos_Div_N
#define Div_p mccf_divpos_Div_p
#define Div_p2 mccf_divpos_Div_p2
{   /* Declarations of f_divpos=DivPos_monitor() SETTING parameters. */
char* filename = mccf_divpos_filename;
MCNUM xmin = mccf_divpos_xmin;
MCNUM xmax = mccf_divpos_xmax;
MCNUM ymin = mccf_divpos_ymin;
MCNUM ymax = mccf_divpos_ymax;
MCNUM xwidth = mccf_divpos_xwidth;
MCNUM yheight = mccf_divpos_yheight;
MCNUM maxdiv_h = mccf_divpos_maxdiv_h;
MCNUM restore_neutron = mccf_divpos_restore_neutron;
MCNUM nx = mccf_divpos_nx;
MCNUM ny = mccf_divpos_ny;
MCNUM nz = mccf_divpos_nz;
int nowritefile = mccf_divpos_nowritefile;
#line 120 "/usr/share/mcstas/2.6.1/tools/Python/mcrun/../mccodelib/../../../monitors/DivPos_monitor.comp"
{
    if (!nowritefile) {
    DETECTOR_OUT_2D(
        "Position-divergence monitor",
        "pos [m]",
        "divergence [deg]",
        xmin, xmax, -maxdiv_h, maxdiv_h,
        nh, ndiv,
        &Div_N[0][0],&Div_p[0][0],&Div_p2[0][0],
        filename);
    }
}
#line 16334 "./ess_extraction.c"
}   /* End of f_divpos=DivPos_monitor() SETTING parameter declarations. */
#undef Div_p2
#undef Div_p
#undef Div_N
#undef ndiv
#undef nh
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* User SAVE code for component 'divhlambda_monitor_f'. */
  SIG_MESSAGE("divhlambda_monitor_f (Save)");
#define mccompcurname  divhlambda_monitor_f
#define mccompcurtype  DivLambda_monitor
#define mccompcurindex 20
#define nL mccdivhlambda_monitor_f_nL
#define nh mccdivhlambda_monitor_f_nh
#define Div_N mccdivhlambda_monitor_f_Div_N
#define Div_p mccdivhlambda_monitor_f_Div_p
#define Div_p2 mccdivhlambda_monitor_f_Div_p2
{   /* Declarations of divhlambda_monitor_f=DivLambda_monitor() SETTING parameters. */
char* filename = mccdivhlambda_monitor_f_filename;
MCNUM xmin = mccdivhlambda_monitor_f_xmin;
MCNUM xmax = mccdivhlambda_monitor_f_xmax;
MCNUM ymin = mccdivhlambda_monitor_f_ymin;
MCNUM ymax = mccdivhlambda_monitor_f_ymax;
MCNUM xwidth = mccdivhlambda_monitor_f_xwidth;
MCNUM yheight = mccdivhlambda_monitor_f_yheight;
MCNUM maxdiv_h = mccdivhlambda_monitor_f_maxdiv_h;
MCNUM Lmin = mccdivhlambda_monitor_f_Lmin;
MCNUM Lmax = mccdivhlambda_monitor_f_Lmax;
MCNUM restore_neutron = mccdivhlambda_monitor_f_restore_neutron;
MCNUM nx = mccdivhlambda_monitor_f_nx;
MCNUM ny = mccdivhlambda_monitor_f_ny;
MCNUM nz = mccdivhlambda_monitor_f_nz;
int nowritefile = mccdivhlambda_monitor_f_nowritefile;
#line 119 "/usr/share/mcstas/2.6.1/tools/Python/mcrun/../mccodelib/../../../monitors/DivLambda_monitor.comp"
{
    if (!nowritefile) {
    DETECTOR_OUT_2D(
        "Wavelength-divergence monitor",
        "Wavelength [AA]",
        "divergence [deg]",
        Lmin, Lmax, -maxdiv_h, maxdiv_h,
        nL, nh,
        &Div_N[0][0],&Div_p[0][0],&Div_p2[0][0],
        filename);
    }
}
#line 16384 "./ess_extraction.c"
}   /* End of divhlambda_monitor_f=DivLambda_monitor() SETTING parameter declarations. */
#undef Div_p2
#undef Div_p
#undef Div_N
#undef nh
#undef nL
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  if (!handle) mcsiminfo_close(); 
} /* end save */
void mcfinally(void) {
  /* User component FINALLY code. */
  mcsiminfo_init(NULL);
  mcsave(mcsiminfo_file); /* save data when simulation ends */

  /* User FINALLY code for component 'origin'. */
  SIG_MESSAGE("origin (Finally)");
#define mccompcurname  origin
#define mccompcurtype  Progress_bar
#define mccompcurindex 1
#define IntermediateCnts mccorigin_IntermediateCnts
#define StartTime mccorigin_StartTime
#define EndTime mccorigin_EndTime
#define CurrentTime mccorigin_CurrentTime
{   /* Declarations of origin=Progress_bar() SETTING parameters. */
char* profile = mccorigin_profile;
MCNUM percent = mccorigin_percent;
MCNUM flag_save = mccorigin_flag_save;
MCNUM minutes = mccorigin_minutes;
#line 133 "/usr/share/mcstas/2.6.1/tools/Python/mcrun/../mccodelib/../../../misc/Progress_bar.comp"
{
  time_t NowTime;
  time(&NowTime);
  fprintf(stdout, "\nFinally [%s: %s]. Time: ", mcinstrument_name, mcdirname ? mcdirname : ".");
  if (difftime(NowTime,StartTime) < 60.0)
    fprintf(stdout, "%g [s] ", difftime(NowTime,StartTime));
  else if (difftime(NowTime,StartTime) > 3600.0)
    fprintf(stdout, "%g [h] ", difftime(NowTime,StartTime)/3660.0);
  else
    fprintf(stdout, "%g [min] ", difftime(NowTime,StartTime)/60.0);
  fprintf(stdout, "\n");
}
#line 16429 "./ess_extraction.c"
}   /* End of origin=Progress_bar() SETTING parameter declarations. */
#undef CurrentTime
#undef EndTime
#undef StartTime
#undef IntermediateCnts
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

    if (!mcNCounter[1]) fprintf(stderr, "Warning: No neutron could reach Component[1] origin\n");
    if (mcAbsorbProp[1]) fprintf(stderr, "Warning: %g events were removed in Component[1] origin=Progress_bar()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[1]);
    if (!mcNCounter[2]) fprintf(stderr, "Warning: No neutron could reach Component[2] source\n");
    if (mcAbsorbProp[2]) fprintf(stderr, "Warning: %g events were removed in Component[2] source=Arm()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[2]);
    if (!mcNCounter[3]) fprintf(stderr, "Warning: No neutron could reach Component[3] source_div\n");
    if (mcAbsorbProp[3]) fprintf(stderr, "Warning: %g events were removed in Component[3] source_div=Source_div()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[3]);
  /* User FINALLY code for component 'psd_monitor_source'. */
  SIG_MESSAGE("psd_monitor_source (Finally)");
#define mccompcurname  psd_monitor_source
#define mccompcurtype  PSD_monitor
#define mccompcurindex 4
#define PSD_N mccpsd_monitor_source_PSD_N
#define PSD_p mccpsd_monitor_source_PSD_p
#define PSD_p2 mccpsd_monitor_source_PSD_p2
{   /* Declarations of psd_monitor_source=PSD_monitor() SETTING parameters. */
int nx = mccpsd_monitor_source_nx;
int ny = mccpsd_monitor_source_ny;
char* filename = mccpsd_monitor_source_filename;
MCNUM xmin = mccpsd_monitor_source_xmin;
MCNUM xmax = mccpsd_monitor_source_xmax;
MCNUM ymin = mccpsd_monitor_source_ymin;
MCNUM ymax = mccpsd_monitor_source_ymax;
MCNUM xwidth = mccpsd_monitor_source_xwidth;
MCNUM yheight = mccpsd_monitor_source_yheight;
MCNUM restore_neutron = mccpsd_monitor_source_restore_neutron;
#line 122 "/usr/share/mcstas/2.6.1/tools/Python/mcrun/../mccodelib/../../../monitors/PSD_monitor.comp"
{
  destroy_darr2d(PSD_N);
  destroy_darr2d(PSD_p);
  destroy_darr2d(PSD_p2);
}
#line 16470 "./ess_extraction.c"
}   /* End of psd_monitor_source=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

    if (!mcNCounter[4]) fprintf(stderr, "Warning: No neutron could reach Component[4] psd_monitor_source\n");
    if (mcAbsorbProp[4]) fprintf(stderr, "Warning: %g events were removed in Component[4] psd_monitor_source=PSD_monitor()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[4]);
  /* User FINALLY code for component 'parabolic_optic_before_guide_v'. */
  SIG_MESSAGE("parabolic_optic_before_guide_v (Finally)");
#define mccompcurname  parabolic_optic_before_guide_v
#define mccompcurtype  FlatEllipse_finite_mirror
#define mccompcurindex 5
#define reflect mccparabolic_optic_before_guide_v_reflect
#define s mccparabolic_optic_before_guide_v_s
#define pTable mccparabolic_optic_before_guide_v_pTable
#define R0 mccparabolic_optic_before_guide_v_R0
#define Qc mccparabolic_optic_before_guide_v_Qc
#define W mccparabolic_optic_before_guide_v_W
#define alpha mccparabolic_optic_before_guide_v_alpha
#define transmit mccparabolic_optic_before_guide_v_transmit
{   /* Declarations of parabolic_optic_before_guide_v=FlatEllipse_finite_mirror() SETTING parameters. */
MCNUM sourceDist = mccparabolic_optic_before_guide_v_sourceDist;
MCNUM LStart = mccparabolic_optic_before_guide_v_LStart;
MCNUM LEnd = mccparabolic_optic_before_guide_v_LEnd;
MCNUM lStart = mccparabolic_optic_before_guide_v_lStart;
MCNUM lEnd = mccparabolic_optic_before_guide_v_lEnd;
MCNUM r_0 = mccparabolic_optic_before_guide_v_r_0;
MCNUM nummirror = mccparabolic_optic_before_guide_v_nummirror;
MCNUM mf = mccparabolic_optic_before_guide_v_mf;
MCNUM mb = mccparabolic_optic_before_guide_v_mb;
MCNUM mirror_width = mccparabolic_optic_before_guide_v_mirror_width;
MCNUM doubleReflections = mccparabolic_optic_before_guide_v_doubleReflections;
#line 211 "FlatEllipse_finite_mirror.comp"
{
    //Mainly Writes Inline Detector Data
    finishSimulation(&s);
}
#line 16511 "./ess_extraction.c"
}   /* End of parabolic_optic_before_guide_v=FlatEllipse_finite_mirror() SETTING parameter declarations. */
#undef transmit
#undef alpha
#undef W
#undef Qc
#undef R0
#undef pTable
#undef s
#undef reflect
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

    if (!mcNCounter[5]) fprintf(stderr, "Warning: No neutron could reach Component[5] parabolic_optic_before_guide_v\n");
    if (mcAbsorbProp[5]) fprintf(stderr, "Warning: %g events were removed in Component[5] parabolic_optic_before_guide_v=FlatEllipse_finite_mirror()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[5]);
  /* User FINALLY code for component 'parabolic_optic_before_guide_h'. */
  SIG_MESSAGE("parabolic_optic_before_guide_h (Finally)");
#define mccompcurname  parabolic_optic_before_guide_h
#define mccompcurtype  FlatEllipse_finite_mirror
#define mccompcurindex 6
#define reflect mccparabolic_optic_before_guide_h_reflect
#define s mccparabolic_optic_before_guide_h_s
#define pTable mccparabolic_optic_before_guide_h_pTable
#define R0 mccparabolic_optic_before_guide_h_R0
#define Qc mccparabolic_optic_before_guide_h_Qc
#define W mccparabolic_optic_before_guide_h_W
#define alpha mccparabolic_optic_before_guide_h_alpha
#define transmit mccparabolic_optic_before_guide_h_transmit
{   /* Declarations of parabolic_optic_before_guide_h=FlatEllipse_finite_mirror() SETTING parameters. */
MCNUM sourceDist = mccparabolic_optic_before_guide_h_sourceDist;
MCNUM LStart = mccparabolic_optic_before_guide_h_LStart;
MCNUM LEnd = mccparabolic_optic_before_guide_h_LEnd;
MCNUM lStart = mccparabolic_optic_before_guide_h_lStart;
MCNUM lEnd = mccparabolic_optic_before_guide_h_lEnd;
MCNUM r_0 = mccparabolic_optic_before_guide_h_r_0;
MCNUM nummirror = mccparabolic_optic_before_guide_h_nummirror;
MCNUM mf = mccparabolic_optic_before_guide_h_mf;
MCNUM mb = mccparabolic_optic_before_guide_h_mb;
MCNUM mirror_width = mccparabolic_optic_before_guide_h_mirror_width;
MCNUM doubleReflections = mccparabolic_optic_before_guide_h_doubleReflections;
#line 211 "FlatEllipse_finite_mirror.comp"
{
    //Mainly Writes Inline Detector Data
    finishSimulation(&s);
}
#line 16557 "./ess_extraction.c"
}   /* End of parabolic_optic_before_guide_h=FlatEllipse_finite_mirror() SETTING parameter declarations. */
#undef transmit
#undef alpha
#undef W
#undef Qc
#undef R0
#undef pTable
#undef s
#undef reflect
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

    if (!mcNCounter[6]) fprintf(stderr, "Warning: No neutron could reach Component[6] parabolic_optic_before_guide_h\n");
    if (mcAbsorbProp[6]) fprintf(stderr, "Warning: %g events were removed in Component[6] parabolic_optic_before_guide_h=FlatEllipse_finite_mirror()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[6]);
    if (!mcNCounter[7]) fprintf(stderr, "Warning: No neutron could reach Component[7] after_optic_source\n");
    if (mcAbsorbProp[7]) fprintf(stderr, "Warning: %g events were removed in Component[7] after_optic_source=Arm()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[7]);
  /* User FINALLY code for component 'psd_monitor_afteropticsource'. */
  SIG_MESSAGE("psd_monitor_afteropticsource (Finally)");
#define mccompcurname  psd_monitor_afteropticsource
#define mccompcurtype  PSD_monitor
#define mccompcurindex 8
#define PSD_N mccpsd_monitor_afteropticsource_PSD_N
#define PSD_p mccpsd_monitor_afteropticsource_PSD_p
#define PSD_p2 mccpsd_monitor_afteropticsource_PSD_p2
{   /* Declarations of psd_monitor_afteropticsource=PSD_monitor() SETTING parameters. */
int nx = mccpsd_monitor_afteropticsource_nx;
int ny = mccpsd_monitor_afteropticsource_ny;
char* filename = mccpsd_monitor_afteropticsource_filename;
MCNUM xmin = mccpsd_monitor_afteropticsource_xmin;
MCNUM xmax = mccpsd_monitor_afteropticsource_xmax;
MCNUM ymin = mccpsd_monitor_afteropticsource_ymin;
MCNUM ymax = mccpsd_monitor_afteropticsource_ymax;
MCNUM xwidth = mccpsd_monitor_afteropticsource_xwidth;
MCNUM yheight = mccpsd_monitor_afteropticsource_yheight;
MCNUM restore_neutron = mccpsd_monitor_afteropticsource_restore_neutron;
#line 122 "/usr/share/mcstas/2.6.1/tools/Python/mcrun/../mccodelib/../../../monitors/PSD_monitor.comp"
{
  destroy_darr2d(PSD_N);
  destroy_darr2d(PSD_p);
  destroy_darr2d(PSD_p2);
}
#line 16600 "./ess_extraction.c"
}   /* End of psd_monitor_afteropticsource=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

    if (!mcNCounter[8]) fprintf(stderr, "Warning: No neutron could reach Component[8] psd_monitor_afteropticsource\n");
    if (mcAbsorbProp[8]) fprintf(stderr, "Warning: %g events were removed in Component[8] psd_monitor_afteropticsource=PSD_monitor()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[8]);
  /* User FINALLY code for component 'guide_gravity_1'. */
  SIG_MESSAGE("guide_gravity_1 (Finally)");
#define mccompcurname  guide_gravity_1
#define mccompcurtype  Guide_gravity
#define mccompcurindex 9
#define GVars mccguide_gravity_1_GVars
#define pTable mccguide_gravity_1_pTable
{   /* Declarations of guide_gravity_1=Guide_gravity() SETTING parameters. */
MCNUM w1 = mccguide_gravity_1_w1;
MCNUM h1 = mccguide_gravity_1_h1;
MCNUM w2 = mccguide_gravity_1_w2;
MCNUM h2 = mccguide_gravity_1_h2;
MCNUM l = mccguide_gravity_1_l;
MCNUM R0 = mccguide_gravity_1_R0;
MCNUM Qc = mccguide_gravity_1_Qc;
MCNUM alpha = mccguide_gravity_1_alpha;
MCNUM m = mccguide_gravity_1_m;
MCNUM W = mccguide_gravity_1_W;
MCNUM nslit = mccguide_gravity_1_nslit;
MCNUM d = mccguide_gravity_1_d;
MCNUM mleft = mccguide_gravity_1_mleft;
MCNUM mright = mccguide_gravity_1_mright;
MCNUM mtop = mccguide_gravity_1_mtop;
MCNUM mbottom = mccguide_gravity_1_mbottom;
MCNUM nhslit = mccguide_gravity_1_nhslit;
MCNUM G = mccguide_gravity_1_G;
MCNUM aleft = mccguide_gravity_1_aleft;
MCNUM aright = mccguide_gravity_1_aright;
MCNUM atop = mccguide_gravity_1_atop;
MCNUM abottom = mccguide_gravity_1_abottom;
MCNUM wavy = mccguide_gravity_1_wavy;
MCNUM wavy_z = mccguide_gravity_1_wavy_z;
MCNUM wavy_tb = mccguide_gravity_1_wavy_tb;
MCNUM wavy_lr = mccguide_gravity_1_wavy_lr;
MCNUM chamfers = mccguide_gravity_1_chamfers;
MCNUM chamfers_z = mccguide_gravity_1_chamfers_z;
MCNUM chamfers_lr = mccguide_gravity_1_chamfers_lr;
MCNUM chamfers_tb = mccguide_gravity_1_chamfers_tb;
MCNUM nelements = mccguide_gravity_1_nelements;
MCNUM nu = mccguide_gravity_1_nu;
MCNUM phase = mccguide_gravity_1_phase;
char* reflect = mccguide_gravity_1_reflect;
#line 562 "/usr/share/mcstas/2.6.1/tools/Python/mcrun/../mccodelib/../../../optics/Guide_gravity.comp"
{
if (GVars.warnings > 100) {
  fprintf(stderr,"%s: warning: neutron has entered guide, but can not exit !\n", GVars.compcurname);
  fprintf(stderr,"%s: warning: This message has been repeated %g times\n", GVars.compcurname, GVars.warnings);
}
}
#line 16660 "./ess_extraction.c"
}   /* End of guide_gravity_1=Guide_gravity() SETTING parameter declarations. */
#undef pTable
#undef GVars
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

    if (!mcNCounter[9]) fprintf(stderr, "Warning: No neutron could reach Component[9] guide_gravity_1\n");
    if (mcAbsorbProp[9]) fprintf(stderr, "Warning: %g events were removed in Component[9] guide_gravity_1=Guide_gravity()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[9]);
    if (!mcNCounter[10]) fprintf(stderr, "Warning: No neutron could reach Component[10] monitor_1\n");
    if (mcAbsorbProp[10]) fprintf(stderr, "Warning: %g events were removed in Component[10] monitor_1=Arm()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[10]);
    if (!mcNCounter[11]) fprintf(stderr, "Warning: No neutron could reach Component[11] divpos_monitor\n");
    if (mcAbsorbProp[11]) fprintf(stderr, "Warning: %g events were removed in Component[11] divpos_monitor=DivPos_monitor()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[11]);
  /* User FINALLY code for component 'psd_monitor_g1'. */
  SIG_MESSAGE("psd_monitor_g1 (Finally)");
#define mccompcurname  psd_monitor_g1
#define mccompcurtype  PSD_monitor
#define mccompcurindex 12
#define PSD_N mccpsd_monitor_g1_PSD_N
#define PSD_p mccpsd_monitor_g1_PSD_p
#define PSD_p2 mccpsd_monitor_g1_PSD_p2
{   /* Declarations of psd_monitor_g1=PSD_monitor() SETTING parameters. */
int nx = mccpsd_monitor_g1_nx;
int ny = mccpsd_monitor_g1_ny;
char* filename = mccpsd_monitor_g1_filename;
MCNUM xmin = mccpsd_monitor_g1_xmin;
MCNUM xmax = mccpsd_monitor_g1_xmax;
MCNUM ymin = mccpsd_monitor_g1_ymin;
MCNUM ymax = mccpsd_monitor_g1_ymax;
MCNUM xwidth = mccpsd_monitor_g1_xwidth;
MCNUM yheight = mccpsd_monitor_g1_yheight;
MCNUM restore_neutron = mccpsd_monitor_g1_restore_neutron;
#line 122 "/usr/share/mcstas/2.6.1/tools/Python/mcrun/../mccodelib/../../../monitors/PSD_monitor.comp"
{
  destroy_darr2d(PSD_N);
  destroy_darr2d(PSD_p);
  destroy_darr2d(PSD_p2);
}
#line 16699 "./ess_extraction.c"
}   /* End of psd_monitor_g1=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

    if (!mcNCounter[12]) fprintf(stderr, "Warning: No neutron could reach Component[12] psd_monitor_g1\n");
    if (mcAbsorbProp[12]) fprintf(stderr, "Warning: %g events were removed in Component[12] psd_monitor_g1=PSD_monitor()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[12]);
    if (!mcNCounter[13]) fprintf(stderr, "Warning: No neutron could reach Component[13] divhlambda_monitor_g1\n");
    if (mcAbsorbProp[13]) fprintf(stderr, "Warning: %g events were removed in Component[13] divhlambda_monitor_g1=DivLambda_monitor()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[13]);
  /* User FINALLY code for component 'parabolic_optic1'. */
  SIG_MESSAGE("parabolic_optic1 (Finally)");
#define mccompcurname  parabolic_optic1
#define mccompcurtype  FlatEllipse_finite_mirror
#define mccompcurindex 14
#define reflect mccparabolic_optic1_reflect
#define s mccparabolic_optic1_s
#define pTable mccparabolic_optic1_pTable
#define R0 mccparabolic_optic1_R0
#define Qc mccparabolic_optic1_Qc
#define W mccparabolic_optic1_W
#define alpha mccparabolic_optic1_alpha
#define transmit mccparabolic_optic1_transmit
{   /* Declarations of parabolic_optic1=FlatEllipse_finite_mirror() SETTING parameters. */
MCNUM sourceDist = mccparabolic_optic1_sourceDist;
MCNUM LStart = mccparabolic_optic1_LStart;
MCNUM LEnd = mccparabolic_optic1_LEnd;
MCNUM lStart = mccparabolic_optic1_lStart;
MCNUM lEnd = mccparabolic_optic1_lEnd;
MCNUM r_0 = mccparabolic_optic1_r_0;
MCNUM nummirror = mccparabolic_optic1_nummirror;
MCNUM mf = mccparabolic_optic1_mf;
MCNUM mb = mccparabolic_optic1_mb;
MCNUM mirror_width = mccparabolic_optic1_mirror_width;
MCNUM doubleReflections = mccparabolic_optic1_doubleReflections;
#line 211 "FlatEllipse_finite_mirror.comp"
{
    //Mainly Writes Inline Detector Data
    finishSimulation(&s);
}
#line 16742 "./ess_extraction.c"
}   /* End of parabolic_optic1=FlatEllipse_finite_mirror() SETTING parameter declarations. */
#undef transmit
#undef alpha
#undef W
#undef Qc
#undef R0
#undef pTable
#undef s
#undef reflect
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

    if (!mcNCounter[14]) fprintf(stderr, "Warning: No neutron could reach Component[14] parabolic_optic1\n");
    if (mcAbsorbProp[14]) fprintf(stderr, "Warning: %g events were removed in Component[14] parabolic_optic1=FlatEllipse_finite_mirror()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[14]);
  /* User FINALLY code for component 'parabolic_optic2'. */
  SIG_MESSAGE("parabolic_optic2 (Finally)");
#define mccompcurname  parabolic_optic2
#define mccompcurtype  FlatEllipse_finite_mirror
#define mccompcurindex 15
#define reflect mccparabolic_optic2_reflect
#define s mccparabolic_optic2_s
#define pTable mccparabolic_optic2_pTable
#define R0 mccparabolic_optic2_R0
#define Qc mccparabolic_optic2_Qc
#define W mccparabolic_optic2_W
#define alpha mccparabolic_optic2_alpha
#define transmit mccparabolic_optic2_transmit
{   /* Declarations of parabolic_optic2=FlatEllipse_finite_mirror() SETTING parameters. */
MCNUM sourceDist = mccparabolic_optic2_sourceDist;
MCNUM LStart = mccparabolic_optic2_LStart;
MCNUM LEnd = mccparabolic_optic2_LEnd;
MCNUM lStart = mccparabolic_optic2_lStart;
MCNUM lEnd = mccparabolic_optic2_lEnd;
MCNUM r_0 = mccparabolic_optic2_r_0;
MCNUM nummirror = mccparabolic_optic2_nummirror;
MCNUM mf = mccparabolic_optic2_mf;
MCNUM mb = mccparabolic_optic2_mb;
MCNUM mirror_width = mccparabolic_optic2_mirror_width;
MCNUM doubleReflections = mccparabolic_optic2_doubleReflections;
#line 211 "FlatEllipse_finite_mirror.comp"
{
    //Mainly Writes Inline Detector Data
    finishSimulation(&s);
}
#line 16788 "./ess_extraction.c"
}   /* End of parabolic_optic2=FlatEllipse_finite_mirror() SETTING parameter declarations. */
#undef transmit
#undef alpha
#undef W
#undef Qc
#undef R0
#undef pTable
#undef s
#undef reflect
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

    if (!mcNCounter[15]) fprintf(stderr, "Warning: No neutron could reach Component[15] parabolic_optic2\n");
    if (mcAbsorbProp[15]) fprintf(stderr, "Warning: %g events were removed in Component[15] parabolic_optic2=FlatEllipse_finite_mirror()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[15]);
    if (!mcNCounter[16]) fprintf(stderr, "Warning: No neutron could reach Component[16] monitor_2\n");
    if (mcAbsorbProp[16]) fprintf(stderr, "Warning: %g events were removed in Component[16] monitor_2=Arm()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[16]);
  /* User FINALLY code for component 'psd_monitor_f_zoom'. */
  SIG_MESSAGE("psd_monitor_f_zoom (Finally)");
#define mccompcurname  psd_monitor_f_zoom
#define mccompcurtype  PSD_monitor
#define mccompcurindex 17
#define PSD_N mccpsd_monitor_f_zoom_PSD_N
#define PSD_p mccpsd_monitor_f_zoom_PSD_p
#define PSD_p2 mccpsd_monitor_f_zoom_PSD_p2
{   /* Declarations of psd_monitor_f_zoom=PSD_monitor() SETTING parameters. */
int nx = mccpsd_monitor_f_zoom_nx;
int ny = mccpsd_monitor_f_zoom_ny;
char* filename = mccpsd_monitor_f_zoom_filename;
MCNUM xmin = mccpsd_monitor_f_zoom_xmin;
MCNUM xmax = mccpsd_monitor_f_zoom_xmax;
MCNUM ymin = mccpsd_monitor_f_zoom_ymin;
MCNUM ymax = mccpsd_monitor_f_zoom_ymax;
MCNUM xwidth = mccpsd_monitor_f_zoom_xwidth;
MCNUM yheight = mccpsd_monitor_f_zoom_yheight;
MCNUM restore_neutron = mccpsd_monitor_f_zoom_restore_neutron;
#line 122 "/usr/share/mcstas/2.6.1/tools/Python/mcrun/../mccodelib/../../../monitors/PSD_monitor.comp"
{
  destroy_darr2d(PSD_N);
  destroy_darr2d(PSD_p);
  destroy_darr2d(PSD_p2);
}
#line 16831 "./ess_extraction.c"
}   /* End of psd_monitor_f_zoom=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

    if (!mcNCounter[17]) fprintf(stderr, "Warning: No neutron could reach Component[17] psd_monitor_f_zoom\n");
    if (mcAbsorbProp[17]) fprintf(stderr, "Warning: %g events were removed in Component[17] psd_monitor_f_zoom=PSD_monitor()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[17]);
  /* User FINALLY code for component 'psd_monitor_f'. */
  SIG_MESSAGE("psd_monitor_f (Finally)");
#define mccompcurname  psd_monitor_f
#define mccompcurtype  PSD_monitor
#define mccompcurindex 18
#define PSD_N mccpsd_monitor_f_PSD_N
#define PSD_p mccpsd_monitor_f_PSD_p
#define PSD_p2 mccpsd_monitor_f_PSD_p2
{   /* Declarations of psd_monitor_f=PSD_monitor() SETTING parameters. */
int nx = mccpsd_monitor_f_nx;
int ny = mccpsd_monitor_f_ny;
char* filename = mccpsd_monitor_f_filename;
MCNUM xmin = mccpsd_monitor_f_xmin;
MCNUM xmax = mccpsd_monitor_f_xmax;
MCNUM ymin = mccpsd_monitor_f_ymin;
MCNUM ymax = mccpsd_monitor_f_ymax;
MCNUM xwidth = mccpsd_monitor_f_xwidth;
MCNUM yheight = mccpsd_monitor_f_yheight;
MCNUM restore_neutron = mccpsd_monitor_f_restore_neutron;
#line 122 "/usr/share/mcstas/2.6.1/tools/Python/mcrun/../mccodelib/../../../monitors/PSD_monitor.comp"
{
  destroy_darr2d(PSD_N);
  destroy_darr2d(PSD_p);
  destroy_darr2d(PSD_p2);
}
#line 16867 "./ess_extraction.c"
}   /* End of psd_monitor_f=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

    if (!mcNCounter[18]) fprintf(stderr, "Warning: No neutron could reach Component[18] psd_monitor_f\n");
    if (mcAbsorbProp[18]) fprintf(stderr, "Warning: %g events were removed in Component[18] psd_monitor_f=PSD_monitor()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[18]);
    if (!mcNCounter[19]) fprintf(stderr, "Warning: No neutron could reach Component[19] f_divpos\n");
    if (mcAbsorbProp[19]) fprintf(stderr, "Warning: %g events were removed in Component[19] f_divpos=DivPos_monitor()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[19]);
    if (!mcNCounter[20]) fprintf(stderr, "Warning: No neutron could reach Component[20] divhlambda_monitor_f\n");
    if (mcAbsorbProp[20]) fprintf(stderr, "Warning: %g events were removed in Component[20] divhlambda_monitor_f=DivLambda_monitor()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[20]);
  mcsiminfo_close(); 
} /* end finally */
#define magnify mcdis_magnify
#define line mcdis_line
#define dashed_line mcdis_dashed_line
#define multiline mcdis_multiline
#define rectangle mcdis_rectangle
#define box mcdis_box
#define circle mcdis_circle
#define cylinder mcdis_cylinder
#define sphere mcdis_sphere
void mcdisplay(void) {
  printf("MCDISPLAY: start\n");
  /* Components MCDISPLAY code. */

  /* MCDISPLAY code for component 'origin'. */
  SIG_MESSAGE("origin (McDisplay)");
  printf("MCDISPLAY: component %s\n", "origin");
#define mccompcurname  origin
#define mccompcurtype  Progress_bar
#define mccompcurindex 1
#define IntermediateCnts mccorigin_IntermediateCnts
#define StartTime mccorigin_StartTime
#define EndTime mccorigin_EndTime
#define CurrentTime mccorigin_CurrentTime
{   /* Declarations of origin=Progress_bar() SETTING parameters. */
char* profile = mccorigin_profile;
MCNUM percent = mccorigin_percent;
MCNUM flag_save = mccorigin_flag_save;
MCNUM minutes = mccorigin_minutes;
#line 147 "/usr/share/mcstas/2.6.1/tools/Python/mcrun/../mccodelib/../../../misc/Progress_bar.comp"
{
  
}
#line 16916 "./ess_extraction.c"
}   /* End of origin=Progress_bar() SETTING parameter declarations. */
#undef CurrentTime
#undef EndTime
#undef StartTime
#undef IntermediateCnts
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'source'. */
  SIG_MESSAGE("source (McDisplay)");
  printf("MCDISPLAY: component %s\n", "source");
#define mccompcurname  source
#define mccompcurtype  Arm
#define mccompcurindex 2
#line 40 "/usr/share/mcstas/2.6.1/tools/Python/mcrun/../mccodelib/../../../optics/Arm.comp"
{
  /* A bit ugly; hard-coded dimensions. */
  
  line(0,0,0,0.2,0,0);
  line(0,0,0,0,0.2,0);
  line(0,0,0,0,0,0.2);
}
#line 16940 "./ess_extraction.c"
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'source_div'. */
  SIG_MESSAGE("source_div (McDisplay)");
  printf("MCDISPLAY: component %s\n", "source_div");
#define mccompcurname  source_div
#define mccompcurtype  Source_div
#define mccompcurindex 3
#define thetah mccsource_div_thetah
#define thetav mccsource_div_thetav
#define sigmah mccsource_div_sigmah
#define sigmav mccsource_div_sigmav
#define tan_h mccsource_div_tan_h
#define tan_v mccsource_div_tan_v
#define p_init mccsource_div_p_init
#define dist mccsource_div_dist
#define focus_xw mccsource_div_focus_xw
#define focus_yh mccsource_div_focus_yh
{   /* Declarations of source_div=Source_div() SETTING parameters. */
MCNUM xwidth = mccsource_div_xwidth;
MCNUM yheight = mccsource_div_yheight;
MCNUM focus_aw = mccsource_div_focus_aw;
MCNUM focus_ah = mccsource_div_focus_ah;
MCNUM E0 = mccsource_div_E0;
MCNUM dE = mccsource_div_dE;
MCNUM lambda0 = mccsource_div_lambda0;
MCNUM dlambda = mccsource_div_dlambda;
MCNUM gauss = mccsource_div_gauss;
MCNUM flux = mccsource_div_flux;
#line 167 "/usr/share/mcstas/2.6.1/tools/Python/mcrun/../mccodelib/../../../sources/Source_div.comp"
{
  
  multiline(5, -xwidth/2.0, -yheight/2.0, 0.0,
                xwidth/2.0, -yheight/2.0, 0.0,
                xwidth/2.0,  yheight/2.0, 0.0,
               -xwidth/2.0,  yheight/2.0, 0.0,
               -xwidth/2.0, -yheight/2.0, 0.0);
  if (dist) {
    dashed_line(0,0,0, -focus_xw/2,-focus_yh/2,dist, 4);
    dashed_line(0,0,0,  focus_xw/2,-focus_yh/2,dist, 4);
    dashed_line(0,0,0,  focus_xw/2, focus_yh/2,dist, 4);
    dashed_line(0,0,0, -focus_xw/2, focus_yh/2,dist, 4);
  }
}
#line 16987 "./ess_extraction.c"
}   /* End of source_div=Source_div() SETTING parameter declarations. */
#undef focus_yh
#undef focus_xw
#undef dist
#undef p_init
#undef tan_v
#undef tan_h
#undef sigmav
#undef sigmah
#undef thetav
#undef thetah
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'psd_monitor_source'. */
  SIG_MESSAGE("psd_monitor_source (McDisplay)");
  printf("MCDISPLAY: component %s\n", "psd_monitor_source");
#define mccompcurname  psd_monitor_source
#define mccompcurtype  PSD_monitor
#define mccompcurindex 4
#define PSD_N mccpsd_monitor_source_PSD_N
#define PSD_p mccpsd_monitor_source_PSD_p
#define PSD_p2 mccpsd_monitor_source_PSD_p2
{   /* Declarations of psd_monitor_source=PSD_monitor() SETTING parameters. */
int nx = mccpsd_monitor_source_nx;
int ny = mccpsd_monitor_source_ny;
char* filename = mccpsd_monitor_source_filename;
MCNUM xmin = mccpsd_monitor_source_xmin;
MCNUM xmax = mccpsd_monitor_source_xmax;
MCNUM ymin = mccpsd_monitor_source_ymin;
MCNUM ymax = mccpsd_monitor_source_ymax;
MCNUM xwidth = mccpsd_monitor_source_xwidth;
MCNUM yheight = mccpsd_monitor_source_yheight;
MCNUM restore_neutron = mccpsd_monitor_source_restore_neutron;
#line 129 "/usr/share/mcstas/2.6.1/tools/Python/mcrun/../mccodelib/../../../monitors/PSD_monitor.comp"
{
  multiline(5,
    (double)xmin, (double)ymin, 0.0,
    (double)xmax, (double)ymin, 0.0,
    (double)xmax, (double)ymax, 0.0,
    (double)xmin, (double)ymax, 0.0,
    (double)xmin, (double)ymin, 0.0);
}
#line 17032 "./ess_extraction.c"
}   /* End of psd_monitor_source=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'parabolic_optic_before_guide_v'. */
  SIG_MESSAGE("parabolic_optic_before_guide_v (McDisplay)");
  printf("MCDISPLAY: component %s\n", "parabolic_optic_before_guide_v");
#define mccompcurname  parabolic_optic_before_guide_v
#define mccompcurtype  FlatEllipse_finite_mirror
#define mccompcurindex 5
#define reflect mccparabolic_optic_before_guide_v_reflect
#define s mccparabolic_optic_before_guide_v_s
#define pTable mccparabolic_optic_before_guide_v_pTable
#define R0 mccparabolic_optic_before_guide_v_R0
#define Qc mccparabolic_optic_before_guide_v_Qc
#define W mccparabolic_optic_before_guide_v_W
#define alpha mccparabolic_optic_before_guide_v_alpha
#define transmit mccparabolic_optic_before_guide_v_transmit
{   /* Declarations of parabolic_optic_before_guide_v=FlatEllipse_finite_mirror() SETTING parameters. */
MCNUM sourceDist = mccparabolic_optic_before_guide_v_sourceDist;
MCNUM LStart = mccparabolic_optic_before_guide_v_LStart;
MCNUM LEnd = mccparabolic_optic_before_guide_v_LEnd;
MCNUM lStart = mccparabolic_optic_before_guide_v_lStart;
MCNUM lEnd = mccparabolic_optic_before_guide_v_lEnd;
MCNUM r_0 = mccparabolic_optic_before_guide_v_r_0;
MCNUM nummirror = mccparabolic_optic_before_guide_v_nummirror;
MCNUM mf = mccparabolic_optic_before_guide_v_mf;
MCNUM mb = mccparabolic_optic_before_guide_v_mb;
MCNUM mirror_width = mccparabolic_optic_before_guide_v_mirror_width;
MCNUM doubleReflections = mccparabolic_optic_before_guide_v_doubleReflections;
#line 217 "FlatEllipse_finite_mirror.comp"
{
    //Enlarge xy-plane when mcdisplay is ran with --zoom
	magnify("xy");

	//Draw xy-axis contour for Conic Surfaces
	int i;
    for (i = 0; i < s.num_c; i++) {
        double step = (s.c[i].ze-s.c[i].zs)/100;
        double cz;
	    for (cz = s.c[i].zs+step; cz <= s.c[i].ze; cz+= step) {
            double rp = rConic(cz-step,s.c[i]);
            double rc = rConic(cz, s.c[i]);

            line(0,rp,cz-step,0,rc,cz);
            line(0,-rp,cz-step,0,-rc,cz);

            line(rp,0,cz-step,rc,0,cz);
            line(-rp,0,cz-step,-rc,0,cz);
        }
    }

    //Draw xy-axis cross hairs for Disks
    for (i = 0; i < s.num_di; i++) {
        line(s.di[i].r0, 0, s.di[i].z0, s.di[i].r1, 0, s.di[i].z0);
        line(-s.di[i].r0, 0, s.di[i].z0, -s.di[i].r1, 0, s.di[i].z0);
        line(0, s.di[i].r0, s.di[i].z0, 0, s.di[i].r1,s.di[i].z0);
        line(0, -s.di[i].r0, s.di[i].z0, 0, -s.di[i].r1,s.di[i].z0);
    }

}
#line 17098 "./ess_extraction.c"
}   /* End of parabolic_optic_before_guide_v=FlatEllipse_finite_mirror() SETTING parameter declarations. */
#undef transmit
#undef alpha
#undef W
#undef Qc
#undef R0
#undef pTable
#undef s
#undef reflect
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'parabolic_optic_before_guide_h'. */
  SIG_MESSAGE("parabolic_optic_before_guide_h (McDisplay)");
  printf("MCDISPLAY: component %s\n", "parabolic_optic_before_guide_h");
#define mccompcurname  parabolic_optic_before_guide_h
#define mccompcurtype  FlatEllipse_finite_mirror
#define mccompcurindex 6
#define reflect mccparabolic_optic_before_guide_h_reflect
#define s mccparabolic_optic_before_guide_h_s
#define pTable mccparabolic_optic_before_guide_h_pTable
#define R0 mccparabolic_optic_before_guide_h_R0
#define Qc mccparabolic_optic_before_guide_h_Qc
#define W mccparabolic_optic_before_guide_h_W
#define alpha mccparabolic_optic_before_guide_h_alpha
#define transmit mccparabolic_optic_before_guide_h_transmit
{   /* Declarations of parabolic_optic_before_guide_h=FlatEllipse_finite_mirror() SETTING parameters. */
MCNUM sourceDist = mccparabolic_optic_before_guide_h_sourceDist;
MCNUM LStart = mccparabolic_optic_before_guide_h_LStart;
MCNUM LEnd = mccparabolic_optic_before_guide_h_LEnd;
MCNUM lStart = mccparabolic_optic_before_guide_h_lStart;
MCNUM lEnd = mccparabolic_optic_before_guide_h_lEnd;
MCNUM r_0 = mccparabolic_optic_before_guide_h_r_0;
MCNUM nummirror = mccparabolic_optic_before_guide_h_nummirror;
MCNUM mf = mccparabolic_optic_before_guide_h_mf;
MCNUM mb = mccparabolic_optic_before_guide_h_mb;
MCNUM mirror_width = mccparabolic_optic_before_guide_h_mirror_width;
MCNUM doubleReflections = mccparabolic_optic_before_guide_h_doubleReflections;
#line 217 "FlatEllipse_finite_mirror.comp"
{
    //Enlarge xy-plane when mcdisplay is ran with --zoom
	magnify("xy");

	//Draw xy-axis contour for Conic Surfaces
	int i;
    for (i = 0; i < s.num_c; i++) {
        double step = (s.c[i].ze-s.c[i].zs)/100;
        double cz;
	    for (cz = s.c[i].zs+step; cz <= s.c[i].ze; cz+= step) {
            double rp = rConic(cz-step,s.c[i]);
            double rc = rConic(cz, s.c[i]);

            line(0,rp,cz-step,0,rc,cz);
            line(0,-rp,cz-step,0,-rc,cz);

            line(rp,0,cz-step,rc,0,cz);
            line(-rp,0,cz-step,-rc,0,cz);
        }
    }

    //Draw xy-axis cross hairs for Disks
    for (i = 0; i < s.num_di; i++) {
        line(s.di[i].r0, 0, s.di[i].z0, s.di[i].r1, 0, s.di[i].z0);
        line(-s.di[i].r0, 0, s.di[i].z0, -s.di[i].r1, 0, s.di[i].z0);
        line(0, s.di[i].r0, s.di[i].z0, 0, s.di[i].r1,s.di[i].z0);
        line(0, -s.di[i].r0, s.di[i].z0, 0, -s.di[i].r1,s.di[i].z0);
    }

}
#line 17169 "./ess_extraction.c"
}   /* End of parabolic_optic_before_guide_h=FlatEllipse_finite_mirror() SETTING parameter declarations. */
#undef transmit
#undef alpha
#undef W
#undef Qc
#undef R0
#undef pTable
#undef s
#undef reflect
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'after_optic_source'. */
  SIG_MESSAGE("after_optic_source (McDisplay)");
  printf("MCDISPLAY: component %s\n", "after_optic_source");
#define mccompcurname  after_optic_source
#define mccompcurtype  Arm
#define mccompcurindex 7
#line 40 "/usr/share/mcstas/2.6.1/tools/Python/mcrun/../mccodelib/../../../optics/Arm.comp"
{
  /* A bit ugly; hard-coded dimensions. */
  
  line(0,0,0,0.2,0,0);
  line(0,0,0,0,0.2,0);
  line(0,0,0,0,0,0.2);
}
#line 17197 "./ess_extraction.c"
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'psd_monitor_afteropticsource'. */
  SIG_MESSAGE("psd_monitor_afteropticsource (McDisplay)");
  printf("MCDISPLAY: component %s\n", "psd_monitor_afteropticsource");
#define mccompcurname  psd_monitor_afteropticsource
#define mccompcurtype  PSD_monitor
#define mccompcurindex 8
#define PSD_N mccpsd_monitor_afteropticsource_PSD_N
#define PSD_p mccpsd_monitor_afteropticsource_PSD_p
#define PSD_p2 mccpsd_monitor_afteropticsource_PSD_p2
{   /* Declarations of psd_monitor_afteropticsource=PSD_monitor() SETTING parameters. */
int nx = mccpsd_monitor_afteropticsource_nx;
int ny = mccpsd_monitor_afteropticsource_ny;
char* filename = mccpsd_monitor_afteropticsource_filename;
MCNUM xmin = mccpsd_monitor_afteropticsource_xmin;
MCNUM xmax = mccpsd_monitor_afteropticsource_xmax;
MCNUM ymin = mccpsd_monitor_afteropticsource_ymin;
MCNUM ymax = mccpsd_monitor_afteropticsource_ymax;
MCNUM xwidth = mccpsd_monitor_afteropticsource_xwidth;
MCNUM yheight = mccpsd_monitor_afteropticsource_yheight;
MCNUM restore_neutron = mccpsd_monitor_afteropticsource_restore_neutron;
#line 129 "/usr/share/mcstas/2.6.1/tools/Python/mcrun/../mccodelib/../../../monitors/PSD_monitor.comp"
{
  multiline(5,
    (double)xmin, (double)ymin, 0.0,
    (double)xmax, (double)ymin, 0.0,
    (double)xmax, (double)ymax, 0.0,
    (double)xmin, (double)ymax, 0.0,
    (double)xmin, (double)ymin, 0.0);
}
#line 17231 "./ess_extraction.c"
}   /* End of psd_monitor_afteropticsource=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'guide_gravity_1'. */
  SIG_MESSAGE("guide_gravity_1 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "guide_gravity_1");
#define mccompcurname  guide_gravity_1
#define mccompcurtype  Guide_gravity
#define mccompcurindex 9
#define GVars mccguide_gravity_1_GVars
#define pTable mccguide_gravity_1_pTable
{   /* Declarations of guide_gravity_1=Guide_gravity() SETTING parameters. */
MCNUM w1 = mccguide_gravity_1_w1;
MCNUM h1 = mccguide_gravity_1_h1;
MCNUM w2 = mccguide_gravity_1_w2;
MCNUM h2 = mccguide_gravity_1_h2;
MCNUM l = mccguide_gravity_1_l;
MCNUM R0 = mccguide_gravity_1_R0;
MCNUM Qc = mccguide_gravity_1_Qc;
MCNUM alpha = mccguide_gravity_1_alpha;
MCNUM m = mccguide_gravity_1_m;
MCNUM W = mccguide_gravity_1_W;
MCNUM nslit = mccguide_gravity_1_nslit;
MCNUM d = mccguide_gravity_1_d;
MCNUM mleft = mccguide_gravity_1_mleft;
MCNUM mright = mccguide_gravity_1_mright;
MCNUM mtop = mccguide_gravity_1_mtop;
MCNUM mbottom = mccguide_gravity_1_mbottom;
MCNUM nhslit = mccguide_gravity_1_nhslit;
MCNUM G = mccguide_gravity_1_G;
MCNUM aleft = mccguide_gravity_1_aleft;
MCNUM aright = mccguide_gravity_1_aright;
MCNUM atop = mccguide_gravity_1_atop;
MCNUM abottom = mccguide_gravity_1_abottom;
MCNUM wavy = mccguide_gravity_1_wavy;
MCNUM wavy_z = mccguide_gravity_1_wavy_z;
MCNUM wavy_tb = mccguide_gravity_1_wavy_tb;
MCNUM wavy_lr = mccguide_gravity_1_wavy_lr;
MCNUM chamfers = mccguide_gravity_1_chamfers;
MCNUM chamfers_z = mccguide_gravity_1_chamfers_z;
MCNUM chamfers_lr = mccguide_gravity_1_chamfers_lr;
MCNUM chamfers_tb = mccguide_gravity_1_chamfers_tb;
MCNUM nelements = mccguide_gravity_1_nelements;
MCNUM nu = mccguide_gravity_1_nu;
MCNUM phase = mccguide_gravity_1_phase;
char* reflect = mccguide_gravity_1_reflect;
#line 571 "/usr/share/mcstas/2.6.1/tools/Python/mcrun/../mccodelib/../../../optics/Guide_gravity.comp"
{

  if (l > 0 && nelements > 0) {
    int i,j,n;
    double x1,x2,x3,x4;
    double y1,y2,y3,y4;
    double nel = (nelements > 11 ? 11 : nelements);


    for (n=0; n<nel; n++)
    {
      double z0, z1;
      z0 =     n*(l/nel);
      z1 = (n+1)*(l/nel);

      for(j = 0; j < nhslit; j++)
      {
        y1 = j*(GVars.h1c+d)         - h1/2.0;
        y2 = j*(GVars.h2c+d)         - h2/2.0;
        y3 = (j+1)*(GVars.h1c+d) - d - h1/2.0;
        y4 = (j+1)*(GVars.h2c+d) - d - h2/2.0;
        for(i = 0; i < nslit; i++)
        {
          x1 = i*(GVars.w1c+d)         - w1/2.0;
          x2 = i*(GVars.w2c+d)         - w2/2.0;
          x3 = (i+1)*(GVars.w1c+d) - d - w1/2.0;
          x4 = (i+1)*(GVars.w2c+d) - d - w2/2.0;
          multiline(5,
                    x1, y1, z0,
                    x2, y2, z1,
                    x2, y4, z1,
                    x1, y3, z0,
                    x1, y1, z0);
          multiline(5,
                    x3, y1, z0,
                    x4, y2, z1,
                    x4, y4, z1,
                    x3, y3, z0,
                    x3, y1, z0);
        }
        line(-w1/2.0, y1, z0, w1/2.0, y1, z0);
        line(-w2/2.0, y2, z1, w2/2.0, y2, z1);
      }
    }

    if (nu || phase) {
      double radius = sqrt(w1*w1+l*l);
      /* cylinder top/center/bottom  */
      circle("xz", 0,-h1/2,l/2,radius);
      circle("xz", 0,0    ,l/2,radius);
      circle("xz", 0, h1/2,l/2,radius);
    }
  }
  else {
    /* A bit ugly; hard-coded dimensions. */

    line(0,0,0,0.2,0,0);
    line(0,0,0,0,0.2,0);
    line(0,0,0,0,0,0.2);
  }

}
#line 17346 "./ess_extraction.c"
}   /* End of guide_gravity_1=Guide_gravity() SETTING parameter declarations. */
#undef pTable
#undef GVars
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'monitor_1'. */
  SIG_MESSAGE("monitor_1 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "monitor_1");
#define mccompcurname  monitor_1
#define mccompcurtype  Arm
#define mccompcurindex 10
#line 40 "/usr/share/mcstas/2.6.1/tools/Python/mcrun/../mccodelib/../../../optics/Arm.comp"
{
  /* A bit ugly; hard-coded dimensions. */
  
  line(0,0,0,0.2,0,0);
  line(0,0,0,0,0.2,0);
  line(0,0,0,0,0,0.2);
}
#line 17368 "./ess_extraction.c"
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'divpos_monitor'. */
  SIG_MESSAGE("divpos_monitor (McDisplay)");
  printf("MCDISPLAY: component %s\n", "divpos_monitor");
#define mccompcurname  divpos_monitor
#define mccompcurtype  DivPos_monitor
#define mccompcurindex 11
#define nh mccdivpos_monitor_nh
#define ndiv mccdivpos_monitor_ndiv
#define Div_N mccdivpos_monitor_Div_N
#define Div_p mccdivpos_monitor_Div_p
#define Div_p2 mccdivpos_monitor_Div_p2
{   /* Declarations of divpos_monitor=DivPos_monitor() SETTING parameters. */
char* filename = mccdivpos_monitor_filename;
MCNUM xmin = mccdivpos_monitor_xmin;
MCNUM xmax = mccdivpos_monitor_xmax;
MCNUM ymin = mccdivpos_monitor_ymin;
MCNUM ymax = mccdivpos_monitor_ymax;
MCNUM xwidth = mccdivpos_monitor_xwidth;
MCNUM yheight = mccdivpos_monitor_yheight;
MCNUM maxdiv_h = mccdivpos_monitor_maxdiv_h;
MCNUM restore_neutron = mccdivpos_monitor_restore_neutron;
MCNUM nx = mccdivpos_monitor_nx;
MCNUM ny = mccdivpos_monitor_ny;
MCNUM nz = mccdivpos_monitor_nz;
int nowritefile = mccdivpos_monitor_nowritefile;
#line 134 "/usr/share/mcstas/2.6.1/tools/Python/mcrun/../mccodelib/../../../monitors/DivPos_monitor.comp"
{
    
    multiline(5, (double)xmin, (double)ymin, 0.0,
                 (double)xmax, (double)ymin, 0.0,
                 (double)xmax, (double)ymax, 0.0,
                 (double)xmin, (double)ymax, 0.0,
                 (double)xmin, (double)ymin, 0.0);
}
#line 17407 "./ess_extraction.c"
}   /* End of divpos_monitor=DivPos_monitor() SETTING parameter declarations. */
#undef Div_p2
#undef Div_p
#undef Div_N
#undef ndiv
#undef nh
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'psd_monitor_g1'. */
  SIG_MESSAGE("psd_monitor_g1 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "psd_monitor_g1");
#define mccompcurname  psd_monitor_g1
#define mccompcurtype  PSD_monitor
#define mccompcurindex 12
#define PSD_N mccpsd_monitor_g1_PSD_N
#define PSD_p mccpsd_monitor_g1_PSD_p
#define PSD_p2 mccpsd_monitor_g1_PSD_p2
{   /* Declarations of psd_monitor_g1=PSD_monitor() SETTING parameters. */
int nx = mccpsd_monitor_g1_nx;
int ny = mccpsd_monitor_g1_ny;
char* filename = mccpsd_monitor_g1_filename;
MCNUM xmin = mccpsd_monitor_g1_xmin;
MCNUM xmax = mccpsd_monitor_g1_xmax;
MCNUM ymin = mccpsd_monitor_g1_ymin;
MCNUM ymax = mccpsd_monitor_g1_ymax;
MCNUM xwidth = mccpsd_monitor_g1_xwidth;
MCNUM yheight = mccpsd_monitor_g1_yheight;
MCNUM restore_neutron = mccpsd_monitor_g1_restore_neutron;
#line 129 "/usr/share/mcstas/2.6.1/tools/Python/mcrun/../mccodelib/../../../monitors/PSD_monitor.comp"
{
  multiline(5,
    (double)xmin, (double)ymin, 0.0,
    (double)xmax, (double)ymin, 0.0,
    (double)xmax, (double)ymax, 0.0,
    (double)xmin, (double)ymax, 0.0,
    (double)xmin, (double)ymin, 0.0);
}
#line 17447 "./ess_extraction.c"
}   /* End of psd_monitor_g1=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'divhlambda_monitor_g1'. */
  SIG_MESSAGE("divhlambda_monitor_g1 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "divhlambda_monitor_g1");
#define mccompcurname  divhlambda_monitor_g1
#define mccompcurtype  DivLambda_monitor
#define mccompcurindex 13
#define nL mccdivhlambda_monitor_g1_nL
#define nh mccdivhlambda_monitor_g1_nh
#define Div_N mccdivhlambda_monitor_g1_Div_N
#define Div_p mccdivhlambda_monitor_g1_Div_p
#define Div_p2 mccdivhlambda_monitor_g1_Div_p2
{   /* Declarations of divhlambda_monitor_g1=DivLambda_monitor() SETTING parameters. */
char* filename = mccdivhlambda_monitor_g1_filename;
MCNUM xmin = mccdivhlambda_monitor_g1_xmin;
MCNUM xmax = mccdivhlambda_monitor_g1_xmax;
MCNUM ymin = mccdivhlambda_monitor_g1_ymin;
MCNUM ymax = mccdivhlambda_monitor_g1_ymax;
MCNUM xwidth = mccdivhlambda_monitor_g1_xwidth;
MCNUM yheight = mccdivhlambda_monitor_g1_yheight;
MCNUM maxdiv_h = mccdivhlambda_monitor_g1_maxdiv_h;
MCNUM Lmin = mccdivhlambda_monitor_g1_Lmin;
MCNUM Lmax = mccdivhlambda_monitor_g1_Lmax;
MCNUM restore_neutron = mccdivhlambda_monitor_g1_restore_neutron;
MCNUM nx = mccdivhlambda_monitor_g1_nx;
MCNUM ny = mccdivhlambda_monitor_g1_ny;
MCNUM nz = mccdivhlambda_monitor_g1_nz;
int nowritefile = mccdivhlambda_monitor_g1_nowritefile;
#line 133 "/usr/share/mcstas/2.6.1/tools/Python/mcrun/../mccodelib/../../../monitors/DivLambda_monitor.comp"
{

    multiline(5, (double)xmin, (double)ymin, 0.0,
                 (double)xmax, (double)ymin, 0.0,
                 (double)xmax, (double)ymax, 0.0,
                 (double)xmin, (double)ymax, 0.0,
                 (double)xmin, (double)ymin, 0.0);
}
#line 17492 "./ess_extraction.c"
}   /* End of divhlambda_monitor_g1=DivLambda_monitor() SETTING parameter declarations. */
#undef Div_p2
#undef Div_p
#undef Div_N
#undef nh
#undef nL
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'parabolic_optic1'. */
  SIG_MESSAGE("parabolic_optic1 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "parabolic_optic1");
#define mccompcurname  parabolic_optic1
#define mccompcurtype  FlatEllipse_finite_mirror
#define mccompcurindex 14
#define reflect mccparabolic_optic1_reflect
#define s mccparabolic_optic1_s
#define pTable mccparabolic_optic1_pTable
#define R0 mccparabolic_optic1_R0
#define Qc mccparabolic_optic1_Qc
#define W mccparabolic_optic1_W
#define alpha mccparabolic_optic1_alpha
#define transmit mccparabolic_optic1_transmit
{   /* Declarations of parabolic_optic1=FlatEllipse_finite_mirror() SETTING parameters. */
MCNUM sourceDist = mccparabolic_optic1_sourceDist;
MCNUM LStart = mccparabolic_optic1_LStart;
MCNUM LEnd = mccparabolic_optic1_LEnd;
MCNUM lStart = mccparabolic_optic1_lStart;
MCNUM lEnd = mccparabolic_optic1_lEnd;
MCNUM r_0 = mccparabolic_optic1_r_0;
MCNUM nummirror = mccparabolic_optic1_nummirror;
MCNUM mf = mccparabolic_optic1_mf;
MCNUM mb = mccparabolic_optic1_mb;
MCNUM mirror_width = mccparabolic_optic1_mirror_width;
MCNUM doubleReflections = mccparabolic_optic1_doubleReflections;
#line 217 "FlatEllipse_finite_mirror.comp"
{
    //Enlarge xy-plane when mcdisplay is ran with --zoom
	magnify("xy");

	//Draw xy-axis contour for Conic Surfaces
	int i;
    for (i = 0; i < s.num_c; i++) {
        double step = (s.c[i].ze-s.c[i].zs)/100;
        double cz;
	    for (cz = s.c[i].zs+step; cz <= s.c[i].ze; cz+= step) {
            double rp = rConic(cz-step,s.c[i]);
            double rc = rConic(cz, s.c[i]);

            line(0,rp,cz-step,0,rc,cz);
            line(0,-rp,cz-step,0,-rc,cz);

            line(rp,0,cz-step,rc,0,cz);
            line(-rp,0,cz-step,-rc,0,cz);
        }
    }

    //Draw xy-axis cross hairs for Disks
    for (i = 0; i < s.num_di; i++) {
        line(s.di[i].r0, 0, s.di[i].z0, s.di[i].r1, 0, s.di[i].z0);
        line(-s.di[i].r0, 0, s.di[i].z0, -s.di[i].r1, 0, s.di[i].z0);
        line(0, s.di[i].r0, s.di[i].z0, 0, s.di[i].r1,s.di[i].z0);
        line(0, -s.di[i].r0, s.di[i].z0, 0, -s.di[i].r1,s.di[i].z0);
    }

}
#line 17560 "./ess_extraction.c"
}   /* End of parabolic_optic1=FlatEllipse_finite_mirror() SETTING parameter declarations. */
#undef transmit
#undef alpha
#undef W
#undef Qc
#undef R0
#undef pTable
#undef s
#undef reflect
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'parabolic_optic2'. */
  SIG_MESSAGE("parabolic_optic2 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "parabolic_optic2");
#define mccompcurname  parabolic_optic2
#define mccompcurtype  FlatEllipse_finite_mirror
#define mccompcurindex 15
#define reflect mccparabolic_optic2_reflect
#define s mccparabolic_optic2_s
#define pTable mccparabolic_optic2_pTable
#define R0 mccparabolic_optic2_R0
#define Qc mccparabolic_optic2_Qc
#define W mccparabolic_optic2_W
#define alpha mccparabolic_optic2_alpha
#define transmit mccparabolic_optic2_transmit
{   /* Declarations of parabolic_optic2=FlatEllipse_finite_mirror() SETTING parameters. */
MCNUM sourceDist = mccparabolic_optic2_sourceDist;
MCNUM LStart = mccparabolic_optic2_LStart;
MCNUM LEnd = mccparabolic_optic2_LEnd;
MCNUM lStart = mccparabolic_optic2_lStart;
MCNUM lEnd = mccparabolic_optic2_lEnd;
MCNUM r_0 = mccparabolic_optic2_r_0;
MCNUM nummirror = mccparabolic_optic2_nummirror;
MCNUM mf = mccparabolic_optic2_mf;
MCNUM mb = mccparabolic_optic2_mb;
MCNUM mirror_width = mccparabolic_optic2_mirror_width;
MCNUM doubleReflections = mccparabolic_optic2_doubleReflections;
#line 217 "FlatEllipse_finite_mirror.comp"
{
    //Enlarge xy-plane when mcdisplay is ran with --zoom
	magnify("xy");

	//Draw xy-axis contour for Conic Surfaces
	int i;
    for (i = 0; i < s.num_c; i++) {
        double step = (s.c[i].ze-s.c[i].zs)/100;
        double cz;
	    for (cz = s.c[i].zs+step; cz <= s.c[i].ze; cz+= step) {
            double rp = rConic(cz-step,s.c[i]);
            double rc = rConic(cz, s.c[i]);

            line(0,rp,cz-step,0,rc,cz);
            line(0,-rp,cz-step,0,-rc,cz);

            line(rp,0,cz-step,rc,0,cz);
            line(-rp,0,cz-step,-rc,0,cz);
        }
    }

    //Draw xy-axis cross hairs for Disks
    for (i = 0; i < s.num_di; i++) {
        line(s.di[i].r0, 0, s.di[i].z0, s.di[i].r1, 0, s.di[i].z0);
        line(-s.di[i].r0, 0, s.di[i].z0, -s.di[i].r1, 0, s.di[i].z0);
        line(0, s.di[i].r0, s.di[i].z0, 0, s.di[i].r1,s.di[i].z0);
        line(0, -s.di[i].r0, s.di[i].z0, 0, -s.di[i].r1,s.di[i].z0);
    }

}
#line 17631 "./ess_extraction.c"
}   /* End of parabolic_optic2=FlatEllipse_finite_mirror() SETTING parameter declarations. */
#undef transmit
#undef alpha
#undef W
#undef Qc
#undef R0
#undef pTable
#undef s
#undef reflect
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'monitor_2'. */
  SIG_MESSAGE("monitor_2 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "monitor_2");
#define mccompcurname  monitor_2
#define mccompcurtype  Arm
#define mccompcurindex 16
#line 40 "/usr/share/mcstas/2.6.1/tools/Python/mcrun/../mccodelib/../../../optics/Arm.comp"
{
  /* A bit ugly; hard-coded dimensions. */
  
  line(0,0,0,0.2,0,0);
  line(0,0,0,0,0.2,0);
  line(0,0,0,0,0,0.2);
}
#line 17659 "./ess_extraction.c"
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'psd_monitor_f_zoom'. */
  SIG_MESSAGE("psd_monitor_f_zoom (McDisplay)");
  printf("MCDISPLAY: component %s\n", "psd_monitor_f_zoom");
#define mccompcurname  psd_monitor_f_zoom
#define mccompcurtype  PSD_monitor
#define mccompcurindex 17
#define PSD_N mccpsd_monitor_f_zoom_PSD_N
#define PSD_p mccpsd_monitor_f_zoom_PSD_p
#define PSD_p2 mccpsd_monitor_f_zoom_PSD_p2
{   /* Declarations of psd_monitor_f_zoom=PSD_monitor() SETTING parameters. */
int nx = mccpsd_monitor_f_zoom_nx;
int ny = mccpsd_monitor_f_zoom_ny;
char* filename = mccpsd_monitor_f_zoom_filename;
MCNUM xmin = mccpsd_monitor_f_zoom_xmin;
MCNUM xmax = mccpsd_monitor_f_zoom_xmax;
MCNUM ymin = mccpsd_monitor_f_zoom_ymin;
MCNUM ymax = mccpsd_monitor_f_zoom_ymax;
MCNUM xwidth = mccpsd_monitor_f_zoom_xwidth;
MCNUM yheight = mccpsd_monitor_f_zoom_yheight;
MCNUM restore_neutron = mccpsd_monitor_f_zoom_restore_neutron;
#line 129 "/usr/share/mcstas/2.6.1/tools/Python/mcrun/../mccodelib/../../../monitors/PSD_monitor.comp"
{
  multiline(5,
    (double)xmin, (double)ymin, 0.0,
    (double)xmax, (double)ymin, 0.0,
    (double)xmax, (double)ymax, 0.0,
    (double)xmin, (double)ymax, 0.0,
    (double)xmin, (double)ymin, 0.0);
}
#line 17693 "./ess_extraction.c"
}   /* End of psd_monitor_f_zoom=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'psd_monitor_f'. */
  SIG_MESSAGE("psd_monitor_f (McDisplay)");
  printf("MCDISPLAY: component %s\n", "psd_monitor_f");
#define mccompcurname  psd_monitor_f
#define mccompcurtype  PSD_monitor
#define mccompcurindex 18
#define PSD_N mccpsd_monitor_f_PSD_N
#define PSD_p mccpsd_monitor_f_PSD_p
#define PSD_p2 mccpsd_monitor_f_PSD_p2
{   /* Declarations of psd_monitor_f=PSD_monitor() SETTING parameters. */
int nx = mccpsd_monitor_f_nx;
int ny = mccpsd_monitor_f_ny;
char* filename = mccpsd_monitor_f_filename;
MCNUM xmin = mccpsd_monitor_f_xmin;
MCNUM xmax = mccpsd_monitor_f_xmax;
MCNUM ymin = mccpsd_monitor_f_ymin;
MCNUM ymax = mccpsd_monitor_f_ymax;
MCNUM xwidth = mccpsd_monitor_f_xwidth;
MCNUM yheight = mccpsd_monitor_f_yheight;
MCNUM restore_neutron = mccpsd_monitor_f_restore_neutron;
#line 129 "/usr/share/mcstas/2.6.1/tools/Python/mcrun/../mccodelib/../../../monitors/PSD_monitor.comp"
{
  multiline(5,
    (double)xmin, (double)ymin, 0.0,
    (double)xmax, (double)ymin, 0.0,
    (double)xmax, (double)ymax, 0.0,
    (double)xmin, (double)ymax, 0.0,
    (double)xmin, (double)ymin, 0.0);
}
#line 17731 "./ess_extraction.c"
}   /* End of psd_monitor_f=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'f_divpos'. */
  SIG_MESSAGE("f_divpos (McDisplay)");
  printf("MCDISPLAY: component %s\n", "f_divpos");
#define mccompcurname  f_divpos
#define mccompcurtype  DivPos_monitor
#define mccompcurindex 19
#define nh mccf_divpos_nh
#define ndiv mccf_divpos_ndiv
#define Div_N mccf_divpos_Div_N
#define Div_p mccf_divpos_Div_p
#define Div_p2 mccf_divpos_Div_p2
{   /* Declarations of f_divpos=DivPos_monitor() SETTING parameters. */
char* filename = mccf_divpos_filename;
MCNUM xmin = mccf_divpos_xmin;
MCNUM xmax = mccf_divpos_xmax;
MCNUM ymin = mccf_divpos_ymin;
MCNUM ymax = mccf_divpos_ymax;
MCNUM xwidth = mccf_divpos_xwidth;
MCNUM yheight = mccf_divpos_yheight;
MCNUM maxdiv_h = mccf_divpos_maxdiv_h;
MCNUM restore_neutron = mccf_divpos_restore_neutron;
MCNUM nx = mccf_divpos_nx;
MCNUM ny = mccf_divpos_ny;
MCNUM nz = mccf_divpos_nz;
int nowritefile = mccf_divpos_nowritefile;
#line 134 "/usr/share/mcstas/2.6.1/tools/Python/mcrun/../mccodelib/../../../monitors/DivPos_monitor.comp"
{
    
    multiline(5, (double)xmin, (double)ymin, 0.0,
                 (double)xmax, (double)ymin, 0.0,
                 (double)xmax, (double)ymax, 0.0,
                 (double)xmin, (double)ymax, 0.0,
                 (double)xmin, (double)ymin, 0.0);
}
#line 17774 "./ess_extraction.c"
}   /* End of f_divpos=DivPos_monitor() SETTING parameter declarations. */
#undef Div_p2
#undef Div_p
#undef Div_N
#undef ndiv
#undef nh
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'divhlambda_monitor_f'. */
  SIG_MESSAGE("divhlambda_monitor_f (McDisplay)");
  printf("MCDISPLAY: component %s\n", "divhlambda_monitor_f");
#define mccompcurname  divhlambda_monitor_f
#define mccompcurtype  DivLambda_monitor
#define mccompcurindex 20
#define nL mccdivhlambda_monitor_f_nL
#define nh mccdivhlambda_monitor_f_nh
#define Div_N mccdivhlambda_monitor_f_Div_N
#define Div_p mccdivhlambda_monitor_f_Div_p
#define Div_p2 mccdivhlambda_monitor_f_Div_p2
{   /* Declarations of divhlambda_monitor_f=DivLambda_monitor() SETTING parameters. */
char* filename = mccdivhlambda_monitor_f_filename;
MCNUM xmin = mccdivhlambda_monitor_f_xmin;
MCNUM xmax = mccdivhlambda_monitor_f_xmax;
MCNUM ymin = mccdivhlambda_monitor_f_ymin;
MCNUM ymax = mccdivhlambda_monitor_f_ymax;
MCNUM xwidth = mccdivhlambda_monitor_f_xwidth;
MCNUM yheight = mccdivhlambda_monitor_f_yheight;
MCNUM maxdiv_h = mccdivhlambda_monitor_f_maxdiv_h;
MCNUM Lmin = mccdivhlambda_monitor_f_Lmin;
MCNUM Lmax = mccdivhlambda_monitor_f_Lmax;
MCNUM restore_neutron = mccdivhlambda_monitor_f_restore_neutron;
MCNUM nx = mccdivhlambda_monitor_f_nx;
MCNUM ny = mccdivhlambda_monitor_f_ny;
MCNUM nz = mccdivhlambda_monitor_f_nz;
int nowritefile = mccdivhlambda_monitor_f_nowritefile;
#line 133 "/usr/share/mcstas/2.6.1/tools/Python/mcrun/../mccodelib/../../../monitors/DivLambda_monitor.comp"
{

    multiline(5, (double)xmin, (double)ymin, 0.0,
                 (double)xmax, (double)ymin, 0.0,
                 (double)xmax, (double)ymax, 0.0,
                 (double)xmin, (double)ymax, 0.0,
                 (double)xmin, (double)ymin, 0.0);
}
#line 17821 "./ess_extraction.c"
}   /* End of divhlambda_monitor_f=DivLambda_monitor() SETTING parameter declarations. */
#undef Div_p2
#undef Div_p
#undef Div_N
#undef nh
#undef nL
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  printf("MCDISPLAY: end\n");
} /* end display */
#undef magnify
#undef line
#undef dashed_line
#undef multiline
#undef rectangle
#undef box
#undef circle
#undef cylinder
#undef sphere
/* end of generated C code ./ess_extraction.c */
