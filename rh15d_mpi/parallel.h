/* ------- file: -------------------------- parallel.h --------------

       Version:       rh2.0, 1-D plane-parallel
       Author:        Tiago Pereira (tiago.pereira@nasa.gov)
       Last modified: Tue Nov 24 10:59:58 2010 --

       --------------------------                      ----------RH-- */

#ifndef __PARALLEL_H__
#define __PARALLEL_H__

#include <mpi.h>

typedef struct {
  int  size, rank, namelen, Ntasks, task, nx, ny;
  int *xnum, *ynum; 
  long backgrrecno; /* record number for background file */
  char name[MPI_MAX_PROCESSOR_NAME];
  MPI_Comm comm;
  MPI_Info info;
} MPI_data;

void init_Background(void);
void Background_p(bool_t analyzeoutput, bool_t equilibria_only);
void close_Background();
void writeBRS_p(void);
void readBRS_p(void);
void init_ncdf_BRS(void);
void close_ncdf_BRS(void);
void writeBRS_ncdf(void);
void init_ncdf_J(void);
void close_ncdf_J(void);
void writeJlambda_ncdf(int nspect, double *J);
void writeJ20_ncdf(int nspect, double *J);
void readJlambda_ncdf(int nspect, double *J);
void readJ20_ncdf(int nspect, double *J);
void initSolution_p(void);

#define ERR(e,r) {printf("(EEE) %s: NetCDF: %s\n", r, nc_strerror(e)); exit(EXIT_FAILURE);}

#endif /* !__PARALLEL_H__ */

/* ---------------------------------------- parallel.h -------------- */
