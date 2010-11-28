/* ------- file: -------------------------- geometry.h --------------

       Version:       rh2.0, 1-D plane-parallel
       Author:        Tiago Pereira (tiago.pereira@nasa.gov)
       Last modified: Tue Nov 23 22:09:58 2010 --

       --------------------------                      ----------RH-- */

#ifndef __GEOMETRY_H__
#define __GEOMETRY_H__

/* --- Define geometric quantities for 1-D plane-parallel version --  */


enum boundcond  {ZERO, THERMALIZED, IRRADIATED, REFLECTIVE};
enum mass_scale {GEOMETRIC, COLUMN_MASS, TAU500};
enum vertical   {TOP=0, BOTTOM};

typedef struct {
  enum     mass_scale  scale;
  enum     boundcond vboundary[2]; 
  int      Ndep, Nrays;
  double  *height, *cmass, *tau_ref, *mux, *muy, *muz, *wmu, *vel,
         **Itop, **Ibottom;
} Geometry;


/* For the NetCDF input file */
typedef struct {
  char     *file_name;
  int       ncid,    nx_id,    ny_id,    nz_id,   nhyd_id;
  int       T_varid, ne_varid, vz_varid, nh_varid; 
  size_t    nx,      ny,       nz,       NHydr;
} NCDF_Atmos_file;


/* --- Associated function prototypes --               -------------- */

void convertScales(Atmosphere *atmos, Geometry *geometry);
void getAngleQuad(Geometry *geometry);
void getBoundary(Geometry *geometry);
void MULTIatmos(Atmosphere *atmos, Geometry *geometry);
void writeGeometry(Geometry *geometry);

void init_ncdf(Atmosphere *atmos, Geometry *geometry, NCDF_Atmos_file *infile);
void readAtmos_ncdf(int xi, int yi, Atmosphere *atmos, Geometry *geometry,
		    NCDF_Atmos_file *infile);
void close_ncdf(Atmosphere *atmos, Geometry *geometry, NCDF_Atmos_file *infile);


/* --- Formal solution related --                      -------------- */

double Feautrier(int nspect, int mu, double *chi, double *S,
		 enum FeautrierOrder order, double *P, double *Psi);
void Piecewise_1D(int nspect, int mu, bool_t to_obs,
		  double *chi, double *S, double *I, double *Psi);
void PiecewiseStokes(int nspect, int mu, bool_t to_obs,
		     double *chi_I, double **S, double **I, double *Psi);


#endif /* !__GEOMETRY_H__ */

/* ---------------------------------------- geometry.h -------------- */
