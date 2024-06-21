/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(obmd,FixObmd)

#else

#ifndef LMP_FIX_OBMD_H
#define LMP_FIX_OBMD_H

#include "fix.h"

namespace LAMMPS_NS {

class FixObmd : public Fix {
 public:
  FixObmd(class LAMMPS *, int, char **);
  ~FixObmd();
  int setmask();
  void init();
  void pre_exchange();
  void post_force(int);
  void setup(int);
  void min_setup(int);
  void min_post_force(int);
  void write_restart(FILE *);
  void restart(char *);
  void *extract(const char *, int &);
  virtual double compute_scalar();
  virtual double compute_vector(int);
 private:
  int ninsert,ntype,nfreq,seed;
  class Region *iregion,*iregion2,*iregion3,*iregion4,*iregion5,*iregion6;
  int globalflag,localflag,maxattempt,rateflag,scaleflag,targetflag;
  int mode,rigidflag,shakeflag,idnext,distflag,orientflag,molflag,usherflag,nearflag;
  double lo,hi,deltasq,nearsq,rate,sigma;
  double vxlo,vxhi,vylo,vyhi,vzlo,vzhi;
  double xlo,xhi,ylo,yhi,zlo,zhi,xmid,ymid,zmid;
  double rx,ry,rz,tx,ty,tz;
  char *idregion;
  char *idregion2;
  char *idregion3;
  char *idregion4;
  char *idregion5;
  char *idregion6;
  char *idrigid,*idshake;

  double tau, alpha, buffer_size, shear_size, nbuf;
  double pxx,pxy,pxz,dpxx,t0_left,t0_right,lambda;
  double mol_len;

  void check_ghosts();
  bigint lastcheck;

  double vcml[3];
  double vcmr[3];

  double vnewl[3];
  double vnewr[3];
  double vnewl_all[3];
  double vnewr_all[3];
  double enl, enl_all;
  double enr, enr_all;
  double energyForce_left[1];
  double energyForce_right[1];


  int me;
  class Molecule **onemols;
  int nmol,natom_max;
  double *molfrac;
  double **coords;
  imageint *imageflags;
  class Fix *fixrigid,*fixshake;
  double oneradius;

  int xstyle,xvar;
  char *xstr;
  int ystyle,yvar;
  char *ystr;
  int zstyle,zvar;
  char *zstr;
  int qstyle,qvar;
  char *qstr;
  int t0_left_style,t0_left_var;
  char *t0_left_str;
  int t0_right_style,t0_right_var;
  char *t0_right_str;

  int nfirst,ninserted;
  tagint maxtag_all,maxmol_all;
  class RanPark *random;

  void find_maxid();
  void options(int, char **);
  double momentumForce_left[3];
  double momentumForce_right[3];
  int *list, *mark;
  int ndeleted, ncount2, nmax, nmax2;

  void add_energy_forces(int, Region*, double *fac, double *vcm);
  void try_inserting(Region*,int, double *vnewl, double *vnewr, double *enl, double *enr);
  void try_deleting(Region*,double *vnewl, double *vnewr, double *enl, double *enr);
  void reg_force(int, Region*,double*,int);
  void reg_force_perp(int, Region*,double*,int);
  double g_par_global(Region*, int);
  double g_perp_global(Region*, int);
  double g_par_local(Region*, double, int);
  double energy(int, int, double *, double *);
  double usher(double **, double, int, int,int&);
  double near_energy(double **,int,int);
  void center_of_mass(int, double **, double *);
  void calc_torque(int, double**, double*, double*, double*);
  int check_proc(double**, double*, double *, double*, int);
  int check_mol_proc(double**, double*, double*, double*, double*, int, int,int&);
  int check_mol_region(Region*, double **, int);
  void vcm_internal_sq(int, double, double *, double *, double *, Region *);
  //double molecule_energy(tagint);
  double xvalue, yvalue, zvalue;
  int varflag;
  //char *xstr, *ystr, *zstr, *estr;
  //int xvar, yvar, zvar, evar, xstyle, ystyle, zstyle, estyle;
  double foriginal[4], foriginal_all[4];
  int force_flag;
  //int ilevel_respa;
  int print_delete_flag;
  int print_insert_flag;
  int stepflag;
  double g_fac, g_fac_inv;
  int maxatom;
  double **sforce;

  char *id_press;
  class Compute *pressure_tensor;

  double mtot;

  double fusher[3];
  double shearForce_left[3];
  double shearForce_right[3];
  double **cutsq;
  class Pair *pair;
  int nattempt;
  double uovlp, dsovlp, eps, ds0, dtheta0, etarget;

};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Invalid atom type in fix deposit command

Self-explanatory.

E: Must specify a region in fix deposit

The region keyword must be specified with this fix.

E: Fix deposit region does not support a bounding box

Not all regions represent bounded volumes.  You cannot use
such a region with the fix deposit command.

E: Fix deposit region cannot be dynamic

Only static regions can be used with fix deposit.

E: Deposition region extends outside simulation box

Self-explanatory.

E: Cannot use fix_deposit unless atoms have IDs

Self-explanatory.

E: Fix deposit molecule must have coordinates

The defined molecule does not specify coordinates.

E: Fix deposit molecule must have atom types

The defined molecule does not specify atom types.

E: Invalid atom type in fix deposit mol command

The atom types in the defined molecule are added to the value
specified in the create_atoms command, as an offset.  The final value
for each atom must be between 1 to N, where N is the number of atom
types.

E: Fix deposit molecule template ID must be same as atom_style template ID

When using atom_style template, you cannot deposit molecules that are
not in that template.

E: Cannot use fix deposit rigid and not molecule

Self-explanatory.

E: Cannot use fix deposit shake and not molecule

Self-explanatory.

E: Cannot use fix deposit rigid and shake

These two attributes are conflicting.

E: Region ID for fix deposit does not exist

Self-explanatory.

E: Fix deposit rigid fix does not exist

UNDOCUMENTED

E: Fix deposit and fix rigid/small not using same molecule template ID

Self-explanatory.

E: Fix deposit shake fix does not exist

Self-explanatory.

E: Fix deposit and fix shake not using same molecule template ID

Self-explanatory.

W: Fix deposit near setting < possible overlap separation %g

This test is performed for finite size particles with a diameter, not
for point particles.  The near setting is smaller than the particle
diameter which can lead to overlaps.

E: Unknown particle distribution in fix deposit

UNDOCUMENTED

W: Particle deposition was unsuccessful

The fix deposit command was not able to insert as many atoms as
needed.  The requested volume fraction may be too high, or other atoms
may be in the insertion region.

E: Too many total atoms

See the setting for bigint in the src/lmptype.h file.

E: New atom IDs exceed maximum allowed ID

See the setting for tagint in the src/lmptype.h file.

E: Molecule template ID for fix deposit does not exist

Self-explanatory.

U: Fix pour rigid fix does not exist

Self-explanatory.

*/
