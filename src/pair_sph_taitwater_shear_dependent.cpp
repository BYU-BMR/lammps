/* ----------------------------------------------------------------------
 LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
 http://lammps.sandia.gov, Sandia National Laboratories
 Steve Plimpton, sjplimp@sandia.gov

 Copyright (2003) Sandia Corporation.  Under the terms of Contract
 DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
 certain rights in this software.  This software is distributed under
 the GNU General Public License.

 See the README file in the top-level LAMMPS directory.
 ------------------------------------------------------------------------- */

#include <math.h>
#include <stdlib.h>
#include "pair_sph_taitwater_shear_dependent.h"
#include "atom.h"
#include "force.h"
#include "comm.h"
#include "neigh_list.h"
#include "memory.h"
#include "error.h"
#include "domain.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairSPHTaitwaterShearDependent::PairSPHTaitwaterShearDependent(LAMMPS *lmp) : Pair(lmp)
{
  restartinfo = 0;
  first = 1;
}

/* ---------------------------------------------------------------------- */

PairSPHTaitwaterShearDependent::~PairSPHTaitwaterShearDependent() {
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);

    memory->destroy(cut);
    memory->destroy(rho0);
    memory->destroy(soundspeed);
    memory->destroy(B);
    memory->destroy(viscosity);
    memory->destroy(mu0);
    memory->destroy(mup1);
    memory->destroy(M);
  }
}

/* ---------------------------------------------------------------------- */

void PairSPHTaitwaterShearDependent::compute(int eflag, int vflag) {
  int i, j, ii, jj, inum, jnum, itype, jtype;
  double xtmp, ytmp, ztmp, delx, dely, delz, fpair, flubx, fluby, flubz, sx, sy, sz, r,vsq, discritize_shear_rate;

  int *ilist, *jlist, *numneigh, **firstneigh;
  double vxtmp, vytmp, vztmp, imass, jmass, fi, fj, fvisc, h, ih, ihsq, velx, vely, velz, localshearRate, viscosity_ij;
  double rsq, tmp, wfd, delVdotDelR, deltaE;
  if (eflag || vflag)
    ev_setup(eflag, vflag);
  else
    evflag = vflag_fdotr = 0;
  
  double **v = atom->vest;
  double **x = atom->x;
  double **f = atom->f;
  double *rho = atom->rho;
  double *mass = atom->mass;
  double *de = atom->de;
  double *drho = atom->drho;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int newton_pair = force->newton_pair;
  
 
  // check consistency of pair coefficients

  if (first) {
    for (i = 1; i <= atom->ntypes; i++) {
      for (j = 1; i <= atom->ntypes; i++) {
        if (cutsq[i][j] > 1.e-32) {
          if (!setflag[i][i] || !setflag[j][j]) {
            if (comm->me == 0) {
              printf(
                  "SPH particle types %d and %d interact with cutoff=%g, but not all of their single particle properties are set.\n",
                  i, j, sqrt(cutsq[i][j]));
            }
          }
        }
      }
    }
    first = 0;
  }

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  //loop for shear rate

   for (ii = 0; ii < inum; ii++) {
     i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    vxtmp = v[i][0];
    vytmp = v[i][1];
    vztmp = v[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    imass = mass[itype];

    localshearRate=1;  // it's better to change this part.

     for (jj = 0; jj < jnum; jj++) {
        j = jlist[jj];
        j &= NEIGHMASK;

        delx = xtmp - x[j][0];
        dely = ytmp - x[j][1];
        delz = ztmp - x[j][2];
        rsq = delx * delx + dely * dely + delz * delz;
        jtype = type[j];
        jmass = mass[jtype];

        velx=vxtmp - v[j][0];
        vely=vytmp - v[j][1];
        velz=vztmp - v[j][2];
        vsq=velx*velx+vely*vely+velz*velz;

        if (rsq < cutsq[itype][jtype]) {
            h = cut[itype][jtype]; // radius of the particles
            ih = 1.0 / h;
            ihsq = ih * ih;
            r=sqrt(rsq);
            wfd = h - r;

            wfd = -25.066903536973515383e0 * wfd * wfd * ihsq * ihsq * ihsq * ih;
            
            //printf("wfd: %f r: %f vsq: %f jmass: %f\n", wfd, r, vsq, jmass);
            discritize_shear_rate= (wfd*r*(sqrt(vsq))*jmass/rho[j]);
            localshearRate += discritize_shear_rate*discritize_shear_rate;
          }
        }
    localshearRate=sqrt(localshearRate)/2;

    //printf("itype: %d jtype: %d mup1[jtype]: %f mu0[jtype]: %f localshearRate: %f M[jtype]: %f\n", itype, jtype, mup1[jtype], mu0[jtype], localshearRate, M[jtype]);
        
    // Morris Viscosity (Morris, 1996);
    viscosity[ii]=mup1[jtype]+ mu0[jtype]/localshearRate*M[jtype]*(1-exp(-M[jtype]*localshearRate));
    //printf("saved_viscosity: %f\n", viscosity[ii]);
          
  }
  // loop over neighbors of my atoms;

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    vxtmp = v[i][0];
    vytmp = v[i][1];
    vztmp = v[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    imass = mass[itype];

    // compute pressure of atom i with Tait EOS
    tmp = rho[i] / rho0[itype];
    fi = tmp * tmp * tmp;
    fi = B[itype] * (fi * fi * tmp - 1.0) / (rho[i] * rho[i]);

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx * delx + dely * dely + delz * delz;
      jtype = type[j];
      jmass = mass[jtype];

      if (rsq < cutsq[itype][jtype]) {
        h = cut[itype][jtype]; // radius of the particles
        ih = 1.0 / h;
        ihsq = ih * ih;

        wfd = h - sqrt(rsq);
        if (domain->dimension == 3) {
          // Lucy Kernel, 3d
          // Note that wfd, the derivative of the weight function with respect to r,
          // is lacking a factor of r.
          // The missing factor of r is recovered by
          // (1) using delV . delX instead of delV . (delX/r) and
          // (2) using f[i][0] += delx * fpair instead of f[i][0] += (delx/r) * fpair
          wfd = -25.066903536973515383e0 * wfd * wfd * ihsq * ihsq * ihsq * ih;
        } else {
          // Lucy Kernel, 2d
          wfd = -19.098593171027440292e0 * wfd * wfd * ihsq * ihsq * ihsq;
        }  


        
        //printf("viscosity[ii]: %f viscosity[jj]: %f\n", viscosity[ii],viscosity[jj]);
        //printf("viscosity[i]: %f viscosity[j]: %f\n", viscosity[i],viscosity[j]);
        viscosity_ij=(viscosity[ii]+viscosity[jj])/2;
        // compute pressure  of atom j with Tait EOS
        tmp = rho[j] / rho0[jtype];
        fj = tmp * tmp * tmp;
        fj = B[jtype] * (fj * fj * tmp - 1.0) / (rho[j] * rho[j]);

        velx=vxtmp - v[j][0];
        vely=vytmp - v[j][1];
        velz=vztmp - v[j][2];
        // define flub
        sx=abs(delx)-2*h;
        sy=abs(dely)-2*h;
        sz=abs(delz)-2*h;
        flubx = (3/2)*3.14*viscosity_ij*h*h*(velx)/(sx) + (soundspeed[jtype]/ pow(sx, 8));
        fluby = (3/2)*3.14*viscosity_ij*h*h*(vely)/(sy) + (soundspeed[jtype]/ pow(sy, 8));
        flubz = (3/2)*3.14*viscosity_ij*h*h*(velz)/(sz) + (soundspeed[jtype]/ pow(sz, 8));
        //printf("flubx: %f fluby: %f flubz: %f\n", flubx, fluby, flubz);
        // dot product of velocity delta and distance vector
        delVdotDelR = delx * velx + dely * vely + delz * velz;

        
        fvisc = 2 * viscosity_ij / (rho[i] * rho[j]);

        fvisc *= imass * jmass * wfd;
        //printf("fvisc: %f\n", fvisc);
        // total pair force & thermal energy increment
        fpair = -imass * jmass * (fi + fj) * wfd;
        deltaE = -0.5 *(fpair * delVdotDelR + fvisc * (velx*velx + vely*vely + velz*velz));

       // printf("testvar= %f, %f \n", delx, dely);
// add all the forces for i particle
        f[i][0] += delx * fpair + velx * fvisc + flubx;
        f[i][1] += dely * fpair + vely * fvisc + fluby;
        f[i][2] += delz * fpair + velz * fvisc + flubz;

        // and change in density
        drho[i] += jmass * delVdotDelR * wfd;

        // change in thermal energy
        de[i] += deltaE;
// apply the third law of motion
        if (newton_pair || j < nlocal) {
          f[j][0] -= delx * fpair + velx * fvisc + flubx;
          f[j][1] -= dely * fpair + vely * fvisc + fluby;
          f[j][2] -= delz * fpair + velz * fvisc + flubz;
          de[j] += deltaE;
          drho[j] += imass * delVdotDelR * wfd;
        }

        if (evflag)
          ev_tally(i, j, nlocal, newton_pair, 0.0, 0.0, fpair, delx, dely, delz);
      }
    }
  }

  if (vflag_fdotr) virial_fdotr_compute();
}


/* ----------------------------------------------------------------------
 allocate all arrays
 ------------------------------------------------------------------------- */

void PairSPHTaitwaterShearDependent::allocate() {
  allocated = 1;
  int n = atom->ntypes;
  int m = atom->natoms;  // sublime???


  memory->create(setflag, n + 1, n + 1, "pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq, n + 1, n + 1, "pair:cutsq");

  memory->create(rho0, n + 1, "pair:rho0");
  memory->create(soundspeed, n + 1, "pair:soundspeed");
  memory->create(B, n + 1, "pair:B");
  memory->create(cut, n + 1, n + 1, "pair:cut");
  memory->create(viscosity, m + 1, "pair:viscosity");  // we wanr an array of viscosity
  memory->create(mu0, n+1, "pair: mu0");
  memory->create(mup1, n+1, "pair: mup1");
  memory->create(M, n+1, "pair: M");
}

/* ----------------------------------------------------------------------
 global settings
 ------------------------------------------------------------------------- */

void PairSPHTaitwaterShearDependent::settings(int narg, char **arg) {
  if (narg != 0)
    error->all(FLERR,
        "Illegal number of setting arguments for pair_style sph/taitwater/shear/dependent");
}

/* ----------------------------------------------------------------------
 set coeffs for one or more type pairs
 ------------------------------------------------------------------------- */

void PairSPHTaitwaterShearDependent::coeff(int narg, char **arg) {        // coeff reads
  if (narg != 8)
    error->all(FLERR,
        "Incorrect args for pair_style sph/taitwater/shear/dependent coefficients");
  if (!allocated)
    allocate();

  int ilo, ihi, jlo, jhi;
  force->bounds(FLERR,arg[0], atom->ntypes, ilo, ihi);
  force->bounds(FLERR,arg[1], atom->ntypes, jlo, jhi);
  double rho0_one = force->numeric(FLERR,arg[2]);
  double soundspeed_one = force->numeric(FLERR,arg[3]);
  double cut_one = force->numeric(FLERR,arg[4]);
  double B_one = soundspeed_one * soundspeed_one * rho0_one / 7.0;
  double mu0_one = force->numeric(FLERR, arg[5]);
  double mup1_one = force->numeric(FLERR, arg[6]);
  double M_one = force->numeric(FLERR, arg[7]);


  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    rho0[i] = rho0_one;
    soundspeed[i] = soundspeed_one;
    B[i] = B_one;
    mu0[i] = mu0_one;
    mup1[i] = mup1_one;
    M[i] = M_one;
    // printf("i: %d mu0[i]: %f mup1[i]: %f M[i]: %f\n", i, mu0[i], mup1[i], M[i]);

    for (int j = MAX(jlo,i); j <= jhi; j++) {
      
      //printf("setting cut[%d][%d] = %f\n", i, j, cut_one);
      cut[i][j] = cut_one;

      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0)
    error->all(FLERR,"Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
 init for one type pair i,j and corresponding j,i
 ------------------------------------------------------------------------- */

double PairSPHTaitwaterShearDependent::init_one(int i, int j) {

  if (setflag[i][j] == 0) {
    error->all(FLERR,"Not all pair sph/taitwater/shear/dependent coeffs are not set");
  }

  cut[j][i] = cut[i][j];
 
  return cut[i][j];
}

/* ---------------------------------------------------------------------- */

double PairSPHTaitwaterShearDependent::single(int i, int j, int itype, int jtype,
    double rsq, double factor_coul, double factor_lj, double &fforce) {
  fforce = 0.0;

  return 0.0;
}
