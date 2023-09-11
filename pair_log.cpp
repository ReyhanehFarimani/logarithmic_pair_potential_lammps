// clang-format off

/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "pair_log.h"

#include <cmath>
#include <cstring>
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "neigh_list.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairLog::PairLog(LAMMPS *lmp) : Pair(lmp) {
    
  writedata = 1;
}

/*---------------------------------------------------------------------- */
PairLog::~PairLog()
{
    if(allocated)
    {
    memory->destroy(setflag);


    memory->destroy(cut);
    memory->destroy(cutsq);
    memory->destroy(const_multi);
    memory->destroy(R);
    memory->destroy(Rsq);
    }
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */
void PairLog::allocate()
{
    
    allocated = 1;
    int n = atom->ntypes;
    int i, j; //iteration dummy variable
    
    memory->create(setflag, n+1, n+1, "pair:setflag");
    for (i = 0; i < n+1 ; i++)
    {
        for (j = 0; j< n+1; j++)
            setflag[i][j] = 0;
    }

    memory->create(cut, n+1, n+1, "pair:cut");
    memory->create(cutsq, n+1, n+1, "pair:cutsq");

    memory->create(const_multi, n+1, n+1, "pair:const_multi");
    memory->create(R, n+1, n+1, "pair:R");
    memory->create(Rsq, n+1, n+1, "pair:Rsq");
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */
void PairLog::settings(int narg, char **arg)
{
    //check for the number of parsed arguements
    if (narg != 1)
        error->all(FLERR, "Illegal pair_style command");

    cut_global = utils::numeric(FLERR,arg[0],false,lmp);

    // reset cutoffs that have been explicitly set
    if(allocated)
    {
        int i,j;
        for(i = 1 ; i<= atom->ntypes; i++)
            for(j = i; j<= atom->ntypes; j++)
                if(setflag[i][j])
                {
                    cut[i][j] = cut_global;
                }
    }
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairLog::coeff(int narg, char **arg)
{
    
    if (narg>5) 
        error->all(FLERR, "Incorrect args for pair coefficients");

    if (!allocated) 
        allocate();
    
    int ilo, ihi, jlo, jhi;
    utils::bounds(FLERR,arg[0],1,atom->ntypes,ilo,ihi,error);
    utils::bounds(FLERR,arg[1],1,atom->ntypes,jlo,jhi,error);
    
    double cut_one = cut_global;
    double const_multi_one;
    double R_one;
    
    cut_one = utils::numeric(FLERR,arg[2],false,lmp);
    const_multi_one = utils::numeric(FLERR,arg[3],false,lmp);
    R_one = utils::numeric(FLERR,arg[4],false,lmp);
    
    if (R_one <= 0.0 || R_one > cut_one)
        error->all(FLERR,"Illegal pair_style command");

    int count =0;
    for (int i = ilo; i <= ihi; i++)
    {
        
        for (int j = MAX(jlo, i); j <= jhi; j++)
        {
            R[i][j] = R_one;
            const_multi[i][j] = const_multi_one;
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

double PairLog::init_one(int i, int j)
{
    
    if(setflag[i][j] == 0) error->all(FLERR,"Error! mix energy!");
    cut[j][i] = cut[i][j];
    cutsq[i][j] = cut[i][j]*cut[i][j];
    cutsq[j][i] = cutsq[i][j];
    R[j][i] = R[i][j];
    Rsq[i][j] = R[i][j]*R[i][j];
    Rsq[j][i] = Rsq[i][j];
    const_multi[j][i] = const_multi[i][j];

    return cut[i][j];
}

/* ---------------------------------------------------------------------- */

void PairLog::compute(int eflag, int vflag)
{
    
    int i, j, ii, jj, inum, jnum, itype, jtype;
    double xtmp, ytmp, ztmp, delx, dely, delz;
    double evdwl, fpair;
    double rsq, r2inv, forcelog, factor_lj;
    int *ilist, *jlist, *numneigh, **firstneigh;    

    evdwl = 0.0;
    ev_init(eflag, vflag);

    double **x = atom->x;
    double **f = atom->f;
    int *type = atom->type;
    int nlocal = atom->nlocal;
    double *special_lj = force->special_lj;
    int newton_pair = force->newton_pair;


    inum = list->inum;
    ilist = list->ilist;
    numneigh = list->numneigh;
    firstneigh = list->firstneigh;

    // loop over neighbors of my atoms

    for (ii = 0; ii < inum; ii++)
    {
        i = ilist[ii];
        xtmp = x[i][0];
        ytmp = x[i][1];
        ztmp = x[i][2];

        itype = type[i];
        jlist = firstneigh[i];
        jnum = numneigh[i];
        
        for(jj = 0; jj< jnum; jj++)
        {
            j = jlist[jj];
            factor_lj = special_lj[sbmask(j)];
            j &= NEIGHMASK;

            delx = xtmp - x[j][0];
            dely = ytmp - x[j][1];
            delz = ztmp - x[j][2];

            jtype = type[j];
               
            rsq = delx * delx + dely * dely + delz * delz;
            if (rsq <= Rsq[itype][jtype])
            {
                delx = delx/R[itype][jtype];
                dely = dely/R[itype][jtype];
                delz = delz/R[itype][jtype];   
                rsq = delx * delx + dely * dely + delz * delz;  
                r2inv = 1.0/rsq;
                forcelog = -const_multi[itype][jtype] * r2inv;
                fpair = factor_lj* forcelog;

                f[i][0] += delx*fpair;
                f[i][1] += dely*fpair;
                f[i][2] += delz*fpair;

                if (newton_pair || j < nlocal)
                {
                    f[j][0] -= delx*fpair;
                    f[j][1] -= dely*fpair;
                    f[j][2] -= delz*fpair;
                }

                if (eflag) 
                {
                    evdwl = const_multi[itype][jtype] * ( -log(rsq)/2 + 1/2.71828/2)*2;
                    evdwl *= factor_lj;
                }
                if (evflag) ev_tally(i, j, nlocal, newton_pair, evdwl, 0.0, fpair, delx, dely, delz);

            }

            else if (rsq <= cutsq[itype][jtype])
            {
                delx = delx/R[itype][jtype];
                dely = dely/R[itype][jtype];
                delz = delz/R[itype][jtype];   
                rsq = delx * delx + dely * dely + delz * delz;  
                r2inv = 1.0/rsq;
                forcelog = -const_multi[itype][jtype] * pow(2.71828, (-rsq));
                fpair = factor_lj* forcelog;

                f[i][0] += delx*fpair;
                f[i][1] += dely*fpair;
                f[i][2] += delz*fpair;

                if (newton_pair || j < nlocal)
                {
                    f[j][0] -= delx*fpair;
                    f[j][1] -= dely*fpair;
                    f[j][2] -= delz*fpair;
                }

                if (eflag) 
                {
                    evdwl = const_multi[itype][jtype] * (pow(2.71828, (-rsq)));
                    evdwl *= factor_lj;
                }
                if (evflag) ev_tally(i, j, nlocal, newton_pair, evdwl, 0.0, fpair, delx, dely, delz);

            }
        }
        
    }
    if (vflag_fdotr) virial_fdotr_compute();
}



/* ----------------------------------------------------------------------
  proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairLog::write_restart(FILE *fp)
{
    write_restart_settings(fp);

    int i, j;
    for(i = 1 ; i < atom->ntypes; i++)
    {
        for (j = 1; j< atom->ntypes; j++)
        {
            fwrite(&setflag[i][j], sizeof(int), 1, fp);
            if (setflag[i][j])
            {
                fwrite(&R[i][j], sizeof(double), 1, fp);
                fwrite(&const_multi[i][j], sizeof(double), 1, fp);
                fwrite(&cut[i][j], sizeof(double), 1, fp);
            }
        }
    }
}

/* ----------------------------------------------------------------------
  proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairLog::read_restart(FILE *fp)
{
    read_restart_settings(fp);
    allocate();

    int i, j;
    int me = comm->me;

    for (i = 0 ; i<= atom->ntypes; i++)
    {
        for (j = 0 ; j<= atom->ntypes; j++)
        {
            if (me == 0)
                utils::sfread(FLERR, &setflag[i][j], sizeof(int), 1, fp, nullptr, error);
            MPI_Bcast(&setflag[i][j] ,1 , MPI_INT, 0 , world);
            if (setflag[i][j])
            {
                if (me ==0)
                {
                    utils::sfread(FLERR,&R[i][j],sizeof(double),1,fp,nullptr,error);
                    utils::sfread(FLERR,&const_multi[i][j],sizeof(double),1,fp,nullptr,error);
                    utils::sfread(FLERR,&cut[i][j],sizeof(double),1,fp,nullptr,error);
                }

                MPI_Bcast(&R[i][j],1,MPI_DOUBLE,0,world);
                MPI_Bcast(&const_multi[i][j],1,MPI_DOUBLE,0,world);
                MPI_Bcast(&cut[i][j],1,MPI_DOUBLE,0,world);
            }
        } 
    }
}


/* ----------------------------------------------------------------------
  proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairLog::write_restart_settings(FILE *fp)
{
  fwrite(&mix_flag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
  proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairLog::read_restart_settings(FILE *fp)
{
    int me = comm->me;
    if (me == 0)
    {
        utils::sfread(FLERR,&mix_flag,sizeof(int),1,fp,nullptr,error);
    }
    MPI_Bcast(&mix_flag,1,MPI_INT,0,world);
}
/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void PairLog::write_data(FILE *fp)
{

    for (int i = 1; i <= atom->ntypes; i++)
        fprintf(fp,"%d %g %g\n",i,const_multi[i][i],R[i][i]);
}
/* ----------------------------------------------------------------------
   proc 0 writes all pairs to data file
------------------------------------------------------------------------- */

void PairLog::write_data_all(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    for (int j = i; j <= atom->ntypes; j++)
      fprintf(fp,"%d %d %g %g %g\n",
              i,j,const_multi[i][j],R[i][j],cut[i][j]);
}
/* ---------------------------------------------------------------------- */

void *PairLog::extract(const char *str, int &dim)
{
  dim = 2;
  if (strcmp(str,"R") == 0) return (void *) R;
  if (strcmp(str,"const_multi") == 0) return (void *) const_multi;
  return nullptr;
}
