/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS
// clang-format off
PairStyle(log,PairLog);
// clang-format on
#else

#ifndef LMP_PAIR_LOG_H
#define LMP_PAIR_LOG_H

#include "pair.h"


namespace LAMMPS_NS
{

class PairLog : public Pair
{
    public:
        PairLog(class LAMMPS *);
        ~PairLog() override;

        void compute(int, int) override;
        void settings(int, char**) override;
        void coeff(int, char**) override;
        double init_one(int, int) override;
        void write_restart(FILE *) override;
        void read_restart(FILE *) override;
        void write_restart_settings(FILE *) override;
        void read_restart_settings(FILE *) override;
        void write_data (FILE *) override;
        void write_data_all(FILE *) override;

        void *extract(const char *, int &) override;

    protected:
        double cut_global, const_multi_global, R_global;
        double **cut;
        double **const_multi, **R; 

        void allocate();

};

}//namespace LAMMPS_NS

#endif
#endif