/* BayesX - Software for Bayesian Inference in
Structured Additive Regression Models.
Copyright (C) 2011  Christiane Belitz, Andreas Brezger,
Thomas Kneib, Stefan Lang, Nikolaus Umlauf

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA. */



#if !defined (FCNONPVARIANCEINCLUDED)

#define FCNONPVARIANCEINCLUDED

#include"../export_type.h"
#include"../values.h"
#include<fstream>
#include"GENERAL_OPTIONS.h"
#include"clstring.h"
#include"FC_nonp.h"
#include"design.h"
#include<cmath>

namespace MCMC
{

//------------------------------------------------------------------------------
//--------------------------- CLASS: FC_nonp_variance --------------------------
//------------------------------------------------------------------------------


class __EXPORT_TYPE FC_nonp_variance  : public FC
  {

  protected:

  FC_nonp * FCnonpp;                         // Pointer to corresponding
                                             // FC_nonp object
  DISTR * likep;                             // Pointer to DISTR obejct
  DESIGN * designp;                          // Pointer to design object

  MASTER_OBJ * masterp;
  unsigned equationnr;

  double a_invgamma;
  double b_invgamma_orig;
  double b_invgamma;
  double lambdastart;

  bool lambdaconst;

  double tildea;
  double tildeb;
  bool cauchy;
  bool wei;
  double scaletau2;

  public:

//----------------------- CONSTRUCTORS, DESTRUCTOR -----------------------------

  // DEFAULT CONSTRUCTOR

  FC_nonp_variance(void);

  // CONSTRUCTOR
  // o    : pointer to GENERAL_OPTIONS object
  // t    : title of the full conditional (for example "fixed effects")
  // fp   : file path for storing sampled parameters

  FC_nonp_variance(MASTER_OBJ * mp,unsigned & enr, GENERAL_OPTIONS * o,DISTR * lp,
           const ST::string & t,
           const ST::string & fp,DESIGN * dp,FC_nonp * FCn,
           vector<ST::string> & op,vector<ST::string> & vn);

  // COPY CONSTRUCTOR

  FC_nonp_variance(const FC_nonp_variance & m);


  // OVERLOADED ASSIGNMENT OPERATOR

  const FC_nonp_variance & operator=(const FC_nonp_variance & m);

  // DESTRUCTOR

  ~FC_nonp_variance()
    {
    }

  // FUNCTION: update
  // TASK: - stores sampled parameters in file 'samplepath'
  //         storing order: first row, second row, ...

  void update(void);

  // FUNCTION: posteriormode
  // TASK: computes the posterior mode

  bool posteriormode(void);

  // FUNCTION: outoptions
  // TASK: writes estimation options (hyperparameters, etc.) to outputstream

  void outoptions(void);

  // FUNCTION: outresults
  // TASK: writes estimation results to logout or into a file

  void outresults(ofstream & out_stata,ofstream & out_R,
                  const ST::string & pathresults);

  // FUNCTION: reset
  // TASK: resets all parameters

  void reset(void);

  void read_options(vector<ST::string> & op,vector<ST::string> & vn);

  // virtual void transform_beta(void);

  };





class __EXPORT_TYPE FC_nonp_variance_varselection  : public FC_nonp_variance
  {

  protected:

  bool singleomega;

  FC FC_psi2;
  FC FC_omega;

  double a_omega;
  double b_omega;

  double v;
  double Q;

  double r;

  datamatrix X;

  double tauold;

  public:

  double omegas;
  FC FC_delta;

//----------------------- CONSTRUCTORS, DESTRUCTOR -----------------------------

  // DEFAULT CONSTRUCTOR

  FC_nonp_variance_varselection(void);

  // CONSTRUCTOR
  // o    : pointer to GENERAL_OPTIONS object
  // t    : title of the full conditional (for example "fixed effects")
  // fp   : file path for storing sampled parameters

  FC_nonp_variance_varselection(MASTER_OBJ * mp,unsigned & enr, GENERAL_OPTIONS * o,DISTR * lp,
          bool so,
           const ST::string & t,
           const ST::string & fp,DESIGN * dp,FC_nonp * FCn,
           vector<ST::string> & op,vector<ST::string> & vn);

  // COPY CONSTRUCTOR

  FC_nonp_variance_varselection(const FC_nonp_variance_varselection & m);


  // OVERLOADED ASSIGNMENT OPERATOR

  const FC_nonp_variance_varselection & operator=(const FC_nonp_variance_varselection & m);

  // DESTRUCTOR

  ~FC_nonp_variance_varselection()
    {
    }

  // FUNCTION: update
  // TASK: - stores sampled parameters in file 'samplepath'
  //         storing order: first row, second row, ...

  void update(void);

  bool posteriormode(void);

  // FUNCTION: outoptions
  // TASK: writes estimation options (hyperparameters, etc.) to outputstream

  void outoptions(void);

  // FUNCTION: outresults
  // TASK: writes estimation results to logout or into a file

  void outresults(ofstream & out_stata,ofstream & out_R,
                  const ST::string & pathresults);

  // FUNCTION: reset
  // TASK: resets all parameters

  void reset(void);

  void read_options(vector<ST::string> & op,vector<ST::string> & vn);

  void get_samples(const ST::string & filename,ofstream & outg) const;


  // virtual void transform_beta(void);

  };



class __EXPORT_TYPE FC_nonp_variance_varselection2  : public FC_nonp_variance
  {

  protected:

  bool singleomega;


  FC FC_psi2;
  FC FC_omega;

  double a_omega;
  double b_omega;

  double v;
  double Q;

//  double scaletau2;

  double r2;

  datamatrix X;

  public:

  double omegas;

  FC FC_delta;

//----------------------- CONSTRUCTORS, DESTRUCTOR -----------------------------

  // DEFAULT CONSTRUCTOR

  FC_nonp_variance_varselection2(void);

  // CONSTRUCTOR
  // o    : pointer to GENERAL_OPTIONS object
  // t    : title of the full conditional (for example "fixed effects")
  // fp   : file path for storing sampled parameters

  FC_nonp_variance_varselection2(MASTER_OBJ * mp,unsigned & enr, GENERAL_OPTIONS * o,DISTR * lp,
           bool so,
            const ST::string & t,
           const ST::string & fp,DESIGN * dp,FC_nonp * FCn,
           vector<ST::string> & op,vector<ST::string> & vn);

  // COPY CONSTRUCTOR

  FC_nonp_variance_varselection2(const FC_nonp_variance_varselection2 & m);


  // OVERLOADED ASSIGNMENT OPERATOR

  const FC_nonp_variance_varselection2 & operator=(const FC_nonp_variance_varselection2 & m);

  // DESTRUCTOR

  ~FC_nonp_variance_varselection2()
    {
    }

  // FUNCTION: update
  // TASK: - stores sampled parameters in file 'samplepath'
  //         storing order: first row, second row, ...

  void update(void);

  bool posteriormode(void);

  // FUNCTION: outoptions
  // TASK: writes estimation options (hyperparameters, etc.) to outputstream

  void outoptions(void);

  // FUNCTION: outresults
  // TASK: writes estimation results to logout or into a file

  void outresults(ofstream & out_stata,ofstream & out_R,
                  const ST::string & pathresults);

  // FUNCTION: reset
  // TASK: resets all parameters

  void reset(void);

  void read_options(vector<ST::string> & op,vector<ST::string> & vn);

  void get_samples(const ST::string & filename,ofstream & outg) const;


  // virtual void transform_beta(void);

  };



class __EXPORT_TYPE FC_varselection_omega  : public FC
  {

  protected:

  double a_omega;
  double b_omega;

  DISTR * likep;                             // Pointer to DISTR obejct

  MASTER_OBJ * masterp;


  public:

  vector<FC_nonp_variance_varselection*> FC_tau2s;
 // vector<FC_nonp_variance_varselection*> FC_tau2s;


//----------------------- CONSTRUCTORS, DESTRUCTOR -----------------------------

  // DEFAULT CONSTRUCTOR

  FC_varselection_omega(void);

  // CONSTRUCTOR
  // o    : pointer to GENERAL_OPTIONS object
  // t    : title of the full conditional (for example "fixed effects")
  // fp   : file path for storing sampled parameters

  FC_varselection_omega(MASTER_OBJ * mp,unsigned & enr, GENERAL_OPTIONS * o,DISTR * lp,
           const ST::string & t);

  // COPY CONSTRUCTOR

  FC_varselection_omega(const FC_varselection_omega & m);


  // OVERLOADED ASSIGNMENT OPERATOR

  const FC_varselection_omega & operator=(const FC_varselection_omega & m);

  // DESTRUCTOR

  ~FC_varselection_omega()
    {
    }

  // FUNCTION: update
  // TASK: - stores sampled parameters in file 'samplepath'
  //         storing order: first row, second row, ...

  void update(void);

//  bool posteriormode(void);
//
//  // FUNCTION: outoptions
//  // TASK: writes estimation options (hyperparameters, etc.) to outputstream
//
//  void outoptions(void);
//
//  // FUNCTION: outresults
//  // TASK: writes estimation results to logout or into a file
//
//  void outresults(ofstream & out_stata,ofstream & out_R,
//                  const ST::string & pathresults);
//
//  // FUNCTION: reset
//  // TASK: resets all parameters
//
//  void reset(void);
//
//  void read_options(vector<ST::string> & op,vector<ST::string> & vn);
//
//  void get_samples(const ST::string & filename,ofstream & outg) const;
//
//
  // virtual void transform_beta(void);

  };




} // end: namespace MCMC

#endif


