
#include "randomeffect.h"
#include "randomeffect_stepwise.h"

namespace MCMC
{


//------------------------------------------------------------------------------
//----------------- class FULLCOND_random --------------------------------------
//------------------------------------------------------------------------------

void FULLCOND_random::compute_lambdavec(vector<double> & lvec,unsigned & number)
  {
  if (get_df_equidist()==true)
     FULLCOND::compute_lambdavec_equi(lvec,number);
  else
     FULLCOND::compute_lambdavec(lvec,number);
  if (randomslope)
    lvec.push_back(-1);
  get_forced();
  if(forced_into==false)
     lvec.push_back(0);
  }


const datamatrix & FULLCOND_random::get_data_forfixedeffects(void)
  {
  // useful for randomslopes only
  if ( (data_forfixed.rows() < data.rows()) && (randomslope==true) )
    {
    data_forfixed=datamatrix(data.rows(),1);
    unsigned i;
    int * workindex = index.getV();
    double * workdata = data.getV();
    for (i=0;i<data.rows();i++,workindex++,workdata++)
      {
      data_forfixed(*workindex,0) = *workdata;
      }
    }

  return data_forfixed;
  }
  

double FULLCOND_random::compute_df(void)
  {
  unsigned i;
  double df=0;
  double * workXX=XX.getV();
  unsigned n;

  if (randomslope && includefixed)
    {
    n = nrpar-1;
    df=1;
    }
  else
    n = nrpar;


  if ((lambdaold1==lambda) && (likep->iwlsweights_constant() == true) )
    {
    df = df_lambdaold1;
    }
  else if ((lambdaold2==lambda) && (likep->iwlsweights_constant() == true) )
    {
    df = df_lambdaold2;
    }
  else
    {
    for(i=0;i<n;i++,workXX++)
      {
      df += (*workXX)/(*workXX+lambda);
      }

    df_lambdaold2 = df_lambdaold1;
    lambdaold2 = lambdaold1;
    df_lambdaold1 = df;
    lambdaold1 = lambda;

    }

  return df;
  }


void FULLCOND_random::update_stepwise(double la)
  {
  lambda = la;
  }


ST::string FULLCOND_random::get_effect(void)
  {

  ST::string h;

  if(randomslope)
    h = datanames[1] + "*" + datanames[0];
  else
    h = datanames[0];

  h = h + "(random,df=" + ST::doubletostring(compute_df(),6) + ",(lambda=" + ST::doubletostring(lambda,6) + "))";

  return h;

  }


void FULLCOND_random::reset_effect(unsigned & pos)
  {

  update_linpred(false);

  unsigned i;
  double * work;
  work = beta.getV();
  for(i=0;i<nrpar;i++,work++)
    *work = 0.0;

  }

/*
//------------------------------------------------------------------------------
//------------------- class FULLCOND_random_gaussian ---------------------------
//------------------------------------------------------------------------------

double FULLCOND_random_gaussian::compute_df(void)
  {
  unsigned i;
  double df=0;
  double * workXX=XX.getV();
  unsigned n;

  if (randomslope && includefixed)
    {
    n = nrpar-1;
    df=1;
    }
  else
    n = nrpar;


  if ((lambdaold1==lambda) && (likep->iwlsweights_constant() == true) )
    {
    df = df_lambdaold1;
    }
  else if ((lambdaold2==lambda) && (likep->iwlsweights_constant() == true) )
    {
    df = df_lambdaold2;
    }
  else
    {
    for(i=0;i<n;i++,workXX++)
      {
      df += (*workXX)/(*workXX+lambda);
      }

    df_lambdaold2 = df_lambdaold1;
    lambdaold2 = lambdaold1;
    df_lambdaold1 = df;
    lambdaold1 = lambda;

    }

  return df;
  }


void FULLCOND_random_gaussian::update_stepwise(double la)
  {
  lambda = la;
  }


ST::string FULLCOND_random_gaussian::get_effect(void)
  {

  ST::string h;

  if(randomslope)
    h = datanames[1] + "*" + datanames[0];
  else
    h = datanames[0];

  h = h + "(random,df=" + ST::doubletostring(compute_df(),6) + ",(lambda=" + ST::doubletostring(lambda,6) + "))";

  return h;

  }


void FULLCOND_random_gaussian::reset_effect(unsigned & pos)
  {

  update_linpred(false);

  unsigned i;
  double * work;
  work = beta.getV();
  for(i=0;i<nrpar;i++,work++)
    *work = 0.0;

  }

//------------------------------------------------------------------------------
//----------------- class FULLCOND_random_nongaussian --------------------------
//------------------------------------------------------------------------------

double FULLCOND_random_nongaussian::compute_df(void)
  {
  unsigned i;
  double df=0;
  double * workXX=XX.getV();
  unsigned n;

  if (randomslope && includefixed)
    {
    n = nrpar-1;
    df=1;
    }
  else
    n = nrpar;


  if ((lambdaold1==lambda) && (likep->iwlsweights_constant() == true) )
    {
    df = df_lambdaold1;
    }
  else if ((lambdaold2==lambda) && (likep->iwlsweights_constant() == true) )
    {
    df = df_lambdaold2;
    }
  else
    {
    for(i=0;i<n;i++,workXX++)
      {
      df += (*workXX)/(*workXX+lambda);
      }

    df_lambdaold2 = df_lambdaold1;
    lambdaold2 = lambdaold1;
    df_lambdaold1 = df;
    lambdaold1 = lambda;

    }

  return df;
  }


void FULLCOND_random_nongaussian::update_stepwise(double la)
  {
  lambda = la;
  }


ST::string FULLCOND_random_nongaussian::get_effect(void)
  {

  ST::string h;

  if(randomslope)
    h = datanames[1] + "*" + datanames[0];
  else
    h = datanames[0];

  h = h + "(random,df=" + ST::doubletostring(compute_df(),6) + ",(lambda=" + ST::doubletostring(lambda,6) + "))";

  return h;

  }


void FULLCOND_random_nongaussian::reset_effect(unsigned & pos)       
  {

  update_linpred(false);

  unsigned i;
  double * work;
  work = beta.getV();
  for(i=0;i<nrpar;i++,work++)
    *work = 0.0;

  }

*/
} // end: namespace MCMC



