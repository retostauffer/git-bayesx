
#include "cox.h"


namespace MCMC
{

//------------------------------------------------------------------------------
//----------------------- CLASS DISTRIBUTION_coxmodel ---------------------------
//------------------------------------------------------------------------------

DISTRIBUTION_coxmodel::DISTRIBUTION_coxmodel(MCMCoptions * o,
                                          const datamatrix & r,
                                          const datamatrix & t,
                                          const datamatrix & dbeg,
                                          const datamatrix & w)

   : DISTRIBUTION(o,r,w)
     {

     unsigned i;

     nrcat = 1; // ändern bei competing risk

     ti = t;

     int_ti=datamatrix(2*t.rows(),1,0.0);
     for(i=0;i<t.rows();i++)
     {
     int_ti(i,0) = t(i,0)-dbeg(i,0);
     int_ti(t.rows()+i,0)=0.0;
     }

     family = "cox";
     scale(0,0) = 1;
     scaleexi sting = false;

     }

DISTRIBUTION_coxmodel::DISTRIBUTION_coxmodel(

                   MCMCoptions * o, const datamatrix & r, const datamatrix & t,
                   const datamatrix & dbeg,
                        const datamatrix & w)
{
     unsigned i;

     nrcat = 1; // ändern bei competing risk

     ti = t;

     int_ti=datamatrix(2*t.rows(),1,0.0);
     for(i=0;i<t.rows();i++)
     {
     int_ti(i,0) = t(i,0)-dbeg(i,0);
     int_ti(t.rows()+i,0)=0.0;
     }

     family = "cox";
     scale(0,0) = 1;
     scaleexisting = false;

}



double DISTRIBUTION_coxmodel::compute_weight(double * linpred, double * weight,const int & i,
                                             const unsigned & col) const
  {
  return  *weight * exp(*linpred)* *(int_ti.getV()+i);
  }



void  DISTRIBUTION_coxmodel::tilde_y(datamatrix & tildey,datamatrix & m,const unsigned & col,
                                     const bool & current,const datamatrix & w)
 {
  unsigned i;

  double * workspline = m.getV();
  double * workresponse = response.getV();
  double * workweight = w.getV();
  double * ywork = tildey.getV();

  for (i=0;i<nrobs;i++,workspline++,ywork++,workresponse++,workweight++)
    {
    if(*workweight == 0.0)
      *ywork = 0.0;  
    else
      *ywork = *workspline + *workresponse/(*workweight)-1.0;
    }

 }


void DISTRIBUTION_coxmodel::compute_iwls(void)
  {
  unsigned i;

  double * linpred = linearpred.getV();
  double * workresponse = response.getV();
  double * workweightiwls = weightiwls.getV();
  double * ywork = tildey.getV();
  double * int_ti_work = int_ti.getV();

  for (i=0;i<nrobs;i++,linpred++,ywork++,workresponse++,workweightiwls++,int_ti_work++)
    {
    *workweightiwls = exp(*linpred)* *int_ti_work;
    *ywork = *linpred + *workresponse/(*workweightiwls)-1.0;
    }

  }


void DISTRIBUTION_coxmodel::compute_IWLS_weight_tildey(double * response,double * linpred,
                              double * weight,const int & i,
                              double * weightiwls,double * tildey,
                              const unsigned & col)
  {
  *weightiwls = exp(*linpred)* *(int_ti.getV() + i);
  *tildey = *response/(*weightiwls)-1.0;
  }


double DISTRIBUTION_coxmodel::compute_IWLS(double * response,double * linpred,double * weight,const int & i,
                              double * weightiwls, double * tildey,
                              bool weightyes,const unsigned & col)
  {
  if (weightyes)
    *weightiwls = exp(*linpred)* *(int_ti.getV() + i);

  *tildey = *response/(*weightiwls)-1.0;

  return *weight * (*response * (*linpred) - *weightiwls);
  }


void DISTRIBUTION_coxmodel::outoptions(void)
  {
  DISTRIBUTION::outoptions();
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTRIBUTION_coxmodel::update(void)
  {
  DISTRIBUTION::update();
  }


bool DISTRIBUTION_coxmodel::posteriormode(void)
  {
  return true;
  }

/*
void DISTRIBUTION_coxmodel::update_predict(void)
  {

  }
*/

void DISTRIBUTION_coxmodel::outresults(void)
  {
  DISTRIBUTION::outresults();
  }

void DISTRIBUTION_coxmodel::compute_deviance(const double * response, const double * weight,
                           const double * mu, double * deviance, double * deviancesat,
                           const datamatrix & scale, const int & i) const
    {

    double help = *mu * *(int_ti.getV()+i);

    *deviance = -2.0 * *weight *( *response * log(*mu) - help );
    *deviancesat = -2.0* *weight * (*response - help + *response * log(help));

    }


} // END: namespace MCMC

//---------------------------------------------------------------------------
#pragma package(smart_init)








