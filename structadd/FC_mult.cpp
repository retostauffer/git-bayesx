
#include "FC_mult.h"


//------------------------------------------------------------------------------
//----------------- CLASS: FC_mult implementation of member functions ----------
//------------------------------------------------------------------------------


namespace MCMC
{

void FC_mult::set_effectp(DESIGN * d,FC_nonp * fp)
  {
  FCnp = fp;
  dp2 = d;
  }

void FC_mult::set_intp(DESIGN * d,FC_nonp * fp)
  {
  dp1 = d;
  FCnp2 = fp;
  }


void FC_mult::set_multeffects(MASTER_OBJ * mp,GENERAL_OPTIONS * o,
                              const ST::string & t, const ST::string & fp,
                              bool sm,bool sstore,bool meane)
  {

  masterp = mp;

  unsigned rows = dp1->Zout.rows()*dp2->Zout.rows();

  samplemult = sm;
  compmeaneffect = meane;

  if (samplemult)
    {
    FCmulteffect = FC(o,t,rows,1,fp,sstore);
    if (compmeaneffect==true)
      FCmulteffect_mean = FC(o,"",rows,1,fp,sstore);
    }

  }



FC_mult::FC_mult(void)
  {
  nosamples = true;
  samplemult = false;
  compmeaneffect = false;
  multexp=false;
  }


FC_mult::FC_mult(bool reu,bool mexp)
     : FC()
  {
  multexp = mexp;
  nosamples = true;
  samplemult = false;
  compmeaneffect = false;
  RE_update = reu;
  }


FC_mult::FC_mult(const FC_mult & m)
  : FC(FC(m))
  {
  masterp = m.masterp;
  multexp = m.multexp;
  FCmulteffect = m.FCmulteffect;
  FCmulteffect_mean = m.FCmulteffect_mean;
  samplemult = m.samplemult;
  compmeaneffect = m.compmeaneffect;
  dp1 = m.dp1;
  dp2 = m.dp2;
  FCnp = m.FCnp;
  FCnp2 = m.FCnp2;
  RE_update = m.RE_update;
  effect = m.effect;
  }



const FC_mult & FC_mult::operator=(const FC_mult & m)
  {

  if (this==&m)
	 return *this;
  FC::operator=(FC(m));
  masterp = m.masterp;
  multexp = m.multexp;
  FCmulteffect = m.FCmulteffect;
  FCmulteffect_mean = m.FCmulteffect_mean;
  samplemult = m.samplemult;
  compmeaneffect = m.compmeaneffect;
  dp1 = m.dp1;
  dp2 = m.dp2;
  FCnp = m.FCnp;
  FCnp2 = m.FCnp2;
  RE_update = m.RE_update;
  effect = m.effect;
  return *this;
  }


void FC_mult::update(void)
  {
  double add;

  if (RE_update)
    {
    add=0;
    }
  else
    {
    add=1;
    }

  dp2->compute_effect(effect,FCnp->beta,MCMC::Function);
  dp1->set_intvar(effect,add);

  if ((RE_update==false) && (samplemult))
    {
    update_multeffect();
    FCmulteffect.acceptance++;
    FCmulteffect.update();
    if (compmeaneffect)
      {
      FCmulteffect_mean.acceptance++;
      FCmulteffect_mean.update();
      }
    }

  }


void FC_mult::update_multeffect(void)
  {
  /*
  FCnp RE
  dp2 RE

  FCnp2 nonl
  dp1 nonl
  */

  unsigned i,j;
  double * mebetap = FCmulteffect.beta.getV();
  double * FCnpbetap = FCnp->beta.getV();
  double * FCnp2betap;

  // TEST
  //  ofstream out("c:\\bayesx\\testh\\results\\beta_RE.res");
  //  FCnp->beta.prettyPrint(out);
  // TEST

  // TEST
  // ofstream out2("c:\\bayesx\\testh\\results\\betaspline.res");
  // FCnp2->beta.prettyPrint(out2);
  // TEST

  if (compmeaneffect)
    {
    double meanlin = masterp->level1_likep->meaneffect-FCnp2->meaneffect;
    double lin;
    double * FCmp = FCmulteffect_mean.beta.getV();
    for (i=0;i<FCnp->beta.rows();i++,FCnpbetap++)
      {
      FCnp2betap = FCnp2->beta.getV();
      for (j=0;j<FCnp2->beta.rows();j++,mebetap++,FCnp2betap++,FCmp++)
        {
        *mebetap = (*FCnpbetap+1)*(*FCnp2betap);
        lin = meanlin + *mebetap;
        masterp->level1_likep->compute_mu(&lin,FCmp);
        }
      }
    }
  else
    {
    for (i=0;i<FCnp->beta.rows();i++,FCnpbetap++)
      {
      FCnp2betap = FCnp2->beta.getV();
      for (j=0;j<FCnp2->beta.rows();j++,mebetap++,FCnp2betap++)
        *mebetap = (*FCnpbetap+1)*(*FCnp2betap);
      }

    }

  }


bool FC_mult::posteriormode(void)
  {

  if (multexp==false)
    {
    double add;

    if (RE_update)
      {
      add=0;
      }
    else
      {
      add=1;
      }

    dp2->compute_effect(effect,FCnp->beta,MCMC::Function);
    dp1->set_intvar(effect,add);

    if ((RE_update==false) && (samplemult))
      {
      update_multeffect();
      bool h = FCmulteffect.posteriormode();

      if (compmeaneffect)
        {
        h = FCmulteffect_mean.posteriormode();
        }

      }
    }
  else
    {

    /*
    FCnp RE
    dp2 RE

    FCnp2 nonl
    dp1 nonl
    */


    double add=0;


    dp2->compute_effect(effect,FCnp->beta,MCMC::Function);

    if (RE_update)
      {
      }
    else
      {
      unsigned i;
      double * effectp = effect.getV();
      for (i=0;i<effect.rows();i++,effectp++)
        *effectp = exp(*effectp);
      }


    dp1->set_intvar(effect,add);

    if ((RE_update==false) && (samplemult))
      {
      update_multeffect();
      bool h = FCmulteffect.posteriormode();

      if (compmeaneffect)
        {
        h = FCmulteffect_mean.posteriormode();
        }

      }

    }

  return true;
  }


void FC_mult::outresults(ofstream & out_stata, ofstream & out_R,
                         const ST::string & pathresults)
  {
  if ((RE_update==false) && (samplemult))
    {
    FCmulteffect.outresults(out_stata,out_R,"");

    if (compmeaneffect)
      FCmulteffect_mean.outresults(out_stata,out_R,"");

    if (pathresults.isvalidfile() != 1)
      {

      FCmulteffect.optionsp->out("    Results are stored in file\n");
      FCmulteffect.optionsp->out("    " +  pathresults + "\n");
      FCmulteffect.optionsp->out("\n");

      ofstream outres(pathresults.strtochar());

      FCmulteffect.optionsp->out("\n");

      unsigned i,j;

      ST::string l1 = ST::doubletostring(FCmulteffect.optionsp->lower1,4);
      ST::string l2 = ST::doubletostring(FCmulteffect.optionsp->lower2,4);
      ST::string u1 = ST::doubletostring(FCmulteffect.optionsp->upper1,4);
      ST::string u2 = ST::doubletostring(FCmulteffect.optionsp->upper2,4);
      l1 = l1.replaceallsigns('.','p');
      l2 = l2.replaceallsigns('.','p');
      u1 = u1.replaceallsigns('.','p');
      u2 = u2.replaceallsigns('.','p');

      outres << "intnr" << "   ";

      /*
      FCnp RE
      dp2 RE

      FCnp2 nonl
      dp1 nonl
     */


      outres << dp1->datanames[0] << "   ";     
      outres << dp2->datanames[0] << "   ";


      outres << "pmean   ";

      if (FCmulteffect.optionsp->samplesize > 1)
        {
        outres << "pqu"  << l1  << "   ";
        outres << "pqu"  << l2  << "   ";
        outres << "pmed   ";
        outres << "pqu"  << u1  << "   ";
        outres << "pqu"  << u2  << "   ";
        outres << "pcat" << FCmulteffect.optionsp->level1 << "   ";
        outres << "pcat" << FCmulteffect.optionsp->level2 << "   ";
        }

      if (compmeaneffect)
        {
        outres << "pmean_mu   ";

        if (FCmulteffect.optionsp->samplesize > 1)
          {
          outres << "pqu"  << l1  << "_mu   ";
          outres << "pqu"  << l2  << "_mu   ";
          outres << "pmed_mu   ";
          outres << "pqu"  << u1  << "_mu   ";
          outres << "pqu"  << u2  << "_mu   ";
          }
        }


      outres << endl;

      double * workmean = FCmulteffect.betamean.getV();
      double * workbetaqu_l1_lower_p = FCmulteffect.betaqu_l1_lower.getV();
      double * workbetaqu_l2_lower_p = FCmulteffect.betaqu_l2_lower.getV();
      double * workbetaqu_l1_upper_p = FCmulteffect.betaqu_l1_upper.getV();
      double * workbetaqu_l2_upper_p = FCmulteffect.betaqu_l2_upper.getV();
      double * workbetaqu50 = FCmulteffect.betaqu50.getV();

      double * workmean_mu;
      double * workbetaqu_l1_lower_p_mu;
      double * workbetaqu_l2_lower_p_mu;
      double * workbetaqu_l1_upper_p_mu;
      double * workbetaqu_l2_upper_p_mu;
      double * workbetaqu50_mu;

      if (compmeaneffect)
        {
        workmean_mu = FCmulteffect_mean.betamean.getV();
        workbetaqu_l1_lower_p_mu = FCmulteffect_mean.betaqu_l1_lower.getV();
        workbetaqu_l2_lower_p_mu = FCmulteffect_mean.betaqu_l2_lower.getV();
        workbetaqu_l1_upper_p_mu = FCmulteffect_mean.betaqu_l1_upper.getV();
        workbetaqu_l2_upper_p_mu = FCmulteffect_mean.betaqu_l2_upper.getV();
        workbetaqu50_mu = FCmulteffect_mean.betaqu50.getV();
        }

      for (i=0;i<FCnp->beta.rows();i++)
        {
        for (j=0;j<FCnp2->beta.rows();j++,workmean++,workbetaqu_l1_lower_p++,
                              workbetaqu_l2_lower_p++,workbetaqu50++,
                              workbetaqu_l1_upper_p++,workbetaqu_l2_upper_p++)
          {
          outres << (i+1) << "   ";
          outres << dp2->effectvalues[i] << "   ";
          outres << dp1->effectvalues[j] << "   ";
          outres << *workmean << "   ";

          if (FCmulteffect.optionsp->samplesize > 1)
            {
            outres << *workbetaqu_l1_lower_p << "   ";
            outres << *workbetaqu_l2_lower_p << "   ";
            outres << *workbetaqu50 << "   ";
            outres << *workbetaqu_l2_upper_p << "   ";
            outres << *workbetaqu_l1_upper_p << "   ";

            if (*workbetaqu_l1_lower_p > 0)
              outres << 1 << "   ";
            else if (*workbetaqu_l1_upper_p < 0)
              outres << -1 << "   ";
            else
              outres << 0 << "   ";

            if (*workbetaqu_l2_lower_p > 0)
              outres << 1 << "   ";
            else if (*workbetaqu_l2_upper_p < 0)
              outres << -1 << "   ";
            else
              outres << 0 << "   ";
            }


          if (compmeaneffect)
            {

            outres << *workmean_mu << "   ";

            if (FCmulteffect.optionsp->samplesize > 1)
              {
              outres << *workbetaqu_l1_lower_p_mu << "   ";
              outres << *workbetaqu_l2_lower_p_mu << "   ";
              outres << *workbetaqu50_mu << "   ";
              outres << *workbetaqu_l2_upper_p_mu << "   ";
              outres << *workbetaqu_l1_upper_p_mu << "   ";

              }

            if ( !( (i == FCnp->beta.rows()-1)  && (j==FCnp2->beta.rows()-1) ) )
              {
              workmean_mu++;
              workbetaqu_l1_lower_p_mu++;
              workbetaqu_l2_lower_p_mu++;
              workbetaqu50_mu++;
              workbetaqu_l1_upper_p_mu++;
              workbetaqu_l2_upper_p_mu++;
              }

            }


          outres << endl;
          }
        }
      }
    }
  }


void FC_mult::reset(void)
  {

  }


} // end: namespace MCMC



