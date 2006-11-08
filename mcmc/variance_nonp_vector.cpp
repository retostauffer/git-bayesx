
#include "first.h"

#include"variance_nonp_vector.h"

namespace MCMC
{

FULLCOND_variance_nonp_vector::FULLCOND_variance_nonp_vector(MCMCoptions * o,
                         FULLCOND_const * p,DISTRIBUTION * d,
                         const vector<double> & a, const vector<double> & b,
                         const ST::string & ti, const ST::string & fp,
                         const ST::string & fr,
                         const unsigned & c)
                         : FULLCOND(o,datamatrix(1,1),ti,1,1,fp)
    {
    fctype = MCMC::variance;
    update_sigma2 = true;
    column = c;
    pathresults = fr;
    Cp = p;
    distrp = d;
    a_invgamma = a;
    b_invgamma = b;
//    priorassumptions.push_back("Inverse gamma prior for variance component with hyperparameters a="
//    +ST::doubletostring(a,6)+ " and b=" + ST::doubletostring(b,6) );
    priorassumptions.push_back("\\\\");

    ST::string path = samplepath.substr(0,samplepath.length()-4)+"_lasso.raw";
    fc_lasso = FULLCOND(o,datamatrix(1,1),Cp->get_title()+"_lasso",1,1,path);
    fc_lasso.setflags(MCMC::norelchange | MCMC::nooutput);

    lasso=1;

    setbeta(Cp->get_variances());
    }

FULLCOND_variance_nonp_vector::FULLCOND_variance_nonp_vector(const FULLCOND_variance_nonp_vector & t)
    : FULLCOND(FULLCOND(t))
  {
  update_sigma2 = t.update_sigma2;
  column = t.column;
  pathresults = t.pathresults;
  Cp = t.Cp;
  distrp = t.distrp;
  a_invgamma = t.a_invgamma;
  b_invgamma = t.b_invgamma;
  fc_lasso = t.fc_lasso;
  lasso = t.lasso;
  }


const FULLCOND_variance_nonp_vector & FULLCOND_variance_nonp_vector::operator=(
                                        const FULLCOND_variance_nonp_vector & t)
  {
  if (this == &t)
    return *this;
  FULLCOND::operator=(FULLCOND(t));
  update_sigma2 = t.update_sigma2;
  column = t.column;
  pathresults = t.pathresults;
  Cp = t.Cp;
  distrp = t.distrp;
  a_invgamma = t.a_invgamma;
  b_invgamma = t.b_invgamma;
  fc_lasso = t.fc_lasso;
  lasso = t.lasso;
  return *this;
  }


void FULLCOND_variance_nonp_vector::update(void)
  {
  acceptance++;

  // Hier kommen die Update-Schritte fuer die Varianzen hin.

  // Ergebnisse zu den Regressionskoeffizienten zurueckschreiben
  Cp->update_variances(beta);

  // Lasso-Parameter abspeichern
  double * lassop = fc_lasso.getbetapointer();
  *lassop = lasso;
  fc_lasso.update();

  FULLCOND::update();
  }


void FULLCOND_variance_nonp_vector::outresults(void)
  {
  FULLCOND::outresults();

  ST::string l1 = ST::doubletostring(lower1,4);
  ST::string l2 = ST::doubletostring(lower2,4);
  ST::string u1 = ST::doubletostring(upper1,4);
  ST::string u2 = ST::doubletostring(upper2,4);

  ST::string nl1 = ST::doubletostring(lower1,4);
  ST::string nl2 = ST::doubletostring(lower2,4);
  ST::string nu1 = ST::doubletostring(upper1,4);
  ST::string nu2 = ST::doubletostring(upper2,4);
  nl1 = nl1.replaceallsigns('.','p');
  nl2 = nl2.replaceallsigns('.','p');
  nu1 = nu1.replaceallsigns('.','p');
  nu2 = nu2.replaceallsigns('.','p');

  ST::string vstr;

  ofstream ou(pathresults.strtochar());

  unsigned i;
  ou << "pmean  pstddev  pqu"  << nl1 << "   pqu" << nl2 << "  pmed pqu" <<
  nu1 << "   pqu" << nu2 << "  pmin  pmax" << endl;
  for(i=0; i<beta.rows(); i++)
    {
    ou << betamean(i,0) << "  ";
    ou << (betavar(i,0)<0.0?0.0:sqrt(betavar(i,0))) << "  ";
    ou << betaqu_l1_lower(i,0) << "  ";
    ou << betaqu_l2_lower(i,0) << "  ";
    ou << betaqu50(i,0) << "  ";
    ou << betaqu_l2_upper(i,0) << "  ";
    ou << betaqu_l1_upper(i,0) << "  ";
    ou << betamin(i,0) << "  ";
    ou << betamax(i,0) << "  " << endl;
    }

  optionsp->out("  Results for the variance components are stored in file\n");
  optionsp->out("  " + pathresults + "\n");

  optionsp->out("\n");

  outresults_lasso();
  }

void FULLCOND_variance_nonp_vector::outresults_lasso(void)
  {

  fc_lasso.outresults();

  ST::string l1 = ST::doubletostring(lower1,4);
  ST::string l2 = ST::doubletostring(lower2,4);
  ST::string u1 = ST::doubletostring(upper1,4);
  ST::string u2 = ST::doubletostring(upper2,4);

  ST::string nl1 = ST::doubletostring(lower1,4);
  ST::string nl2 = ST::doubletostring(lower2,4);
  ST::string nu1 = ST::doubletostring(upper1,4);
  ST::string nu2 = ST::doubletostring(upper2,4);
  nl1 = nl1.replaceallsigns('.','p');
  nl2 = nl2.replaceallsigns('.','p');
  nu1 = nu1.replaceallsigns('.','p');
  nu2 = nu2.replaceallsigns('.','p');

  ST::string vstr;

  optionsp->out("\n");
  optionsp->out("  f_ridge_lasso \n");
  optionsp->out("\n");
  optionsp->out("\n");

  vstr = "  Mean:         ";
  vstr = vstr + ST::string(' ',20-vstr.length())
              + ST::doubletostring(fc_lasso.get_betamean(0,0),6);
  optionsp->out(vstr + ST::string(' ',40-vstr.length()) + "\n");

  vstr = "  Std. dev.:    ";
  optionsp->out(vstr + ST::string(' ',20-vstr.length()) +
    ST::doubletostring((fc_lasso.get_betavar(0,0)<0.0?0.0:sqrt(fc_lasso.get_betavar(0,0))),6) + "\n");

  vstr = "  " + l1 + "% Quantile: ";
  vstr = vstr + ST::string(' ',20-vstr.length())
              + ST::doubletostring(fc_lasso.get_beta_lower1(0,0),6);
  optionsp->out(vstr + ST::string(' ',40-vstr.length()) + "\n");

  vstr = "  " + l2 + "% Quantile: ";
  vstr = vstr + ST::string(' ',20-vstr.length())
              + ST::doubletostring(fc_lasso.get_beta_lower2(0,0),6);
  optionsp->out(vstr + ST::string(' ',40-vstr.length()) + "\n");

  vstr = "  50% Quantile: ";
  vstr = vstr + ST::string(' ',20-vstr.length())
              + ST::doubletostring(fc_lasso.get_betaqu50(0,0),6);
  optionsp->out(vstr + ST::string(' ',40-vstr.length()) + "\n");

  vstr = "  " + u1 + "% Quantile: ";
  vstr = vstr + ST::string(' ',20-vstr.length())
              + ST::doubletostring(fc_lasso.get_beta_upper2(0,0),6);
  optionsp->out(vstr + ST::string(' ',40-vstr.length()) + "\n");

  vstr = "  " + u2 + "% Quantile: ";
  vstr = vstr + ST::string(' ',20-vstr.length())
              + ST::doubletostring(fc_lasso.get_beta_upper1(0,0),6);
  optionsp->out(vstr + ST::string(' ',40-vstr.length()) + "\n");

  optionsp->out("\n");

  ST::string lasso_pathresults = pathresults.substr(0,pathresults.length()-7) + "lasso.res";

  ofstream ou(lasso_pathresults.strtochar());

  ou << "pmean  pstddev  pqu"  << nl1 << "   pqu" << nl2 << "  pmed pqu" <<
  nu1 << "   pqu" << nu2 << "  pmin  pmax" << endl;
  ou << fc_lasso.get_betamean(0,0) << "  ";
  ou << (fc_lasso.get_betavar(0,0)<0.0?0.0:sqrt(fc_lasso.get_betavar(0,0))) << "  ";
  ou << fc_lasso.get_beta_lower1(0,0) << "  ";
  ou << fc_lasso.get_beta_lower2(0,0) << "  ";
  ou << fc_lasso.get_betaqu50(0,0) << "  ";
  ou << fc_lasso.get_beta_upper2(0,0) << "  ";
  ou << fc_lasso.get_beta_upper1(0,0) << "  ";
  ou << fc_lasso.get_betamin(0,0) << "  ";
  ou << fc_lasso.get_betamax(0,0) << "  " << endl;

  optionsp->out("  Results for the lasso parameter are also stored in file\n");
  optionsp->out("  " + lasso_pathresults + "\n");

  optionsp->out("\n");
  }


void FULLCOND_variance_nonp_vector::outoptions(void)
  {
  optionsp->out("  Hyperprior a for variance parameter: " +
                   ST::doubletostring(a_invgamma[0]) + "\n" );
  optionsp->out("  Hyperprior b for variance parameter: " +
                   ST::doubletostring(b_invgamma[0]) + "\n" );
  optionsp->out("\n");
  }


} // end: namespace MCMC





