
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
  return *this;
  }


void FULLCOND_variance_nonp_vector::update(void)
  {
  acceptance++;

  // Hier kommen die Update-Schritte fuer die Varianzen hin.

  vector<double> vars(beta.rows());
  unsigned i;
  for(i=0; i<beta.rows(); i++)
    vars[i] = beta(i,0);
  Cp->update_variances(vars);

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





