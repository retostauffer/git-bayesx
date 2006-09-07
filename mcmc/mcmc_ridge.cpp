
#include "first.h"

#include "mcmc_ridge.h"

//------------------------------------------------------------------------------
//------------ CLASS: FULLCOND_ridge implementation of member functions --------
//------------------------------------------------------------------------------


namespace MCMC
{

//----------------------- CONSTRUCTORS, DESTRUCTOR -----------------------------

FULLCOND_ridge::FULLCOND_ridge(MCMCoptions * o, DISTRIBUTION * dp,
                const datamatrix & d, const ST::string & t,
                const ST::string & fs, const ST::string & fr,
                const vector<double> & vars, const unsigned & c)
  : FULLCOND(o,d,t,d.cols(),1,fs)
  {
  unsigned i;//, j;

  variances = vars;
  nrpar = d.cols();
  likep = dp;

  linold = datamatrix(likep->get_nrobs(),1,0);
  mu1 = datamatrix(likep->get_nrobs(),1,0);

  column = c;

/*  XX = datamatrix(d.cols(),1,0);
  for(i=0; i<d.rows(); i++)
    {
    for(j=0; j<d.cols(); j++)
      {
      XX(j,0) += d(i,j)*d(i,j);
      }
    }*/

  datamatrix XX = d.transposed()*d;
  for(i=0; i<nrpar; i++)
    {
    XX(i,i) = XX(i,i) + 1/variances[i];
    }
  X2 = XX.cinverse()*(d.transposed());
  X1 = (XX.cinverse()).root();

  setbeta(d.cols(),1,0);

  pathresult = fr;
  pathcurrent = fr;
  }

// COPY CONSTRUCTOR

FULLCOND_ridge::FULLCOND_ridge(const FULLCOND_ridge & m)
  : FULLCOND(FULLCOND(m))
  {
  variances = m.variances;
  X1 = m.X1;
  X2 = m.X2;
  linold = m.linold;
  mu1 = m.mu1;
  likep = m.likep;
  }

// OVERLOADED ASSIGNMENT OPERATOR

const FULLCOND_ridge & FULLCOND_ridge::operator=(const FULLCOND_ridge & m)
  {
  if (this == &m)
    return *this;
  FULLCOND::operator=(FULLCOND(m));

  variances = m.variances;
  X1 = m.X1;
  X2 = m.X2;
  linold = m.linold;
  mu1 = m.mu1;
  likep = m.likep;

  return *this;
  }

//-------------------------- UPDATE and related methods-------------------------

  // FUNCTION: update
  // TASK: - stores sampled parameters in file 'samplepath'
  //         storing order: first row, second row, ...

void FULLCOND_ridge::update(void)
  {
  FULLCOND::update();

  likep->substr_linearpred_m(linold,column);

  likep->compute_respminuslinpred(mu1,column);

  beta.mult(X2,mu1);
  beta+= sqrt(likep->get_scale(column))*X1*rand_normvek(nrpar);

  linold.mult(data,beta);
  likep->add_linearpred_m(linold,column);

  acceptance++;

  transform = likep->get_trmult(column);
  
  if  (optionsp->get_nriter() == optionsp->get_iterations())
    {
    FULLCOND::outresults();
    }
  }

  // FUNCTION: outresults
  // TASK: - write results to output window and files

void FULLCOND_ridge::outresults(void)
  {
  FULLCOND::outresults();

  ofstream outp(pathcurrent.strtochar());

  ST::string l1 = ST::doubletostring(lower1,4);
  ST::string l2 = ST::doubletostring(lower2,4);
  ST::string u1 = ST::doubletostring(upper1,4);
  ST::string u2 = ST::doubletostring(upper2,4);
  l1 = l1.replaceallsigns('.','p');
  l2 = l2.replaceallsigns('.','p');
  u1 = u1.replaceallsigns('.','p');
  u2 = u2.replaceallsigns('.','p');

  vector<ST::string> resultstable(6);

  outp << "paramnr varname pmean pstd pqu" << l1 << " pqu" << l2 <<
            " pmed pqu" << u1 << " pqu" << u2 << " pcat" << level1
            << " pcat" << level2 << endl;

  unsigned i;

  optionsp->out("\n");

  ST::string l;
  int maxvarnamelength = 0;
  int len;

  for(i=0;i<nrpar;i++)
    {
    len = datanames[i].length();
    if (len > maxvarnamelength)
      maxvarnamelength = len;
    }

  if (maxvarnamelength>10)
    l = ST::string(' ',maxvarnamelength-6);
  else
    l = "  ";

    ST::string help =  ST::doubletostring(lower1,4) + "% quant.";
    ST::string levell = help + ST::string(' ',15-help.length());
    help = ST::doubletostring(upper2,4) + "% quant.";
    ST::string levelu = help + ST::string(' ',15-help.length());

    optionsp->out("  Variable" + l +
                  "mean           " +
                  "Std. Dev.      " +
                  levell +
                  "median         " +
                  levelu + "\n");

    ST::string mean;
    ST::string std;
    ST::string qu10;
    ST::string qu50;
    ST::string qu90;

    double m,stddouble;

    unsigned nsp;

    for (i=0;i<nrpar;i++)
      {

      if (maxvarnamelength  > 10)
        nsp = 2+maxvarnamelength-datanames[i].length();
      else
        nsp = 10-datanames[i].length();

      m= betamean(i,0);

      if (betavar(i,0) == 0)
        stddouble = 0;
      else
        stddouble = sqrt(betavar(i,0));

      outp << (i+1) << "   ";
      outp << datanames[i] << "   ";
      outp << m << "   ";
      outp << stddouble << "   ";
      outp << betaqu_l1_lower(i,0) << "   ";
      outp << betaqu_l2_lower(i,0) << "   ";
      outp << betaqu50(i,0) << "   ";
      outp << betaqu_l2_upper(i,0) << "   ";
      outp << betaqu_l1_upper(i,0) << "   ";
      if (betaqu_l1_lower(i,0) > 0)
        outp << "1   ";
      else if (betaqu_l1_upper(i,0) < 0)
        outp << "-1   ";
      else
        outp << "0   ";

      if (betaqu_l2_lower(i,0) > 0)
        outp << "1   ";
      else if (betaqu_l2_upper(i,0) < 0)
        outp << "-1   ";
      else
        outp << "0   ";

      outp << endl;

      optionsp->out(ST::outresults(nsp,datanames[i],m,
                      stddouble,betaqu_l1_lower(i,0),
                      betaqu50(i,0),betaqu_l1_upper(i,0)) + "\n");


      char hchar = '_';
      ST::string hstring = "\\_";
      resultstable[0] = datanames[i].insert_string_char(hchar,hstring);
      resultstable[1] = ST::doubletostring(m,6);
      resultstable[2] = ST::doubletostring(stddouble,6);
      resultstable[3] = ST::doubletostring(betaqu_l1_lower(i,0),6);
      resultstable[4] = ST::doubletostring(betaqu50(i,0),6);
      resultstable[5] = ST::doubletostring(betaqu_l1_upper(i,0),6);

      results_latex.push_back(ST::make_latextable(resultstable));

      }

    optionsp->out("\n");

    optionsp->out("  Results for shrinked effects are also stored in file\n");
    optionsp->out("  " + pathcurrent + "\n");

    optionsp->out("\n");

  }


} // end: namespace MCMC

