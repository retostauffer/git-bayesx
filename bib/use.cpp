


#include"use.h"


//------------------------------------------------------------------------------
//-------------- CLASS use: implementation of member functions -----------------
//------------------------------------------------------------------------------


const use & use::operator=(const use & u)
  {
  if (this == &u)
	 return *this;
  errormessages = u.errormessages;
  notext = u.notext;
  usingtext = u.usingtext;
  return *this;
  }


//------------------------------------------------------------------------------
//----------- CLASS usePathRead: implementation of member functions ------------
//------------------------------------------------------------------------------


const usePathRead & usePathRead::operator=(const usePathRead & u)
  {
  if (this == &u)
	 return *this;
  errormessages = u.errormessages;
  path = u.path;
  notext = u.notext;
  return *this;
  }


void usePathRead::parse(const ST::string & usetext)
  {

  path = "";
  errormessages.clear();
  notext = true;

  if (usetext.length() > 0)
	 {
	 notext = false;
	 if (usetext.isexistingfile() == 1)
		errormessages.push_back(
		"ERROR: file " + usetext + " could not be opened for reading\n");
	 if (errormessages.empty())
		path = usetext;
	 }

  }


//------------------------------------------------------------------------------
//----------- CLASS usePathWrite: implementation of member functions -----------
//------------------------------------------------------------------------------


const usePathWrite & usePathWrite::operator=(const usePathWrite & u)
  {
  if (this == &u)
	 return *this;
  errormessages = u.errormessages;
  notext = u.notext;
  path = u.path;
  alreadyexisting = u.alreadyexisting;
  return *this;
  }


void usePathWrite::parse(const ST::string & usetext)
  {

  path = "";
  errormessages.clear();
  notext = true;

  if (usetext.length() > 0)
	 {
	 notext = false;
	 int k = usetext.isvalidfile();

	 if (k == 1)
		{
		errormessages.push_back("ERROR: file " + usetext +
		" could not be opened for writing\n");
		alreadyexisting = false;
		}
	 else if (k == 0)
		alreadyexisting = false;
	 else
		alreadyexisting = true;

	 if (errormessages.empty())
		path = usetext;
	 }

  }


//------------------------------------------------------------------------------
//------------ CLASS useDataset: implementation of member functions ------------
//------------------------------------------------------------------------------


useDataset::useDataset(const useDataset & u)
  {
  errormessages = u.errormessages;
  notext = u.notext;
  datasets = u.datasets;
  }


const useDataset & useDataset::operator=(const useDataset & u)
  {
  if (this == &u)
	 return *this;
  errormessages = u.errormessages;
  notext = u.notext;
  datasets = u.datasets;
  return *this;
  }


void useDataset::parse(const ST::string & usetext)
  {

  errormessages.clear();
  notext = true;

  if (usetext.length() > 0)
	 {
	 notext = false;
	 if (! datasets->empty())
		{
		int i = 0;
		bool existing = false;
		while ( (i < datasets->size()) && (existing == false) )
		  {
		  if (usetext == (*((*datasets)[i])).getname())
			 {
			 existing = true;
			 datasetpointer = (*datasets)[i];
			 }
		  i++;
		  }
		if (existing == false)
		  errormessages.push_back(
		  "ERROR: dataset " + usetext + " is not existing\n");

		}
	 else
		errormessages.push_back(
		"ERROR: dataset " + usetext + " is not existing\n");
	 }

  }

