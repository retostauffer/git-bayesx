//---------------------------------------------------------------------------
#include <vcl.h>
#pragma hdrstop

#include <Registry.hpp>
#include <FileCtrl.hpp>

#include "statwin_haupt.h"
#include "StatwinFrame.h"
#include "command.h"
#include "StatResults.h"
#include "StatReview.h"
#include "StatObjects.h"
#include "dir.h"
//---------------------------------------------------------------------------
#pragma package(smart_init)
#pragma resource "*.dfm"
Thauptformular *hauptformular;
//---------------------------------------------------------------------------

__fastcall Thauptformular::Thauptformular(TComponent* Owner)
    : TForm(Owner)
  {

  char path[100];
  getcurdir(0,path);
  char disk = 'A'+getdisk();
  ST::string d(disk,1);

  defaultpath = d + ":\\" + path;

  bool error = false;

  if(DirectoryExists((defaultpath+"\\temp").strtochar()))
    {
    AnsiString temp = (defaultpath+"\\temp\\test").strtochar();
    ForceDirectories(temp);
    if(DirectoryExists(temp))
      {
      rmdir((defaultpath+"\\temp\\test").strtochar());
      }
    else
      {
      out("ERROR: No permission to write to " + defaultpath + "\\temp\n");
      error = true;
      }
    }
  else
    {
    AnsiString temp = (defaultpath+"\\temp").strtochar();
    ForceDirectories(temp);
    if(!DirectoryExists(temp))
      {
      error = true;
      }
    }

  if(DirectoryExists((defaultpath+"\\output").strtochar()))
    {
    AnsiString output = (defaultpath+"\\output\\test").strtochar();
    ForceDirectories(output);
    if(DirectoryExists(output))
      {
      rmdir((defaultpath+"\\output\\test").strtochar());
      }
    else
      {
      out("ERROR: No permission to write to " + defaultpath + "\\output\n");
      error = true;
      }
    }
  else
    {
    AnsiString output = (defaultpath+"\\output").strtochar();
    ForceDirectories(output);
    if(!DirectoryExists(output))
      {
      error = true;
      }
    }

  if(error==true)
    {
    out("ERROR: No permission to write to " + defaultpath + "\n");
    out("  Specify a new default directory using the defaultpath command\n");
    out("  Type for example: defaultpath=c:\\temp");
    }

  logfileopen = false;
  input = &cin;

  objecttyps.reserve(10);

  objecttyps.push_back("dataset");
#if defined(INCLUDE_MCMC)
  objecttyps.push_back("bayesreg");
#endif
  objecttyps.push_back("map");
#if defined(INCLUDE_DAG)
  objecttyps.push_back("dag");
#endif
#if defined(INCLUDE_REML)
  objecttyps.push_back("remlreg");
#endif
#if defined(INCLUDE_STEP)
  objecttyps.push_back("stepwisereg");
#endif
  objecttyps.push_back("mcmcreg");  


  delim = '\r';
  commandedit->Tag = false;

  }
//---------------------------------------------------------------------------


void Thauptformular::out(const ST::string & c)
  {
#if defined(BORLAND_OUTPUT_WINDOW)
  ST::string sh = c;
  sh = sh.replaceallsigns('\n',' ');
  if (!Frame->suppoutput)
    Results->ResultsRichEdit->Lines->Append(sh.strtochar());
  if (logout.is_open())
    logout << c << flush;
#else
  cout << c << flush;
  if (logout.is_open())
    logout << c << flush;
#endif
  }


void Thauptformular::out(const vector<ST::string> & m)
  {
  for (unsigned i=0;i<m.size();i++)
	 out(m[i]);
  }


void Thauptformular::dropobjects(ST::string name, ST::string type)
  {

  int recognized = 0;
  unsigned i=0;

  if (type == "dataset")
	 {
	 while ( (i < dataobjects.size()) && (recognized == 0) )
		{

		if ( name == dataobjects[i].getname())
		  {

		  dataobjects.erase(dataobjects.begin()+i,dataobjects.begin()+i+1);
		  recognized = 1;
		  }
		i++;
		}
	 } // end: type dataset
#if defined(INCLUDE_MCMC)
  else if (type == "bayesreg")
	 {
	 while ( (i < bayesregobjects.size()) && (recognized == 0) )
		{

		if ( name == bayesregobjects[i].getname())
		  {

		  bayesregobjects.erase(bayesregobjects.begin()+i,bayesregobjects.begin()+i+1);
		  recognized = 1;
		  }
		i++;
		}
	 }  // end: type bayesreg
#endif
  else if (type == "mcmcreg")
	 {
	 while ( (i < mcmcregobjects.size()) && (recognized == 0) )
		{

		if ( name == mcmcregobjects[i].getname())
		  {

		  mcmcregobjects.erase(mcmcregobjects.begin()+i,mcmcregobjects.begin()+i+1);
		  recognized = 1;
		  }
		i++;
		}
	 }  // end: type bayesreg
#if defined(INCLUDE_REML)
  else if (type == "remlreg")
	 {
	 while ( (i < remlregobjects.size()) && (recognized == 0) )
		{

		if ( name == remlregobjects[i].getname())
		  {

		  remlregobjects.erase(remlregobjects.begin()+i,remlregobjects.begin()+i+1);
		  recognized = 1;
		  }
		i++;
		}
	 }  // end: type remlreg
#endif
#if defined(INCLUDE_STEP)
  else if (type == "stepwisereg")
	 {
	 while ( (i < stepwiseregobjects.size()) && (recognized == 0) )
		{

		if ( name == stepwiseregobjects[i].getname())
		  {

		  stepwiseregobjects.erase(stepwiseregobjects.begin()+i,
                  stepwiseregobjects.begin()+i+1);
		  recognized = 1;
		  }
		i++;
		}
	 }  // end: type stepwisereg
#endif
  else if (type == "map")
	 {
	 while ( (i < mapobjects.size()) && (recognized == 0) )
		{

		if ( name == mapobjects[i].getname())
		  {

		  mapobjects.erase(mapobjects.begin()+i,mapobjects.begin()+i+1);
		  recognized = 1;
		  }
		i++;
		}
	 } // end: type map
#if defined(INCLUDE_DAG)
  else if (type == "dag")
     {
	 while ( (i < dagobjects.size()) && (recognized == 0) )
		{

		if ( name == dagobjects[i].getname())
		  {

		  dagobjects.erase(dagobjects.begin()+i,dagobjects.begin()+i+1);
		  recognized = 1;
		  }
		i++;
		}

     }
#endif

 /*
  else if (type == "diseasemap")
     {
	 while ( (i < diseasemapobjects.size()) && (recognized == 0) )
		{

		if ( name == diseasemapobjects[i].getname())
		  {

		  diseasemapobjects.erase(diseasemapobjects.begin()+i,
          diseasemapobjects.begin()+i+1);
		  recognized = 1;
		  }
		i++;
		}

     }
 */


  adjustobjects();

  }


bool Thauptformular::alreadyexisting(const ST::string & name)
  {

  unsigned i = 0;
  bool existing = false;
  while ( (i < objects.size()) && (existing == false) )
	 {
	 if (name == objects[i]->getname())
		existing = true;
	 i++;
	 }
  return existing;
  }


ST::string Thauptformular::create(const ST::string & in)
  {

  ST::string name;
  vector<ST::string> token = in.strtoken(" ");

  if (token.size() != 2)
	 {
	 errormessages.push_back("ERROR: invalid creation of new object\n");
	 return "";
	 }
  else // token.size() == 2
	 {

	 if (token[1].isvarname() == 1)
		{
		errormessages.push_back("ERROR: invalid object name\n");
		return "";
		}
	 else
		{
		if (alreadyexisting(token[1]) == true)
		  errormessages.push_back(
		  "ERROR: object " + token[1] + " is already existing\n");
		else
		  {

		  if (token[0] == "dataset")
			 {
			 dataobject newobject(token[1],&logout,input);
			 dataobjects.push_back(newobject);
			 }
#if defined(INCLUDE_MCMC)
		  else if (token[0] == "bayesreg")
			 {
			 bayesreg newobject(token[1],&logout,input,defaultpath,&objects);
			 bayesregobjects.push_back(newobject);
			 }
#endif
		  else if (token[0] == "mcmcreg")
			 {
			 superbayesreg newobject(token[1],&logout,input,defaultpath,&objects);
			 mcmcregobjects.push_back(newobject);
			 }
#if defined(INCLUDE_REML)
		  else if (token[0] == "remlreg")
			 {
			 remlreg newobject(token[1],&logout,input,defaultpath,&objects);
			 remlregobjects.push_back(newobject);
			 }
#endif             
#if defined(INCLUDE_STEP)
		  else if (token[0] == "stepwisereg")
			 {
			 stepwisereg newobject(token[1],&logout,input,defaultpath,&objects);
			 stepwiseregobjects.push_back(newobject);
			 }
#endif			 
		  else if (token[0] == "map")
			 {
			 mapobject newobject(token[1],&logout,input,defaultpath,&objects);
			 mapobjects.push_back(newobject);
			 }
#if defined(INCLUDE_DAG)
           else if (token[0] == "dag")
             {
             dagobject newobject(token[1],&logout,input,defaultpath,&objects);
             dagobjects.push_back(newobject);
             }
#endif
             
//           else if (token[0] == "diseasemap")
//             {
//             diseaseobj newobject(token[1],&logout,input,defaultpath,&objects);
//             diseasemapobjects.push_back(newobject);
//             }



		  adjustobjects();


		  }
		return token[1];
		}
	 } // end: if token.size() == 2
  }


void Thauptformular::parseexisting(const ST::string & objectname,const ST::string & com)
  {

  int recognized = 0;
  unsigned i=0;


  while ( (i < objects.size()) && (recognized == 0) )
	 {

	 if ( objectname == objects[i]->getname())
		{
		objects[i]->parse(com);
		errormessages = objects[i]->geterrormessages();
		recognized = 1;
		}
	 i++;
	 }

  if (recognized == 0)
	 errormessages.push_back(
		"ERROR: object " + objectname + " is not existing\n");

  }


void Thauptformular::adjustobjects(void)
  {
  objects.erase(objects.begin(),objects.end());
  unsigned i;

  for (i=0;i<dataobjects.size();i++)
	 objects.push_back(&dataobjects[i]);

#if defined(INCLUDE_MCMC)
  for (i=0;i<bayesregobjects.size();i++)
	 objects.push_back(&bayesregobjects[i]);
#endif

  for (i=0;i<mcmcregobjects.size();i++)
	 objects.push_back(&mcmcregobjects[i]);

#if defined(INCLUDE_REML)
  for (i=0;i<remlregobjects.size();i++)
	 objects.push_back(&remlregobjects[i]);
#endif     

#if defined(INCLUDE_STEP)
  for (i=0;i<stepwiseregobjects.size();i++)
	 objects.push_back(&stepwiseregobjects[i]);
#endif

  for (i=0;i<mapobjects.size();i++)
	 objects.push_back(&mapobjects[i]);

#if defined(INCLUDE_DAG)
  for (i=0;i<dagobjects.size();i++)
	 objects.push_back(&dagobjects[i]);
#endif

//  for (i=0;i<diseasemapobjects.size();i++)
//	 objects.push_back(&diseasemapobjects[i]);


  }


bool Thauptformular::parse(ST::string & in)
  {

  errormessages.clear();

  ST::string objectname;
  ST::string firsttoken = in.getFirstToken(" ,.=");
  int pointpos = in.checksign('.');

  if (firsttoken.length() > 0 && firsttoken[0] != '%')
	 {
	 if ( (firsttoken == "quit") || (firsttoken == "exit") )
		return true;
	 else if (firsttoken == "delimiter")
		{
		vector<ST::string> token = in.strtoken(" =");
		if (token.size() != 3)
		  errormessages.push_back("ERROR: invalid syntax\n");
		else if (token[1] != "=")
		  errormessages.push_back("ERROR: \"=\" expected\n");
		else
		  {
		  if (token[2] == "return")
			 delim = '\r';
		  else if (token[2] == ";")
            delim = token[2][0];
		  else
            errormessages.push_back("ERROR: invalid delimiter symbol\n");
		  }
		return false;
		} // end: delimiter
	 else if (firsttoken == "usefile")
		{

		vector<ST::string> token = in.strtoken(" ");
		if (token.size() < 2)
		  errormessages.push_back("ERROR: filename expected\n");
		else if (token.size() > 2)
		  errormessages.push_back("ERROR: invalid syntax\n");

		if (errormessages.empty())
		  {
		  ST::string path = token[1];
		  if (path.isexistingfile() == 1)
			 errormessages.push_back("ERROR: file " + path +
                         " could not be opened\n");
		  else
			 {
			 ST::string in;
			 ifstream infile;
			 input = &infile;
			 ST::open(infile,path);
			 while (! infile.eof())
				{

                Application->ProcessMessages();
                if (Frame->stop)
                  {
                  out("\n");
                  out("USER BREAK\n");
                  out("\n");
                  }
                if (Frame->stop)
                  break;

                if (Frame->pause)
                  {
                  out("\n");
                  out("PROGRAM PAUSED\n");
                  out("Click CONTINUE to proceed\n");
                  out("\n");

                  while (Frame->pause)
                   {
                   Application->ProcessMessages();
                   }

                   out("CONTINUED\n");
                   out("\n");

                   }

				if (delim != '\r')
                  {
                  bool end = false;
                  ST::string help="";
                  in = "";
                  while ( (! infile.eof()) && (end == false) )
                    {
                    ST::getline(infile,10000,help,'\n');
                    help = help.eatwhitespace();
                    help = help.eatallcarriagereturns();
                    if (help.length() > 0 && help[0] != '%')
                      {
                      in = in + " " + help;
                      if (help[help.length()-1] == ';')
                        end = true;
                      }
                    }

                  if (in.length() > 0)
                    {
                    if (in[in.length()-1] == delim)
                      in = in.deletesign(in.length()-1);
                    in = in.replaceallsigns('\n',' ');
                    }

//                  ST::getline(infile,10000,in,delim);

                  }
                else
                  ST::getline(infile,10000,in,'\n');

	            in = in.eatwhitespace();
                in = in.eatallcarriagereturns();

                if (in.length() > 0 && in[0] != '%')
                  {
				  out("> " + in + "\n");
                  Review->ReviewListBox->Items->Add(in.strtochar());

                  if(Review->ReviewListBox->Items->Count > 100)
                    Review->ReviewListBox->Items->Delete(0);

                  }
				parse(in);
				} // end: while (! infile.eof())

			 }

		  }

		input = &cin;

        if (!Frame->stop)
          {
		  out(errormessages);
          errormessages.clear();
          }

		return false;

		}
     else if (firsttoken == "saveoutput")
        {

		model m;
		simpleoption replace("replace",false);
        vector<ST::string> types;
        types.push_back("rtf");
        types.push_back("txt");
        stroption type("type",types,"rtf");
		optionlist saveoutoptions;
		saveoutoptions.push_back(&replace);
        saveoutoptions.push_back(&type);
		usePathWrite uw;

		command saveoutput("saveoutput",&m,&saveoutoptions,&uw,notallowed,notallowed,
        notallowed,notallowed,optional,required);

		saveoutput.parse(in);
		errormessages = saveoutput.geterrormessages();
		if (errormessages.empty())
		  {

		  if ((replace.getvalue() == false) && (uw.isexisting() == true))
            {
            out("ERROR: file " + uw.getPath() +
                                    " is already existing\n");
            }
		  else
            {

            if (type.getvalue() == "rtf")
              Results->ResultsRichEdit->PlainText = false;
            else
              Results->ResultsRichEdit->PlainText = true;

            Results->ResultsRichEdit->Lines->SaveToFile(uw.getPath().strtochar());
            Results->ResultsRichEdit->Modified = false;
            Results->ResultsRichEdit->Tag = true;
            Results->Caption = uw.getPath().strtochar();

            }
          }
        else
          {
          out(errormessages);
          errormessages.clear();
          }

        return false;
        }
     else if (firsttoken == "clearoutput")
       {
       Results->ResultsRichEdit->Clear();
       }
	 else if (firsttoken == "logopen")
		{
		model m;
		simpleoption replace("replace",false);
		optionlist logoptions;
		logoptions.push_back(&replace);
		usePathWrite uw;

		command logopen("logopen",&m,&logoptions,&uw,notallowed,notallowed,
        notallowed,notallowed,optional,required);

		logopen.parse(in);
		errormessages = logopen.geterrormessages();
		if (logfileopen == true)
		  errormessages.push_back("ERROR: log-file is already open\n");
		if (errormessages.empty())
		  {
		  logfileopen = true;
		  logfilepath = uw.getPath();
		  if ((replace.getvalue() == false) && (uw.isexisting() == true))
            {
            ST::open(logout,logfilepath,ios::app);
            }
		  else
            {
            ST::open(logout,logfilepath);
            }
		  }
		else
          {
		  out(errormessages);
          errormessages.clear();
          }
		return false;
		}  // end: logopen
	 else if (firsttoken == "logclose")
		{
		if (logfileopen == false)
		  {
		  errormessages.push_back("ERROR: currently no log-file open\n");
		  out(errormessages);
          errormessages.clear();
		  }
		else
		  {
		  logfileopen = false;
		  logout.close();
		  out("NOTE: log-file " + logfilepath + " closed\n");
		  }
		return false;
		}  // end: logclose
     else if (firsttoken == "defaultpath")
        {
        vector<ST::string> token = in.strtoken(" =");
        if (token.size() != 3)
          errormessages.push_back("ERROR: invalid syntax\n");
        else if (token[1] != "=")
          errormessages.push_back("ERROR: \"=\" expected\n");
        else if (!DirectoryExists(token[2].strtochar()))
          errormessages.push_back("ERROR: path " + token[2] + " does not exist\n");
        else
           {

           defaultpath = token[2];

           bool error = false;

           if(DirectoryExists((defaultpath+"\\temp").strtochar()))
             {
             AnsiString temp = (defaultpath+"\\temp\\test").strtochar();
             ForceDirectories(temp);
             if(DirectoryExists(temp))
               {
               rmdir((defaultpath+"\\temp\\test").strtochar());
               }
             else
               {
               out("ERROR: No permission to write to " + defaultpath + "\\temp\n");
               error = true;
               }
             }
           else
             {
             AnsiString temp = (defaultpath+"\\temp").strtochar();
             ForceDirectories(temp);
             if(!DirectoryExists(temp))
               {
               error = true;
               }
             }

           if(DirectoryExists((defaultpath+"\\output").strtochar()))
             {
             AnsiString output = (defaultpath+"\\output\\test").strtochar();
             ForceDirectories(output);
             if(DirectoryExists(output))
               {
               rmdir((defaultpath+"\\output\\test").strtochar());
               }
             else
               {
               out("ERROR: No permission to write to " + defaultpath + "\\output\n");
               error = true;
               }
             }
           else
             {
             AnsiString output = (defaultpath+"\\output").strtochar();
             ForceDirectories(output);
             if(!DirectoryExists(output))
               {
               error = true;
               }
             }

           if(error==true)
             {
             out("ERROR: No permission to write to " + defaultpath + "\n");
         	 out("  Specify a new default directory using the defaultpath command\n");
         	 out("  Type for example: defaultpath=c:\\temp");
             }

           }
	    out(errormessages);
        errormessages.clear();
        return false;
        }
	 else if (firsttoken == "drop")
		{

		modelStandard m;
		optionlist dropoptions;
		usePathWrite uw;

		command drop("drop",&m,&dropoptions,&uw,required,notallowed,
							 notallowed,notallowed,notallowed,notallowed);

		drop.parse(in);
		errormessages = drop.geterrormessages();
		vector<ST::string> objectnames =  m.getModelVarnamesAsVector();
		if (objectnames.size() == 0)
		  errormessages.push_back("ERROR: objectlist required\n");

		if (errormessages.empty())
		  {

		  int j;
		  for (j=0;j<objectnames.size();j++)
			 {

			 int recognized = 0;
			 int i=0;
			 while ( (i < objects.size()) && (recognized == 0) )
				{

				if ( objectnames[j] == objects[i]->getname())
				  {
				  ST::string type = objects[i]->gettype();
				  dropobjects(objectnames[j],type);
				  recognized = 1;
				  }
				i++;
				}

			  if (recognized == 0)
				 errormessages.push_back(
				 "ERROR: object " + objectnames[j] + " is not existing\n");

			 } // end: for (j=0;j<objectnames.size();j++)

		  }  // end: if (errormessages.empty())

		out(errormessages);
        errormessages.clear();
		return false;
		} // end: drop
	 else if (firsttoken.isinlist(objecttyps) >= 0)      // create a new object
		{


		if (pointpos == -1)
		  objectname = create(in);
		else
          objectname = create(in.substr(0,pointpos));


		if ( (errormessages.empty()) && (pointpos > 0) )
		  {
		  if (in.length()-1-pointpos <= 0)
			 errormessages.push_back("ERROR: invalid syntax\n");
		  else
			parseexisting(objectname,in.substr(pointpos+1,in.length()-1-pointpos));
		  }

		out(errormessages);
        errormessages.clear();
		return false;

		}               // end: create a new object
	 else     // existing object
		{
		if (pointpos != firsttoken.length())
		  errormessages.push_back("ERROR: invalid syntax\n");
		else
		  if (in.length() > pointpos+1)
			parseexisting(firsttoken,in.substr(pointpos+1,in.length()-pointpos-1));
		  else
			 errormessages.push_back("ERROR: invalid syntax\n");

		out(errormessages);
        errormessages.clear();
		return false;

		}

	 }  // end: if (firsttoken.length() > 0)
  else                                                  // empty command
    {
    if (in.length() > 0 && in[0] != '%')
      out("ERROR: invalid syntax\n");
    return false;
    }

  }


void __fastcall Thauptformular::commandeditKeyPress(TObject *Sender, char &Key)
{

  cmd = commandedit->Text.c_str();
  if (Key == delim)     // Return
    {
    commandedit->Tag = true;
    }

}
//---------------------------------------------------------------------------


void __fastcall Thauptformular::commandeditKeyUp(TObject *Sender,
      WORD &Key, TShiftState Shift)
{

  if (commandedit->Tag == true)
    {
//    ST::string inp = commandedit->Text.c_str();
    commandedit->Clear();
    Results->Show();
//    parsecommand(inp);
    parsecommand(cmd);
    hauptformular->Show();
    commandedit->Tag = false;
    }

}
//---------------------------------------------------------------------------


void Thauptformular::parsecommand(ST::string & inp)
  {

  Frame->stop = false;

  if (Frame->pause)
    {
    out("\n");
    out("PROGRAM PAUSED\n");
    out("Click CONTINUE to proceed\n");
    out("\n");

    while (Frame->pause)
      {
      Application->ProcessMessages();
      }

    out("CONTINUED\n");
    out("\n");

    }

  bool stop = false;

  inp = inp.deleteallsigns(delim);
  inp = inp.replaceallsigns('\n',' ');
  inp = inp.replaceallsigns('\r',' ');
  inp = inp.eatwhitespace();

  char* INP = inp.strtochar();

  if (strcmp(INP,"quit") != 0 && strcmp(INP,"exit") != 0)
    {
    if (inp.length() > 0 && inp[0] != '%' )
      {
      out("> " + inp + "\n");
      Review->ReviewListBox->Items->Add(INP);

      if(Review->ReviewListBox->Items->Count > 100)
        Review->ReviewListBox->Items->Delete(0);

      }
    }

  Frame->processruning = true;

  hauptformular->Enabled = false;
//  Review->Enabled = false;
  Objectbrowser->Enabled = false;
//  Frame->FMWReview->Enabled = false;
  Frame->FMWObjectbrowser->Enabled = false;
  Frame->FMWCommand->Enabled = false;
//  Frame->FileSave->Enabled = false;
//  Frame->FileSaveas->Enabled = false;
//  Frame->FilePrint->Enabled = false;
  Frame->FileOpen->Enabled = false;

  stop = parse(inp);

  Frame->processruning = false;

//  Frame->FileSave->Enabled = true;
//  Frame->FileSaveas->Enabled = true;
//  Frame->FilePrint->Enabled = true;
  Frame->FileOpen->Enabled = true;
//  Frame->FMWReview->Enabled = true;
  Frame->FMWObjectbrowser->Enabled = true;
  Frame->FMWCommand->Enabled = true;
  Objectbrowser->Enabled = true;
  hauptformular->Enabled = true;
//  Review->Enabled = true;

  if (stop)
    {
    Frame->writetoregistry();
    if(CheckOutputSaved() != IDCANCEL)
      {
      Application->Terminate();
      }
    }

  }
//---------------------------------------------------------------------------


void __fastcall Thauptformular::FormCreate(TObject *Sender)
{

  TRegistry * myRegistry;
  myRegistry = new TRegistry;
  if (myRegistry->OpenKey("Software\\BayesX", false))
    {
    if (myRegistry->ValueExists("hauptformularleft"))
      Left = myRegistry->ReadInteger("hauptformularleft");
    else
      Left = 0;
    if (myRegistry->ValueExists("hauptformulartop"))
      Top = myRegistry->ReadInteger("hauptformulartop");
    else
      Top = 0;
    if (myRegistry->ValueExists("hauptformularheight"))
      Height = myRegistry->ReadInteger("hauptformularheight");
    else
      Height = 90;
    if (myRegistry->ValueExists("hauptformularwidth"))
      Width = myRegistry->ReadInteger("hauptformularwidth");
    else
      Width = Frame->ClientWidth - 4;
    }
  else
    {
    Left = 0;
    Top = 0;
    Height = 90;
    Width = Frame->ClientWidth - 4;
    }
  myRegistry->CloseKey();
  delete myRegistry;

  commandedit->Width = hauptformular->ClientWidth;
  commandedit->Clear();

  for(unsigned i=0;i<hauptformular->get_objecttype().size();i++)
    Objectbrowser->ObjectbrowserListBoxType->Items->Add(hauptformular->
      get_objecttype()[i].strtochar());

}
//---------------------------------------------------------------------------


int Thauptformular::CheckOutputSaved(void)
{

  int button;
  if (Results->ResultsRichEdit->Modified == true)
    {
    button = Application->MessageBox("Output has been changed!\nSave changes?", "Statwin Output",
      MB_YESNOCANCEL + MB_ICONEXCLAMATION + MB_DEFBUTTON1);
    if(button == IDYES)
      {
      if (Results->ResultsRichEdit->Tag == true)
        {
        Results->ResultsRichEdit->Lines->SaveToFile(Frame->SaveDialog->FileName);
        Results->ResultsRichEdit->Modified = false;
        }
      else if (Frame->SaveDialog->Execute())
        {
        Results->ResultsRichEdit->Lines->SaveToFile(Frame->SaveDialog->FileName);
        Results->ResultsRichEdit->Modified = false;
        Results->ResultsRichEdit->Tag = true;
        }
      }
    }

  return button;

}


bool Thauptformular::breakcommand(void)
  {
  Application->ProcessMessages();
  if (Frame->stop)
    {
    out("USER BREAK\n");
    return true;
    }


  if (Frame->pause)
    {
    out("PROGRAM PAUSED\n");
    out("Click CONTINUE to proceed\n");

    while (Frame->pause)
      {
      Application->ProcessMessages();
      }

    out("CONTINUED\n");
    out("\n");

    }

  return false;
  }
