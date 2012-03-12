//---------------------------------------------------------------------------
#include <vcl.h>
#pragma hdrstop

#include <Registry.hpp>
#include <Printers.hpp>

#include "StatwinFrame.h"
#include "statwin_haupt.h"
#include "StatReview.h"
#include "StatResults.h"
#include "StatObjects.h"
//#include "sicher.h"
//#include "describe_map.h"

//---------------------------------------------------------------------------
#pragma package(smart_init)
#pragma resource "*.dfm"
TFrame *Frame;
//---------------------------------------------------------------------------
__fastcall TFrame::TFrame(TComponent* Owner)
    : TForm(Owner)
{
stop = false;
pause = false;
suppoutput=false;
processruning = false;
}
//---------------------------------------------------------------------------


void __fastcall TFrame::FrameMenuWindowClick(TObject *Sender)
{
  FMWCommand->Checked = hauptformular->Active;
  FMWReview->Checked = Review->Active;
  FMWOutput->Checked = Results->Active;
  FMWObjectbrowser->Checked = Objectbrowser->Active;
}
//---------------------------------------------------------------------------

void __fastcall TFrame::FMWCommandClick(TObject *Sender)
{
  hauptformular->WindowState = wsNormal;
  hauptformular->Show();
}
//---------------------------------------------------------------------------

void __fastcall TFrame::FMWOutputClick(TObject *Sender)
{
  Results->WindowState = wsNormal;
  Results->Show();
}
//---------------------------------------------------------------------------

void __fastcall TFrame::FMWReviewClick(TObject *Sender)
{
  Review->WindowState = wsNormal;
  Review->Show();
}
//---------------------------------------------------------------------------

void __fastcall TFrame::FMWObjectbrowserClick(TObject *Sender)
{

  Objectbrowser->ObjectbrowserListBoxType->Clear();
  Objectbrowser->WindowState = wsNormal;
  Objectbrowser->Show();
  for(unsigned i=0;i<hauptformular->get_objecttype().size();i++)
    Objectbrowser->ObjectbrowserListBoxType->Items->Add(hauptformular->
      get_objecttype()[i].strtochar());

}
//---------------------------------------------------------------------------

void __fastcall TFrame::FormCloseQuery(TObject *Sender, bool &CanClose)
{

  int button;
  if(processruning)
    {
//    Application->MessageBox("You can not quit the program when a process is running!",
//             "Exit query", MB_OK + MB_ICONSTOP);
//    CanClose = false;

    int answer;
    answer = Application->MessageBox(
        "Are you sure? Do you want\nto interrupt the current process\nand exit the program?",
        "Break", MB_YESNO + MB_ICONQUESTION + MB_DEFBUTTON2);
    if (answer == IDYES)
      {
      Frame->stop = true;
      CanClose = true;
      }
    else
      CanClose = false;

    }
  else
    {
    button = hauptformular->CheckOutputSaved();
    if(button != IDCANCEL)
      CanClose = true;
    else
      CanClose = false;
    }
}
//---------------------------------------------------------------------------

void __fastcall TFrame::FileExitClick(TObject *Sender)
{
// sicherform->ShowModal();
  if(processruning)
    {
//    Application->MessageBox("You can not quit the program when a process is running!",
//             "Exit query", MB_OK + MB_ICONSTOP);
    int answer;
    answer = Application->MessageBox(
        "Are you sure? Do you want\nto interrupt the current process\nand exit the program?",
        "Break", MB_YESNO + MB_ICONQUESTION + MB_DEFBUTTON2);
    if (answer == IDYES)
      Frame->stop = true;
    }
  else
    {
    writetoregistry();
    if(hauptformular->CheckOutputSaved() != IDCANCEL)
      {
      Application->Terminate();
      }
    }
}
//---------------------------------------------------------------------------


void __fastcall TFrame::FormClose(TObject *Sender, TCloseAction &Action)
{
  writetoregistry();
}
//---------------------------------------------------------------------------

void __fastcall TFrame::FormCreate(TObject *Sender)
{

  Application->OnDeactivate = AppDeactivate;
  TRegistry * myRegistry;
  myRegistry = new TRegistry;
  if (myRegistry->OpenKey("Software\\BayesX", false))
    {
    if (myRegistry->ValueExists("left"))
      Left = myRegistry->ReadInteger("left");
    if (myRegistry->ValueExists("top"))
      Top = myRegistry->ReadInteger("top");
    if (myRegistry->ValueExists("height"))
      Height = myRegistry->ReadInteger("height");
    if (myRegistry->ValueExists("width"))
      Width = myRegistry->ReadInteger("width");
    }
  myRegistry->CloseKey();
  delete myRegistry;

}
//---------------------------------------------------------------------------


void __fastcall TFrame::breakbuttonClick(TObject *Sender)
{

int answer;
answer = Application->MessageBox("Are you sure? Do you want\nto interrupt the current process?",
             "Break", MB_YESNO + MB_ICONQUESTION + MB_DEFBUTTON2);
if (answer == IDYES)
  Frame->stop = true;

//  sicherform->ShowModal();
}
//---------------------------------------------------------------------------

void __fastcall TFrame::pausebuttonClick(TObject *Sender)
{

if (pause==true)
  {
  pause = false;
  pausebutton->Caption = "PAUSE";
  }
else
  {
  pause = true;
  pausebutton->Caption = "CONTINUE";
  }


}
//---------------------------------------------------------------------------


void TFrame::writetoregistry(void)
  {

  TRegistry * myRegistry;
  myRegistry = new TRegistry;
  if (myRegistry->OpenKey("Software\\BayesX", true))
    {
    if(WindowState == wsNormal)
      {
      myRegistry->WriteInteger("left",Left);
      myRegistry->WriteInteger("top",Top);
      myRegistry->WriteInteger("height",Height);
      myRegistry->WriteInteger("width",Width);
      }
    if(Review->WindowState == wsNormal)
      {
      myRegistry->WriteInteger("Reviewleft",Review->Left);
      myRegistry->WriteInteger("Reviewtop",Review->Top);
      myRegistry->WriteInteger("Reviewheight",Review->Height);
      myRegistry->WriteInteger("Reviewwidth",Review->Width);
      }
    if(Objectbrowser->WindowState == wsNormal)
      {
      myRegistry->WriteInteger("Objectbrowserleft",Objectbrowser->Left);
      myRegistry->WriteInteger("Objectbrowsertop",Objectbrowser->Top);
      myRegistry->WriteInteger("Objectbrowserheight",Objectbrowser->Height);
      myRegistry->WriteInteger("Objectbrowserwidth",Objectbrowser->Width);
      }
    if(Results->WindowState == wsNormal)
      {
      myRegistry->WriteInteger("Resultsleft",Results->Left);
      myRegistry->WriteInteger("Resultstop",Results->Top);
      myRegistry->WriteInteger("Resultsheight",Results->Height);
      myRegistry->WriteInteger("Resultswidth",Results->Width);
      }
    if(hauptformular->WindowState == wsNormal)
      {
      myRegistry->WriteInteger("hauptformularleft",hauptformular->Left);
      myRegistry->WriteInteger("hauptformulartop",hauptformular->Top);
      myRegistry->WriteInteger("hauptformularheight",hauptformular->Height);
      myRegistry->WriteInteger("hauptformularwidth",hauptformular->Width);
      }
    }

//  myRegistry->WriteInteger("Splitter",Objectbrowser->ObjectbrowserListBoxType->Width);

  myRegistry->CloseKey();
  delete myRegistry;

  }
void __fastcall TFrame::FormDestroy(TObject *Sender)
{
//  writetoregistry();
}
//---------------------------------------------------------------------------



void __fastcall TFrame::FileSaveClick(TObject *Sender)
{
  if (Results->ResultsRichEdit->Tag == true){
    if (Results->ResultsRichEdit->Modified == true){
      Results->ResultsRichEdit->Lines->SaveToFile(SaveDialog->FileName);
      Results->ResultsRichEdit->Modified = false;
    }
  }
  else
    if (SaveDialog->Execute())
      {

      if (SaveDialog->FilterIndex == 1)
        Results->ResultsRichEdit->PlainText = false;
      else
        Results->ResultsRichEdit->PlainText = true;

      Results->ResultsRichEdit->Lines->SaveToFile(SaveDialog->FileName);
      Results->ResultsRichEdit->Modified = false;
      Results->ResultsRichEdit->Tag = true;
      Results->Caption = SaveDialog->FileName;
    }
}
//---------------------------------------------------------------------------


void __fastcall TFrame::FileSaveasClick(TObject *Sender)
{
  if (SaveDialog->Execute())
    {

    if (SaveDialog->FilterIndex == 1)
      Results->ResultsRichEdit->PlainText = false;
    else
      Results->ResultsRichEdit->PlainText = true;

    Results->ResultsRichEdit->Lines->SaveToFile(SaveDialog->FileName);
    Results->ResultsRichEdit->Modified = false;
    Results->ResultsRichEdit->Tag = true;
    Results->Caption = SaveDialog->FileName;
  }
}
//---------------------------------------------------------------------------

void __fastcall TFrame::FileOpenClick(TObject *Sender)
{
  if (OpenDialog->Execute()){
    Results->ResultsRichEdit->Lines->LoadFromFile(OpenDialog->FileName);
    Results->ResultsRichEdit->Modified = false;
    Results->ResultsRichEdit->ReadOnly = OpenDialog->Options.Contains(ofReadOnly);
    Results->ResultsRichEdit->Lines->Append("");
    Results->ResultsRichEdit->Tag = true;
  }
}
//---------------------------------------------------------------------------

void __fastcall TFrame::FilePrintClick(TObject *Sender)
{
  if (PrintDialog->Execute())
    Results->ResultsRichEdit->Print(PrintDialog->Name);
}
//---------------------------------------------------------------------------

void __fastcall TFrame::EditBoldClick(TObject *Sender)
{
  if (EditBold->Checked == true)
    Results->ResultsRichEdit->SelAttributes->Style =
      Results->ResultsRichEdit->SelAttributes->Style >> fsBold;
  else
    Results->ResultsRichEdit->SelAttributes->Style =
      Results->ResultsRichEdit->SelAttributes->Style << fsBold;
}
//---------------------------------------------------------------------------

void __fastcall TFrame::EditItalicClick(TObject *Sender)
{
  if (EditItalic->Checked == true)
    Results->ResultsRichEdit->SelAttributes->Style =
      Results->ResultsRichEdit->SelAttributes->Style >> fsItalic;
  else
    Results->ResultsRichEdit->SelAttributes->Style =
      Results->ResultsRichEdit->SelAttributes->Style << fsItalic;
}
//---------------------------------------------------------------------------

void __fastcall TFrame::EditUnderlineClick(TObject *Sender)
{
  if (EditUnderline->Checked == true)
    Results->ResultsRichEdit->SelAttributes->Style =
      Results->ResultsRichEdit->SelAttributes->Style >> fsUnderline;
  else
    Results->ResultsRichEdit->SelAttributes->Style =
      Results->ResultsRichEdit->SelAttributes->Style << fsUnderline;

}
//---------------------------------------------------------------------------

void __fastcall TFrame::EditClick(TObject *Sender)
{
  EditBold->Checked = Results->ResultsRichEdit->SelAttributes->Style.Contains(fsBold);
  EditItalic->Checked = Results->ResultsRichEdit->SelAttributes->Style.Contains(fsItalic);
  EditUnderline->Checked = Results->ResultsRichEdit->SelAttributes->Style.Contains(fsUnderline);
}
//---------------------------------------------------------------------------

void __fastcall TFrame::ClearOutputWindowClick(TObject *Sender)
{
  Results->ResultsRichEdit->Clear();
}
//---------------------------------------------------------------------------
void __fastcall TFrame::AppDeactivate(TObject *Sender)
{
//  Application->Minimize();
}
void __fastcall TFrame::supproutputClick(TObject *Sender)
{

if (suppoutput==true)
  {
  suppoutput = false;
  supproutput->Caption = "SUPPRESS OUTPUT";
  }
else
  {
  suppoutput = true;
  supproutput->Caption = "SHOW OUTPUT";
  }

}
//---------------------------------------------------------------------------

