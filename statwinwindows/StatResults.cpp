//---------------------------------------------------------------------------
#include <vcl.h>
#pragma hdrstop

#include <Registry.hpp>

#include "StatResults.h"
#include "statwin_haupt.h"
#include "StatwinFrame.h"
#include "StatReview.h"
#include "StatObjects.h"
//#include "sicher.h"
//---------------------------------------------------------------------------
#pragma package(smart_init)
#pragma resource "*.dfm"
TResults *Results;
//---------------------------------------------------------------------------
__fastcall TResults::TResults(TComponent* Owner)
    : TForm(Owner)
{

}
//---------------------------------------------------------------------------
void __fastcall TResults::FormCreate(TObject *Sender)
{

  TRegistry * myRegistry;
  myRegistry = new TRegistry;
  if (myRegistry->OpenKey("Software\\BayesX", false))
    {
    if (myRegistry->ValueExists("Resultsleft"))
      Left = myRegistry->ReadInteger("Resultsleft");
    else
      Left = 0;
    if (myRegistry->ValueExists("Resultstop"))
      Top = myRegistry->ReadInteger("Resultstop");
    else
      Top = 90;
    if (myRegistry->ValueExists("Resultsheight"))
      Height = myRegistry->ReadInteger("Resultsheight");
    else
      Height = Frame->ClientHeight - 90 - 4;
    if (myRegistry->ValueExists("Resultswidth"))
      Width = myRegistry->ReadInteger("Resultswidth");
    else
      Width = (Frame->ClientWidth - 4) * 3 / 4;
    }
  else
    {
    Left = 0;
    Top = 90;
    Height = Frame->ClientHeight - 90 - 4;
    Width = (Frame->ClientWidth - 4) * 3 / 4;
    }
  myRegistry->CloseKey();
  delete myRegistry;

  ResultsRichEdit->Clear();
  ResultsRichEdit->Tag = false;

}
//---------------------------------------------------------------------------



















