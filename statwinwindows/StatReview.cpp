//---------------------------------------------------------------------------
#include <vcl.h>
#pragma hdrstop

#include <Registry.hpp>

#include "StatReview.h"
#include "statwin_haupt.h"
#include "StatwinFrame.h"
//---------------------------------------------------------------------------
#pragma package(smart_init)
#pragma resource "*.dfm"
TReview *Review;
//---------------------------------------------------------------------------
__fastcall TReview::TReview(TComponent* Owner)
    : TForm(Owner)
{
}
//---------------------------------------------------------------------------
void __fastcall TReview::FormCreate(TObject *Sender)
{
  TRegistry * myRegistry;
  myRegistry = new TRegistry;
  if (myRegistry->OpenKey("Software\\BayesX", false))
    {
/*    if(myRegistry->ValueExists("ReviewWS"))
      {
      if(myRegistry->ReadInteger("ReviewWS") == 0)
        {
        WindowState = wsNormal;
*/
        if (myRegistry->ValueExists("Reviewleft"))
          Left = myRegistry->ReadInteger("Reviewleft");
        else
          Left = (Frame->ClientWidth - 4) * 3 / 4;
        if (myRegistry->ValueExists("Reviewtop"))
          Top = myRegistry->ReadInteger("Reviewtop");
        else
          Top = 90;
        if (myRegistry->ValueExists("Reviewheight"))
          Height = myRegistry->ReadInteger("Reviewheight");
        else
          Height = (Frame->ClientHeight - 90 - 4);
        if (myRegistry->ValueExists("Reviewwidth"))
          Width = myRegistry->ReadInteger("Reviewwidth");
        else
          Width = (Frame->ClientWidth - 4) / 4;
/*        }
      else if(myRegistry->ReadInteger("ReviewWS") == 1)
        WindowState = wsMinimized;
      else if(myRegistry->ReadInteger("ReviewWS") == 2)
        WindowState = wsMaximized;
      }
*/
    }
  else
    {
    Left = (Frame->ClientWidth - 4) * 3 / 4;
    Top = 90;
    Height = (Frame->ClientHeight - 90 - 4);
    Width = (Frame->ClientWidth - 4) / 4;
    }
  myRegistry->CloseKey();
  delete myRegistry;

}
//---------------------------------------------------------------------------


void __fastcall TReview::ReviewListBoxClick(TObject *Sender)
{
  hauptformular->commandedit->Clear();
  hauptformular->commandedit->Lines->Append(AnsiString(ReviewListBox->Items->Strings[ReviewListBox->ItemIndex]));
  hauptformular->Show();
}
//---------------------------------------------------------------------------

void __fastcall TReview::ReviewListBoxDblClick(TObject *Sender)
{
  hauptformular->commandedit->Clear();
  hauptformular->commandedit->Lines->Append(AnsiString(ReviewListBox->Items->Strings[ReviewListBox->ItemIndex]));
  hauptformular->Show();
}
//---------------------------------------------------------------------------


