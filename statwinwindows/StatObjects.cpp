//---------------------------------------------------------------------------
#include <vcl.h>
#pragma hdrstop

#include <Registry.hpp>

#include "StatObjects.h"
#include "StatwinFrame.h"
#include "statwin_haupt.h"
#include "StatReview.h"
//---------------------------------------------------------------------------
#pragma package(smart_init)
#pragma resource "*.dfm"
TObjectbrowser *Objectbrowser;
//---------------------------------------------------------------------------
__fastcall TObjectbrowser::TObjectbrowser(TComponent* Owner)
    : TForm(Owner)
{
}
//---------------------------------------------------------------------------
void __fastcall TObjectbrowser::FormCreate(TObject *Sender)
{

  TRegistry * myRegistry;
  myRegistry = new TRegistry;

  if (myRegistry->OpenKey("Software\\BayesX", false))
    {
    if (myRegistry->ValueExists("Splitter"))
      {
      Splitter->Left = myRegistry->ReadInteger("Splitter");
      ObjectbrowserListBoxType->Width = myRegistry->ReadInteger("Splitter");
      }
    if (myRegistry->ValueExists("Objectbrowserleft"))
      Left = myRegistry->ReadInteger("Objectbrowserleft");
    else
      Left = (Frame->ClientWidth - 4) / 2;
    if (myRegistry->ValueExists("Objectbrowsertop"))
      Top = myRegistry->ReadInteger("Objectbrowsertop");
    else
      Top = 90 + (Frame->ClientHeight - 90 - 4) / 2;
    if (myRegistry->ValueExists("Objectbrowserheight"))
      Height = myRegistry->ReadInteger("Objectbrowserheight");
    else
      Height = (Frame->ClientHeight - 90 - 4) / 2;
    if (myRegistry->ValueExists("Objectbrowserwidth"))
      Width = myRegistry->ReadInteger("Objectbrowserwidth");
    else
      Width = (Frame->ClientWidth - 4) / 2;
    }
  else
    {
    Left = (Frame->ClientWidth - 4) / 2;
    Top = 90 + (Frame->ClientHeight - 90 - 4) / 2;
    Height = (Frame->ClientHeight - 90 - 4) / 2;
    Width = (Frame->ClientWidth - 4) / 2;
    }

  myRegistry->CloseKey();
  delete myRegistry;

}
//---------------------------------------------------------------------------

void __fastcall TObjectbrowser::ObjectbrowserListBoxTypeClick(
      TObject *Sender)
{
  ObjectbrowserListBoxObjects->Clear();
  unsigned i;
  for(i=0;i<hauptformular->get_objects().size();i++)
    {
    if(strcmp(ObjectbrowserListBoxType->Items->
         Strings[ObjectbrowserListBoxType->ItemIndex].c_str(),
         hauptformular->get_objects()[i]->gettype().strtochar()) == 0)
      {
      ObjectbrowserListBoxObjects->Items->Add(hauptformular->
        get_objects()[i]->getname().strtochar());
      }
    }

}
//---------------------------------------------------------------------------


void __fastcall TObjectbrowser::FormClick(TObject *Sender)
{

  for(unsigned i=0;i<hauptformular->get_objecttype().size();i++)
    Objectbrowser->ObjectbrowserListBoxType->Items->Add(hauptformular->
      get_objecttype()[i].strtochar());

}
//---------------------------------------------------------------------------

void __fastcall TObjectbrowser::ObjectbrowserListBoxObjectsKeyDown(
      TObject *Sender, WORD &Key, TShiftState Shift)
{
  if(Key == VK_DELETE)
    {
    if(ObjectbrowserListBoxObjects->ItemIndex >= 0)
      hauptformular->parsecommand("drop "+ST::string(ObjectbrowserListBoxObjects->Items->
         Strings[ObjectbrowserListBoxObjects->ItemIndex].c_str()));
// display remaining objects in the objectbrowser window
    ObjectbrowserListBoxObjects->Clear();
    unsigned i;
    for(i=0;i<hauptformular->get_objects().size();i++)
      {
      if(strcmp(ObjectbrowserListBoxType->Items->
         Strings[ObjectbrowserListBoxType->ItemIndex].c_str(),
         hauptformular->get_objects()[i]->gettype().strtochar()) == 0)
        {
        ObjectbrowserListBoxObjects->Items->Add(hauptformular->
          get_objects()[i]->getname().strtochar());
        }
      }
    }
}
//---------------------------------------------------------------------------

void __fastcall TObjectbrowser::ObjectbrowserListBoxObjectsDblClick(
      TObject *Sender)
{
  hauptformular->parsecommand(ST::string(ObjectbrowserListBoxObjects->Items->
         Strings[ObjectbrowserListBoxObjects->ItemIndex].c_str()) + ".describe");
}
//---------------------------------------------------------------------------


