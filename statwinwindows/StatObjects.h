//---------------------------------------------------------------------------
#ifndef StatObjectsH
#define StatObjectsH
//---------------------------------------------------------------------------
#include <Classes.hpp>
#include <Controls.hpp>
#include <StdCtrls.hpp>
#include <Forms.hpp>
#include <ExtCtrls.hpp>
//---------------------------------------------------------------------------
class TObjectbrowser : public TForm
{
__published:	// Komponenten, die von der IDE verwaltet werden
    TListBox *ObjectbrowserListBoxType;
    TListBox *ObjectbrowserListBoxObjects;
    TSplitter *Splitter;
    void __fastcall FormCreate(TObject *Sender);
    void __fastcall ObjectbrowserListBoxTypeClick(TObject *Sender);
    void __fastcall FormClick(TObject *Sender);
    void __fastcall ObjectbrowserListBoxObjectsKeyDown(TObject *Sender,
          WORD &Key, TShiftState Shift);
    void __fastcall ObjectbrowserListBoxObjectsDblClick(TObject *Sender);
    
private:	// Benutzerdeklarationen
public:		// Benutzerdeklarationen
    __fastcall TObjectbrowser(TComponent* Owner);
};
//---------------------------------------------------------------------------
extern PACKAGE TObjectbrowser *Objectbrowser;
//---------------------------------------------------------------------------
#endif
