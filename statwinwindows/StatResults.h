//---------------------------------------------------------------------------
#ifndef StatResultsH
#define StatResultsH
//---------------------------------------------------------------------------
#include <Classes.hpp>
#include <Controls.hpp>
#include <StdCtrls.hpp>
#include <Forms.hpp>
#include <ComCtrls.hpp>
#include <Menus.hpp>
#include <Dialogs.hpp>
//---------------------------------------------------------------------------
class TResults : public TForm
{
__published:	// Komponenten, die von der IDE verwaltet werden
    TRichEdit *ResultsRichEdit;
    void __fastcall FormCreate(TObject *Sender);
    
    
    
    
    
private:	// Benutzerdeklarationen
public:		// Benutzerdeklarationen
bool stop;
    __fastcall TResults(TComponent* Owner);
};
//---------------------------------------------------------------------------
extern PACKAGE TResults *Results;
//---------------------------------------------------------------------------
#endif
