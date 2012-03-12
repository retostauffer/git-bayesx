//---------------------------------------------------------------------------
#ifndef StatReviewH
#define StatReviewH
//---------------------------------------------------------------------------
#include <Classes.hpp>
#include <Controls.hpp>
#include <StdCtrls.hpp>
#include <Forms.hpp>
//---------------------------------------------------------------------------
class TReview : public TForm
{
__published:	// Komponenten, die von der IDE verwaltet werden
    TListBox *ReviewListBox;
    void __fastcall FormCreate(TObject *Sender);
    void __fastcall ReviewListBoxClick(TObject *Sender);
    void __fastcall ReviewListBoxDblClick(TObject *Sender);
    
private:	// Benutzerdeklarationen
public:		// Benutzerdeklarationen
    __fastcall TReview(TComponent* Owner);
};
//---------------------------------------------------------------------------
extern PACKAGE TReview *Review;
//---------------------------------------------------------------------------
#endif
