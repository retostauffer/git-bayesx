//---------------------------------------------------------------------------
#ifndef describe_datasetH
#define describe_datasetH
//---------------------------------------------------------------------------
#include <Classes.hpp>
#include <Controls.hpp>
#include <StdCtrls.hpp>
#include <Forms.hpp>
#include <DBGrids.hpp>
#include <Grids.hpp>
#include <data.h>
//---------------------------------------------------------------------------
class Tdatasetform : public TForm
{
__published:	// Komponenten, die von der IDE verwaltet werden
    TStringGrid *StringGrid;
    TScrollBar *VertScrollBar;
    TScrollBar *HorzScrollBar;
    void __fastcall FormActivate(TObject *Sender);
    
    
    
    void __fastcall VertScrollBarChange(TObject *Sender);
    void __fastcall FormResize(TObject *Sender);
    
    void __fastcall FormClose(TObject *Sender, TCloseAction &Action);
    
        void __fastcall HorzScrollBarChange(TObject *Sender);
    
private:	// Benutzerdeklarationen
public:		// Benutzerdeklarationen
dataset * datap;
ST::string dataname;
    __fastcall Tdatasetform(TComponent* Owner);
void fillGrid(void);    
};
//---------------------------------------------------------------------------
extern PACKAGE Tdatasetform *datasetform;
//---------------------------------------------------------------------------
#endif
