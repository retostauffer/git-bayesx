//---------------------------------------------------------------------------
#ifndef StatwinFrameH
#define StatwinFrameH
//---------------------------------------------------------------------------
#include <Classes.hpp>
#include <Controls.hpp>
#include <StdCtrls.hpp>
#include <Forms.hpp>
#include <Menus.hpp>
#include <Dialogs.hpp>
#include <ExtCtrls.hpp>
//---------------------------------------------------------------------------
class TFrame : public TForm
{
__published:	// Komponenten, die von der IDE verwaltet werden
    TMainMenu *FrameMenu;
    TMenuItem *FrameMenuWindow;
    TMenuItem *FrameMenuHelp;
    TMenuItem *FMWCommand;
    TMenuItem *FMWOutput;
    TMenuItem *FMWReview;
    TMenuItem *FMWObjectbrowser;
    TMenuItem *File;
    TMenuItem *FileExit;
    TMenuItem *Edit;
    TMenuItem *EditBold;
    TMenuItem *EditItalic;
    TMenuItem *EditUnderline;
    TMenuItem *FileOpen;
    TMenuItem *FileSave;
    TMenuItem *FileSaveas;
    TMenuItem *FilePrint;
    TBevel *Bevel1;
    TPanel *Panel1;
    TButton *breakbutton;
    TButton *pausebutton;
    TOpenDialog *OpenDialog;
    TSaveDialog *SaveDialog;
    TPrintDialog *PrintDialog;
    TMenuItem *ClearOutputWindow;
    TButton *supproutput;
    void __fastcall FrameMenuWindowClick(TObject *Sender);
    void __fastcall FMWCommandClick(TObject *Sender);
    void __fastcall FMWOutputClick(TObject *Sender);
    void __fastcall FMWReviewClick(TObject *Sender);
    void __fastcall FMWObjectbrowserClick(TObject *Sender);
    void __fastcall FormCloseQuery(TObject *Sender, bool &CanClose);
    void __fastcall FileExitClick(TObject *Sender);
    
    
    void __fastcall FormClose(TObject *Sender, TCloseAction &Action);
    void __fastcall FormCreate(TObject *Sender);
    void __fastcall breakbuttonClick(TObject *Sender);
    void __fastcall pausebuttonClick(TObject *Sender);
    void __fastcall FormDestroy(TObject *Sender);
    void __fastcall FileSaveClick(TObject *Sender);
    
    void __fastcall FileSaveasClick(TObject *Sender);
    void __fastcall FileOpenClick(TObject *Sender);
    void __fastcall FilePrintClick(TObject *Sender);
    void __fastcall EditBoldClick(TObject *Sender);
    void __fastcall EditItalicClick(TObject *Sender);
    void __fastcall EditUnderlineClick(TObject *Sender);
    void __fastcall EditClick(TObject *Sender);
    
    void __fastcall ClearOutputWindowClick(TObject *Sender);
    void __fastcall supproutputClick(TObject *Sender);
private:	// Benutzerdeklarationen
public:		// Benutzerdeklarationen
bool stop;
bool pause;
bool suppoutput;
bool processruning;
void writetoregistry(void);
    __fastcall TFrame(TComponent* Owner);
    void __fastcall AppDeactivate(TObject *Sender);
};
//---------------------------------------------------------------------------
extern PACKAGE TFrame *Frame;
//---------------------------------------------------------------------------
#endif
