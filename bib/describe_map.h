//---------------------------------------------------------------------------
#ifndef describe_mapH
#define describe_mapH
//---------------------------------------------------------------------------
#include <Classes.hpp>
#include <Controls.hpp>
#include <StdCtrls.hpp>
#include <Forms.hpp>
#include <ExtCtrls.hpp>
#include <map.h>
//---------------------------------------------------------------------------
class Tmapform : public TForm
{
__published:	// Komponenten, die von der IDE verwaltet werdenvoid __fastcall FormShow(TObject *Sender);void __fastcall FormResize(TObject *Sender);void __fastcall FormResize(TObject *Sender);
    void __fastcall FormPaint(TObject *Sender);
    
    void __fastcall FormResize(TObject *Sender);
private:	// Benutzerdeklarationen
public:		// Benutzerdeklarationen
MAP::map * mapinfo;
ST::string mapname;
    __fastcall Tmapform(TComponent* Owner);
void drawmap(void);    
};
//---------------------------------------------------------------------------
extern PACKAGE Tmapform *mapform;
//---------------------------------------------------------------------------
#endif
