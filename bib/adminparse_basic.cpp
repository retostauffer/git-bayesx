
#include "adminparse_basic.h"

administrator_basic::administrator_basic(void)
  {

  pause = false;
  stop = false;
  processrunning = false;
  suppressoutput = false;

  }
//---------------------------------------------------------------------------

bool administrator_basic::breakcommand(void)
  {

  if(stop)
    return true;

  if(pause)
    {

    while(pause)
      {
      if(stop)
        {
        return true;
        }
      }

    }

  return false;
  }







