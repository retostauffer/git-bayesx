
#include "first.h"

#include "statmat_penalty.h"

namespace STATMAT_PENALTY
{

statmatrix<double> Kmrf(const MAP::map & m)
  {
  unsigned r,i,j,n;
  r = m.get_nrregions();

  statmatrix<double>res(r,r,0);

  for(i=0; i<r; i++)
    {
    res(i,i) = m.get_weightssum(i);
    for(j=0; j<m.get_neighbors()[i].size(); j++)
      {
      n=m.get_neighbors()[i][j];
      res(i,n) = -m.get_weights()[i][j];
      res(n,i) = res(i,n);
      }
    }

  return res;
  }

statmatrix<double> K2dim_pspline(const unsigned & nknots)
  {
  statmatrix<double> res(nknots*nknots,nknots*nknots,0);
  unsigned i,j;

// 4 corners:
  res(0,0)=2;
  res(0,1)=-1;
  res(0,nknots)=-1;

  res(nknots-1,nknots-1)=2;
  res(nknots-1,nknots-2)=-1;
  res(nknots-1,2*nknots-1)=-1;

  res((nknots-1)*nknots,(nknots-1)*nknots)=2;
  res((nknots-1)*nknots,(nknots-2)*nknots)=-1;
  res((nknots-1)*nknots,(nknots-1)*nknots+1)=-1;

  res(nknots*nknots-1,nknots*nknots-1)=2;
  res(nknots*nknots-1,(nknots-1)*nknots-1)=-1;
  res(nknots*nknots-1,nknots*nknots-2)=-1;

// upper and lower edge

  for(i=1; i<nknots-1; i++)
    {
    res(i,i)=3;
    res(i,i-1)=-1;
    res(i,i+1)=-1;
    res(i,i+nknots)=-1;
    }

  for(i=(nknots-1)*nknots+1; i<nknots*nknots-1; i++)
    {
    res(i,i)=3;
    res(i,i-1)=-1;
    res(i,i+1)=-1;
    res(i,i-nknots)=-1;
    }

// rest

  for(i=2; i<nknots; i++)
    {
    res((i-1)*nknots,(i-1)*nknots)=3;
    res((i-1)*nknots,(i-2)*nknots)=-1;
    res((i-1)*nknots,(i-1)*nknots+1)=-1;
    res((i-1)*nknots,i*nknots)=-1;

    for(j=1; j<nknots-1; j++)
      {
      res((i-1)*nknots+j,(i-1)*nknots+j)=4;
      res((i-1)*nknots+j,(i-2)*nknots+j)=-1;
      res((i-1)*nknots+j,(i-1)*nknots+j-1)=-1;
      res((i-1)*nknots+j,(i-1)*nknots+j+1)=-1;
      res((i-1)*nknots+j,i*nknots+j)=-1;
      }

    res(i*nknots-1,i*nknots-1)=3;
    res(i*nknots-1,(i-1)*nknots-1)=-1;
    res(i*nknots-1,i*nknots-2)=-1;
    res(i*nknots-1,(i+1)*nknots-1)=-1;
    }
  return res;
  }

//}


statmatrix<double> K2dim_pspline_rw2(const unsigned & nknots, const unsigned & ox, const unsigned & oy)
  {
  statmatrix<double> res(nknots*nknots,nknots*nknots,0);
//  unsigned i,j;

  statmatrix<double>  Dx = diffmat(ox,nknots);       // difference matrix of order ox
  statmatrix<double> DDx = Dx.transposed() * Dx;     // D'D

  statmatrix<double>  Dy = diffmat(oy,nknots);       // difference matrix of order oy
  statmatrix<double> DDy = Dy.transposed() * Dy;

  statmatrix<double> I = statmatrix<double>::diag(nknots,1);  // identity matrix

  statmatrix<double> Px = kronecker(I,DDx);      // penalty matrix of x-direction
  statmatrix<double> Py = kronecker(DDy,I);      // penalty matrix of y-direction
  res = Px + Py;
  return res;
  }

statmatrix<double> K2dim_pspline_biharmonic(const unsigned & nknots)
  {
  unsigned i,j;
  statmatrix<double> res(nknots*nknots,nknots*nknots,0);

// 4 corners

// upper left
res(0,0) = 4;
res(0,1) = res(0,nknots) = -4;
res(0,2) = res(0,2*nknots) = 1;
res(0,nknots+1) = 2;

// upper right
res(nknots-1,nknots-1) = 4;
res(nknots-1,nknots-2) = res(nknots-1,2*nknots-1) = -4;
res(nknots-1,nknots-3) = res(nknots-1,3*nknots-1) = 1;
res(nknots-1,2*nknots-2) = 2;

// lower left
res(nknots*(nknots-1),nknots*(nknots-1)) = 4;
res(nknots*(nknots-1),nknots*(nknots-1)+1) = res(nknots*(nknots-1),nknots*(nknots-2)) = -4;
res(nknots*(nknots-1),nknots*(nknots-1)+2) = res(nknots*(nknots-1),nknots*(nknots-3)) = 1;
res(nknots*(nknots-1),nknots*(nknots-2)+1) = 2;

// lower right
res(nknots*nknots-1,nknots*nknots-1) = 4;
res(nknots*nknots-1,nknots*nknots-2) = res(nknots*nknots-1,nknots*(nknots-1)-1) = -4;
res(nknots*nknots-1,nknots*nknots-3) = res(nknots*nknots-1,nknots*(nknots-2)-1) = 1;
res(nknots*nknots-1,nknots*(nknots-1)-2) = 2;

// 4 second order corners

// upper left
res(nknots+1,nknots+1) = 18;
res(nknots+1,nknots+2) = res(nknots+1,2*nknots+1) = -8;
res(nknots+1,nknots) = res(nknots+1,1) = -6;
res(nknots+1,0) = res(nknots+1,2) = res(nknots+1,2*nknots) = res(nknots+1,2*nknots+2) = 2;
res(nknots+1,nknots+3) = res(nknots+1,3*nknots+1) = 1;

// upper right
res(2*nknots-2,2*nknots-2) = 18;
res(2*nknots-2,2*nknots-3) = res(2*nknots-2,3*nknots-2) = -8;
res(2*nknots-2,2*nknots-1) = res(2*nknots-2,nknots-2) = -6;
res(2*nknots-2,nknots-3) = res(2*nknots-2,nknots-1) = res(2*nknots-2,3*nknots-3) = res(2*nknots-2,3*nknots-1) = 2;
res(2*nknots-2,4*nknots-2) = res(2*nknots-2,2*nknots-4) = 1;

// lower left
res(nknots*(nknots-2)+1,nknots*(nknots-2)+1) = 18;
res(nknots*(nknots-2)+1,nknots*(nknots-2)+2) = res(nknots*(nknots-2)+1,nknots*(nknots-3)+1) = -8;
res(nknots*(nknots-2)+1,nknots*(nknots-2)) = res(nknots*(nknots-2)+1,nknots*(nknots-1)+1) = -6;
res(nknots*(nknots-2)+1,nknots*(nknots-3)) = res(nknots*(nknots-2)+1,nknots*(nknots-3)+2) = res(nknots*(nknots-2)+1,nknots*(nknots-1)) = res(nknots*(nknots-2)+1,nknots*(nknots-1)+2) = 2;
res(nknots*(nknots-2)+1,nknots*(nknots-2)+3) = res(nknots*(nknots-2)+1,nknots*(nknots-4)+1) = 1;

// lower right
res(nknots*(nknots-1)-2,nknots*(nknots-1)-2) = 18;
res(nknots*(nknots-1)-2,nknots*(nknots-1)-3) = res(nknots*(nknots-1)-2,nknots*(nknots-2)-2) = -8;
res(nknots*(nknots-1)-2,nknots*(nknots-1)-1) = res(nknots*(nknots-1)-2,nknots*nknots-2) = -6;
res(nknots*(nknots-1)-2,nknots*(nknots-2)-3) = res(nknots*(nknots-1)-2,nknots*(nknots-2)-1) = res(nknots*(nknots-1)-2,nknots*nknots-3) = res(nknots*(nknots-1)-2,nknots*nknots-1) = 2;
res(nknots*(nknots-1)-2,nknots*(nknots-3)-2) = res(nknots*(nknots-1)-2,nknots*(nknots-1)-4) = 1;

// 8 next to corner edges

// upper left
res(1,1) = 10;
res(1,2) = -6;
res(1,nknots+1) = -6;
res(1,0) = -4;
res(1,nknots) = res(1,nknots+2) = 2;
res(1,3) = res(1,2*nknots+1) = 1;

res(nknots,nknots) = 10;
res(nknots,2*nknots) = -6;
res(nknots,nknots+1) = -6;
res(nknots,0) = -4;
res(nknots,1) = res(nknots,2*nknots+1) = 2;
res(nknots,nknots+2) = res(nknots,3*nknots) = 1;

// upper right
res(nknots-2,nknots-2) = 10;
res(nknots-2,nknots-3) = -6;
res(nknots-2,2*nknots-2) = -6;
res(nknots-2,nknots-1) = -4;
res(nknots-2,2*nknots-3) = res(nknots-2,2*nknots-1) = 2;
res(nknots-2,nknots-4) = res(nknots-2,3*nknots-2) = 1;

res(2*nknots-1,2*nknots-1) = 10;
res(2*nknots-1,3*nknots-1) = -6;
res(2*nknots-1,2*nknots-2) = -6;
res(2*nknots-1,nknots-1) = -4;
res(2*nknots-1,nknots-2) = res(2*nknots-1,3*nknots-2) = 2;
res(2*nknots-1,2*nknots-3) = res(2*nknots-1,4*nknots-1) = 1;

// lower left
res(nknots*(nknots-2),nknots*(nknots-2)) = 10;
res(nknots*(nknots-2),nknots*(nknots-3)) = -6;
res(nknots*(nknots-2),nknots*(nknots-2)+1) = -6;
res(nknots*(nknots-2),nknots*(nknots-1)) = -4;
res(nknots*(nknots-2),nknots*(nknots-3)+1) = res(nknots*(nknots-2),nknots*(nknots-1)+1) = 2;
res(nknots*(nknots-2),nknots*(nknots-2)+2) = res(nknots*(nknots-2),nknots*(nknots-4)) = 1;

res(nknots*(nknots-1)+1,nknots*(nknots-1)+1) = 10;
res(nknots*(nknots-1)+1,nknots*(nknots-1)+2) = -6;
res(nknots*(nknots-1)+1,nknots*(nknots-2)+1) = -6;
res(nknots*(nknots-1)+1,nknots*(nknots-1)) = -4;
res(nknots*(nknots-1)+1,nknots*(nknots-2)) = res(nknots*(nknots-1)+1,nknots*(nknots-2)+2) = 2;
res(nknots*(nknots-1)+1,nknots*(nknots-1)+3) = res(nknots*(nknots-1)+1,nknots*(nknots-3)+1) = 1;

// lower right
res(nknots*(nknots-1)-1,nknots*(nknots-1)-1) = 10;
res(nknots*(nknots-1)-1,nknots*(nknots-2)-1) = -6;
res(nknots*(nknots-1)-1,nknots*(nknots-1)-2) = -6;
res(nknots*(nknots-1)-1,nknots*nknots-1) = -4;
res(nknots*(nknots-1)-1,nknots*nknots-2) = res(nknots*(nknots-1)-1,nknots*(nknots-2)-2) = 2;
res(nknots*(nknots-1)-1,nknots*(nknots-1)-3) = res(nknots*(nknots-1)-1,nknots*(nknots-3)-1) = 1;

res(nknots*nknots-2,nknots*nknots-2) = 10;
res(nknots*nknots-2,nknots*nknots-3) = -6;
res(nknots*nknots-2,nknots*(nknots-1)-2) = -6;
res(nknots*nknots-2,nknots*nknots-1) = -4;
res(nknots*nknots-2,nknots*(nknots-1)-3) = res(nknots*nknots-2,nknots*(nknots-1)-1) = 2;
res(nknots*nknots-2,nknots*nknots-4) = res(nknots*nknots-2,nknots*(nknots-2)-2) = 1;

for(i=2; i<nknots-2; i++)
  {
  // real edges

  // upper edge
  res(i,i) = 11;
  res(i,i-1) = res(i,i+1) = res(i,nknots+i) = -6;
  res(i,nknots+i-1) = res(i,nknots+i+1) = 2;
  res(i,i-2) = res(i,i+2) = res(i,2*nknots+i) = 1;

  // left edge
  res(i*nknots,i*nknots) = 11;
  res(i*nknots,i*nknots+1) = res(i*nknots,(i-1)*nknots) = res(i*nknots,(i+1)*nknots) = -6;
  res(i*nknots,(i-1)*nknots+1) = res(i*nknots,(i+1)*nknots+1) = 2;
  res(i*nknots,i*nknots+2) = res(i*nknots,(i-2)*nknots) = res(i*nknots,(i+2)*nknots) = 1;

  // right edge
  res((i+1)*nknots-1,(i+1)*nknots-1) = 11;
  res((i+1)*nknots-1,i*nknots-1) = res((i+1)*nknots-1,(i+1)*nknots-2) = res((i+1)*nknots-1,(i+2)*nknots-1) = -6;
  res((i+1)*nknots-1,i*nknots-2) = res((i+1)*nknots-1,(i+2)*nknots-2) = 2;
  res((i+1)*nknots-1,(i+1)*nknots-3) = res((i+1)*nknots-1,(i-1)*nknots-1) = res((i+1)*nknots-1,(i+3)*nknots-1) = 1;

  // lower edge
  res(nknots*(nknots-1)+i,nknots*(nknots-1)+i) = 11;
  res(nknots*(nknots-1)+i,nknots*(nknots-1)+i-1) = res(nknots*(nknots-1)+i,nknots*(nknots-1)+i+1) = res(nknots*(nknots-1)+i,nknots*(nknots-2)+i) = -6;
  res(nknots*(nknots-1)+i,nknots*(nknots-2)+i-1) = res(nknots*(nknots-1)+i,nknots*(nknots-2)+i+1) = 2;
  res(nknots*(nknots-1)+i,nknots*(nknots-1)+i-2) = res(nknots*(nknots-1)+i,nknots*(nknots-1)+i+2) = res(nknots*(nknots-1)+i,nknots*(nknots-3)+i) = 1;

  // second order edges

  // upper edge
  res(nknots+i,nknots+i) = 19;
  res(nknots+i,nknots+i-1) = res(nknots+i,nknots+i+1) = res(nknots+i,2*nknots+i) = -8;
  res(nknots+i,i) = -6;
  res(nknots+i,i-1) = res(nknots+i,i+1) = res(nknots+i,2*nknots+i-1) = res(nknots+i,2*nknots+i+1) = 2;
  res(nknots+i,nknots+i-2) = res(nknots+i,nknots+i+2) = res(nknots+i,3*nknots+i) = 1;

  // left edge
  res(i*nknots+1,i*nknots+1) = 19;
  res(i*nknots+1,i*nknots+2) = res(i*nknots+1,(i-1)*nknots+1) = res(i*nknots+1,(i+1)*nknots+1) = -8;
  res(i*nknots+1,i*nknots) = -6;
  res(i*nknots+1,(i-1)*nknots) = res(i*nknots+1,(i-1)*nknots+2) = res(i*nknots+1,(i+1)*nknots) = res(i*nknots+1,(i+1)*nknots+2) = 2;
  res(i*nknots+1,i*nknots+3) = res(i*nknots+1,(i-2)*nknots+1) = res(i*nknots+1,(i+2)*nknots+1) = 1;

  // right edge
  res((i+1)*nknots-2,(i+1)*nknots-2) = 19;
  res((i+1)*nknots-2,(i+1)*nknots-3) = res((i+1)*nknots-2,i*nknots-2) = res((i+1)*nknots-2,(i+2)*nknots-2) = -8;
  res((i+1)*nknots-2,(i+1)*nknots-1) = -6;
  res((i+1)*nknots-2,i*nknots-3) = res((i+1)*nknots-2,i*nknots-1) = res((i+1)*nknots-2,(i+2)*nknots-3) = res((i+1)*nknots-2,(i+2)*nknots-1) = 2;
  res((i+1)*nknots-2,(i+1)*nknots-4) = res((i+1)*nknots-2,(i-1)*nknots-2) = res((i+1)*nknots-2,(i+3)*nknots-2) = 1;

  // lower edge
  res(nknots*(nknots-2)+i,nknots*(nknots-2)+i) = 19;
  res(nknots*(nknots-2)+i,nknots*(nknots-2)+i-1) = res(nknots*(nknots-2)+i,nknots*(nknots-2)+i+1) = res(nknots*(nknots-2)+i,nknots*(nknots-3)+i) = -8;
  res(nknots*(nknots-2)+i,nknots*(nknots-1)+i) = -6;
  res(nknots*(nknots-2)+i,nknots*(nknots-3)+i-1) = res(nknots*(nknots-2)+i,nknots*(nknots-3)+i+1) = res(nknots*(nknots-2)+i,nknots*(nknots-1)+i-1) = res(nknots*(nknots-2)+i,nknots*(nknots-1)+i+1) = 2;
  res(nknots*(nknots-2)+i,nknots*(nknots-2)+i-2) = res(nknots*(nknots-2)+i,nknots*(nknots-2)+i+2) = res(nknots*(nknots-2)+i,nknots*(nknots-4)+i) = 1;
  }

// interior

for(i=2; i<nknots-2; i++)
  {
  for(j=2; j<nknots-2; j++)
    {
    res(i*nknots+j,i*nknots+j) = 20;
    res(i*nknots+j,i*nknots+j-1) = res(i*nknots+j,i*nknots+j+1) = res(i*nknots+j,(i-1)*nknots+j) = res(i*nknots+j,(i+1)*nknots+j) = -8;
    res(i*nknots+j,(i-1)*nknots+j-1) = res(i*nknots+j,(i-1)*nknots+j+1) = res(i*nknots+j,(i+1)*nknots+j-1) = res(i*nknots+j,(i+1)*nknots+j+1) = 2;
    res(i*nknots+j,i*nknots+j-2) = res(i*nknots+j,i*nknots+j+2) = res(i*nknots+j,(i-2)*nknots+j) = res(i*nknots+j,(i+2)*nknots+j) = 1;
    }
  }

ofstream out1("c:\\temp\\K.raw");
res.prettyPrint(out1);
out1.close();

/*
// Corners
  res(0,1) = 1;
  res(0,nknots) = 1;

  res(nknots-1,nknots-2) = 1;
  res(nknots-1,2*nknots-1) = 1;

  res(nknots*(nknots-1),nknots*(nknots-2)) = 1;
  res(nknots*(nknots-1),nknots*(nknots-1)+1) = 1;

  res(nknots*nknots-1,nknots*(nknots-1)-1) = 1;
  res(nknots*nknots-1,nknots*nknots-2) = 1;

// edges

  for(i=1; i<nknots-1; i++)
    {
    res(i,i-1) = 1;
    res(i,i+1) = 1;
    res(i,nknots+i) = 1;

    res(nknots*(nknots-1)+i,nknots*(nknots-1)+i-1) = 1;
    res(nknots*(nknots-1)+i,nknots*(nknots-1)+i+1) = 1;
    res(nknots*(nknots-1)+i,nknots*(nknots-2)+i) = 1;

    res(i*nknots,(i-1)*nknots) = 1;
    res(i*nknots,i*nknots+1) = 1;
    res(i*nknots,(i+1)*nknots) = 1;

    res((i+1)*nknots-1,i*nknots-1) = 1;
    res((i+1)*nknots-1,(i+1)*nknots-2) = 1;
    res((i+1)*nknots-1,(i+2)*nknots-1) = 1;
    }

// interior

  for(i=1; i<nknots-1; i++)
    {
    for(j=1; j<nknots-1; j++)
      {
      res(i*nknots+j ,i*nknots+j-1) = 1;
      res(i*nknots+j ,i*nknots+j+1) = 1;
      res(i*nknots+j ,(i-1)*nknots+j) = 1;
      res(i*nknots+j ,(i+1)*nknots+j) = 1;
      }
    }

// diagonal

  for(i=0; i<nknots*nknots; i++)
    {
    res(i,i) = -res.sum(i);
    }

  res = res.transposed()*res;
  */
  return res;
  }


}


