#include <statmat_penalty.h>

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

}

/*
statmatrix<double> K2dim_pspline_rw2(const unsigned & nknots)
  {
  statmatrix<double> res(nknots*nknots,nknots*nknots,0);
//  unsigned i,j;

  statmatrix<double> D = matrix.diffmat(2,nknots);  ? // difference matrix of 2nd order

  statmatrix<double> DD = matrix.mult(t(D),D) ?       // D'D
  
  statmatrix<double> P1 = matrix.kronecker(I,t(D)*D); ?   // penalty matrix of x-direction
  statmatrix<double> P2 = matrix.kronecker(t(D)*D,I); ?   // penalty matrix of y-direction
  res = P1 + P2;

  return res;
  }
