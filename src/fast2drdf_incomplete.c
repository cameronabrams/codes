#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define SDIM 2

typedef struct CELL * cp;
typedef struct CELL {
   int * mem;
   int n;
   int i;
   int * nb;
   double ll[SDIM];
   double ur[SDIM];
   cp next;
} cell;

cell * new_cells ( double L[], double rcut ) {
   int nc[SDIM];
   double sp[SDIM];
   double ll[SDIM];
   double ur[SDIM];
   int i,j,k,ind;
   int tc=0;
   cell * C = NULL, * c = NULL;

   for (i=0;i<SDIM;i++) {
      nc[i] = (int)(L[i]/rcut);
      sp[i] = L[i]/nc[i];
      tc+=nc[i];
   }
   ind=0;
#if SDIM > 0
   for (i=0;i<nc[0];i++) {
     ll[0]=i*sp[0];
     ur[0]=(i+1)*sp[0];
#if SDIM > 1
   	for (j=0;j<nc[1];j++) {
        ll[1]=j*sp[1];
        ur[1]=(j+1)*sp[1];
#if SDIM > 2
	for (k=0;k<nc[2];k++) {
           ll[2]=k*sp[2];
           ur[2]=(k+1)*sp[2];
#endif
#endif
#endif
	   if (c) {
	      c->next=new_cell(ind,ll,ur);
	   } else {
	      C=new_cell(ind,ll,ur);
	      c=C;
	   }
	   ind++;
	   c=c->next;
#if SDIM > 2
	}
#if SDIM > 1
     }
#if SDIM > 0
   }
#endif
#endif
#endif
   return C;
}

cell * cellify ( cell * C, double r[][SDIM], int N, double L[], double rcut ) {
   if (!C) {
      C = new_cells(L,rcut);
   } else {
      empty_cells(C);
   }
   fill_cells(C,r,L,rcut);
   return C;
}

int main ( ) {
   
}
