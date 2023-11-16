#include "bits.h"

char __sprint_as_bits__[MAXBITSTRINGLENGTH];
char * sprint_as_bits ( void *a, uint s )
{
  int i,j;
  char * p=__sprint_as_bits__;
    
  uint mask=0;
  //printf("a 0x%8x contains %i bytes\n",a,s);
  for (i=s-1;i>-1;i--)
  {
    //printf("byte %i is %x\n",i,(*((char*)(a)+i)));
    /* for byte i, make a mask */
    mask=1<<7;
    for (j=0;j<8;j++)
    {
      //printf(" > bit %i (mask %x) is %c\n",j,mask,(*((char*)(a)+i))&mask?'1':'0');
      p+=sprintf(p,"%c",(*((char*)(a)+i))&mask?'1':'0');
      if (!((j+1)%4)) p+=sprintf(p," ");
      mask>>=1;
    }
  }

  return __sprint_as_bits__;
}

/* setbit
 * sets the b'th bit of a bit string *a of size s, 
 * assuming the leftmost bit is bit number `1'. */
void setbit ( uint * a, uint b, uint s )
{
   *a=(*a)|(1<<(s-b));
}

/* clearbit
 * clears the b'th bit of a bit string *a of size s, 
 * assuming the leftmost bit is bit number `1'. */
void clearbit ( uint * a, uint b, uint s )
{
   *a=(*a)&~(1<<(s-b));
}

/* flipbit
 * flips the b'th bit of a bit string *a of size s, 
 * assuming the leftmost bit is bit number `1'. */
void flipbit ( uint * a, uint b, uint s )
{
   (*a)=~((*a)^~(1<<(s-b)));
}
