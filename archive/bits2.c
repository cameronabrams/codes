/*#ifndef CGFLAGS_H
#define CGFLAGS_H

#include "libgeneric/itypes.h"
#include "libpolycfg/configuration.h"
#include <string.h>
#include "cgglobal.h"
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define INLINE inline static

typedef unsigned int uint32;

// if cgtype is SUPERATOM, then "group" identifies a mapping group;
// if cgtype is ATOM, then "group" is not used

char __sprint_as_bits__[70];
char * sprint_as_bits ( void *a, uint s )
{
  int i,j;
  char * p=__sprint_as_bits__;
    
  uint mask=0;
  //printf("a 0x%8x contains %i bytes\n",a,s);
  for (i=s-1;i>-1;i--) {
    //printf("byte %i is %x\n",i,(*((char*)(a)+i)));
    /* for byte i, make a mask */
    mask=1<<7;
    for (j=0;j<8;j++) {
      //printf(" > bit %i (mask %x) is %c\n",j,mask,(*((char*)(a)+i))&mask?'1':'0');
      p+=sprintf(p,"%c",(*((char*)(a)+i))&mask?'1':'0');
      if (!((j+1)%4)) p+=sprintf(p," ");
      mask>>=1;
    }
  }

  return __sprint_as_bits__;
}

/**********************************************************************************/

enum { SUPERATOM, ATOM_O, ATOM };

/* We are storing FIVE pieces of information in one 32-bit unsigned integer:
   1. bits 32 to 31 store the 'active' flag, which can have value 0 to 3;
   2. bits 30 to 29 store the 'cgtype' flag, which can have value 0 to 3;
   3. bits 28 to 25 store the 'group' flag, which can have values 0 to 15;
   4. bits 25 to 21 store the 'mark' flag, which can have values 0 to 3;
   5. bits 20 to 1  store the 'partner' flag, which can have values 0 to 4194303
*/

#define ACTIVE_BITS 30
#define CGTYPE_BITS 28
#define GROUP_BITS  24
#define MARK_BITS   22

#define ACTIVE_MASK ((1 << ACTIVE_BITS) -1)
#define CGTYPE_MASK ((1 << CGTYPE_BITS) -1)
#define GROUP_MASK  ((1 << GROUP_BITS) -1)
#define MARK_MASK   ((1 << MARK_BITS) -1)

INLINE uint32 CGFlag_str2ptype ( char * str ) {
  if (str) {
    if (!strcmp(str,"SUPERATOM")) return SUPERATOM;
    if (!strcmp(str,"ATOM_O")) return ATOM_O;
    if (!strcmp(str,"ATOM")) return ATOM;
  }
  return SUPERATOM;
}

INLINE void CGFlag_sprintf_ptype ( char * s, uint32 ptype ) {
  if (s) {
    if (ptype==SUPERATOM) sprintf(s,"SUPERATOM");
    else if (ptype==ATOM_O) sprintf(s,"ATOM_O");
    else if (ptype==ATOM) sprintf(s,"ATOM");
    else sprintf(s,"UNKNOWN");
  }
}

INLINE uint32 CGUser_GetActiveFlag ( uint32 * flag )
{
  return (uint32) (*flag >> ACTIVE_BITS);
}

INLINE uint32 CGUser_GetCGType ( uint32 * flag )
{
  return (uint32) ( ( *flag & ACTIVE_MASK ) >> CGTYPE_BITS );
}

INLINE uint32 CGUser_GetGroup ( uint32 * flag )
{
  return (uint32) (( *flag & CGTYPE_MASK ) >> GROUP_BITS);
}

INLINE uint32 CGUser_GetMark ( uint32 * flag )
{
  return (uint32) (( *flag & GROUP_MASK ) >> MARK_BITS);
}

INLINE uint32 CGUser_GetPartner ( uint32 * flag )
{
  return (uint32) ( *flag & MARK_MASK );
}

INLINE uint32 CGUser_SetActiveFlag ( uint32 * flag, uint32 val )
{
  if (val >> (32-ACTIVE_BITS)) {
    fprintf(stderr,"# warning: value %u for flag 'active' exceeds limit %u\n",
	    val,~(~0<<(32-ACTIVE_BITS)));
  }
  *flag = 
    ( val << ACTIVE_BITS ) |
    ( CGUser_GetCGType(flag) << CGTYPE_BITS ) |
    ( CGUser_GetGroup(flag) << GROUP_BITS ) |
    ( CGUser_GetMark(flag) << MARK_BITS) |
    ( CGUser_GetPartner(flag) );

  return *flag;
}

INLINE uint32 CGUser_SetCGType ( uint32 * flag, uint32 val )
{
  if (val >> (ACTIVE_BITS-CGTYPE_BITS)) {
    fprintf(stderr,"# warning: value %u for flag 'cgtype' exceeds limit %u\n",
	    val,~(~0<<(ACTIVE_BITS-CGTYPE_BITS)));
  } 
  *flag = 
    ( CGUser_GetActiveFlag(flag) << ACTIVE_BITS ) |
    ( val << CGTYPE_BITS ) |
    ( CGUser_GetGroup(flag) << GROUP_BITS ) |
    ( CGUser_GetMark(flag) << MARK_BITS) |
    ( CGUser_GetPartner(flag) );

  return *flag;
}

INLINE uint32 CGUser_SetGroup ( uint32 * flag, uint32 val )
{
  if (val >> (CGTYPE_BITS-GROUP_BITS)) {
    fprintf(stderr,"# warning: value %u for flag 'group' exceeds limit %u\n",
	    val,~(~0<<(CGTYPE_BITS-GROUP_BITS)));
  } 
  *flag = 
    ( CGUser_GetActiveFlag(flag) << ACTIVE_BITS ) |
    ( CGUser_GetCGType(flag) << CGTYPE_BITS ) |
    ( val << GROUP_BITS ) |
    ( CGUser_GetMark(flag) << MARK_BITS) |
    ( CGUser_GetPartner(flag) );

  return *flag;
}

INLINE uint32 CGUser_SetMark ( uint32 * flag, uint32 val )
{
  if (val >> (GROUP_BITS-MARK_BITS)) {
    fprintf(stderr,"# warning: value %u for flag 'mark' exceeds limit %u\n",
	    val,~(~0<<(GROUP_BITS-MARK_BITS)));
  } 
  *flag =
    ( CGUser_GetActiveFlag(flag) << ACTIVE_BITS ) |
    ( CGUser_GetCGType(flag) << CGTYPE_BITS ) |
    ( CGUser_GetGroup(flag) << GROUP_BITS ) |
    ( val << MARK_BITS) |
    ( CGUser_GetPartner(flag) );

  return *flag;
}

INLINE uint32 CGUser_SetPartner (uint32 * flag, uint32 val )
{
  if (val >> (MARK_BITS)) {
    fprintf(stderr,"# warning: value %u for flag 'partner' exceeds limit %u\n",
	    val,~(~0<<(MARK_BITS)));
  } 
  *flag =
    ( CGUser_GetActiveFlag(flag) << ACTIVE_BITS ) |
    ( CGUser_GetCGType(flag) << CGTYPE_BITS ) |
    ( CGUser_GetGroup(flag) << GROUP_BITS ) |
    ( CGUser_GetMark(flag) << MARK_BITS) |
    ( val );

  return *flag;
}

/*INLINE uint32 CGUser_nInactives ( TConfig * cfg ) {

  PRIVATE(TLocalConfig,lCfg,cfg->hLocalCfg);
  int i;
  uint32 cnt=0;
  uint32 active;

  for(i = 0;i < lCfg->nLocal;i++) {
    active=CGUser_GetActiveFlag(&lCfg->label[i].user);
    if (active) cnt++; // anything other than 0 means that the atom is INACTIVE
  }
  return cnt;
}
*/

int main (int argc, char * argv[]) {

  uint a=1, t=1, g=1, m=1,p=1,i, flag=0;
  
  for (i=1;i<argc;i++) {
    if (!strcmp(argv[i],"-a")) a=atoi(argv[++i]);
    else if (!strcmp(argv[i],"-t")) t=atoi(argv[++i]);
    else if (!strcmp(argv[i],"-g")) g=atoi(argv[++i]);
    else if (!strcmp(argv[i],"-m")) m=atoi(argv[++i]);
    else if (!strcmp(argv[i],"-p")) p=atoi(argv[++i]);
  }

  /*  a=ACTIVE_MASK;
  fprintf(stdout,"a %s\n",sprint_as_bits((uint*)&a,sizeof(uint)));
  t=CGTYPE_MASK;
  fprintf(stdout,"t %s\n",sprint_as_bits((uint*)&t,sizeof(uint)));
  g=GROUP_MASK; 
  fprintf(stdout,"g %s\n",sprint_as_bits((uint*)&g,sizeof(uint)));
  m=MARK_MASK;
  fprintf(stdout,"m %s\n",sprint_as_bits((uint*)&m,sizeof(uint)));
  */
  CGUser_SetActiveFlag(&flag,a);
  fprintf(stdout,"set active flag to %u; get returns %u\n",a,CGUser_GetActiveFlag(&flag));
  fprintf(stdout,"f %s\n",sprint_as_bits((uint*)&flag,sizeof(uint)));
  CGUser_SetCGType(&flag,t);
  fprintf(stdout,"set cgtype to %u; get returns %u\n",t,CGUser_GetCGType(&flag));
  fprintf(stdout,"f %s\n",sprint_as_bits((uint*)&flag,sizeof(uint)));
  CGUser_SetGroup(&flag,g);
  fprintf(stdout,"set group to %u; get returns %u\n",g,CGUser_GetGroup(&flag));
  fprintf(stdout,"f %s\n",sprint_as_bits((uint*)&flag,sizeof(uint)));
  CGUser_SetMark(&flag,m);
  fprintf(stdout,"set mark to %u; get returns %u\n",m,CGUser_GetMark(&flag));
  fprintf(stdout,"f %s\n",sprint_as_bits((uint*)&flag,sizeof(uint)));
  CGUser_SetPartner(&flag,p);
  fprintf(stdout,"set partner to %u; get returns %u\n",p,CGUser_GetPartner(&flag));
  fprintf(stdout,"f %s\n",sprint_as_bits((uint*)&flag,sizeof(uint)));
}


