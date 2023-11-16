#include <stdio.h>
#include <string.h>

typedef unsigned int       uint32;
#define INLINE inline static

enum { SUPERATOM, ATOM_O, ATOM, BOTH, NULL_CGFLAG };
static char * CGFLAGNAMES[NULL_CGFLAG] = {"SUPERATOM", "ATOM", "ATOM_O", "BOTH"};

#define ACTIVE_BITS 30
#define CGTYPE_BITS 28
#define GROUP_BITS  24
#define MARK_BITS   22
#define THERMO_BITS 20

#define ACTIVE_MASK ((1 << ACTIVE_BITS) -1)
#define CGTYPE_MASK ((1 << CGTYPE_BITS) -1)
#define GROUP_MASK  ((1 << GROUP_BITS) -1)
#define MARK_MASK   ((1 << MARK_BITS) -1)
#define THERMO_MASK ((1 << THERMO_BITS) -1)

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

INLINE uint32 CGUser_GetThermo ( uint32 * flag )
{
  return (uint32) (( *flag & MARK_MASK ) >> THERMO_BITS);
}

INLINE uint32 CGUser_GetPartner ( uint32 * flag )
{
  return (uint32) ( *flag & THERMO_MASK );
}

INLINE void CGUser_fprintf ( FILE * fp, uint32 * flag )
{
  fprintf(fp,"a(%u) t(%u) g(%u) m(%u) th(%u) p(%u)\n",
	  CGUser_GetActiveFlag(flag),CGUser_GetCGType(flag),
	  CGUser_GetGroup(flag),CGUser_GetMark(flag),
	  CGUser_GetThermo(flag),CGUser_GetPartner(flag));
  fflush(fp);
}

INLINE uint32 CGUser_SetActiveFlag ( uint32 * flag, uint32 val )
{
  if (val >> (32-ACTIVE_BITS)) {
    fprintf(stderr,"# Error: value %u for flag 'active' exceeds limit %u.  Using 0.\n",
	    val,~(~0<<(32-ACTIVE_BITS)));
    val = 0;    
  }
  *flag = 
    ( val << ACTIVE_BITS ) |
    ( CGUser_GetCGType(flag) << CGTYPE_BITS ) |
    ( CGUser_GetGroup(flag) << GROUP_BITS ) |
    ( CGUser_GetMark(flag) << MARK_BITS) |
    ( CGUser_GetThermo(flag) << THERMO_BITS) |
    ( CGUser_GetPartner(flag) );

  return *flag;
}

INLINE uint32 CGUser_SetCGType ( uint32 * flag, uint32 val )
{
  if (val >> (ACTIVE_BITS-CGTYPE_BITS)) {
    fprintf(stderr,"# warning: value %u for flag 'cgtype' exceeds limit %u. Using 0.\n",
	    val,~(~0<<(ACTIVE_BITS-CGTYPE_BITS)));
    val = 0;    
  } 
  *flag = 
    ( CGUser_GetActiveFlag(flag) << ACTIVE_BITS ) |
    ( val << CGTYPE_BITS ) |
    ( CGUser_GetGroup(flag) << GROUP_BITS ) |
    ( CGUser_GetMark(flag) << MARK_BITS) |
    ( CGUser_GetThermo(flag) << THERMO_BITS) |
    ( CGUser_GetPartner(flag) );

  return *flag;
}

INLINE uint32 CGUser_SetGroup ( uint32 * flag, uint32 val )
{
  if (val >> (CGTYPE_BITS-GROUP_BITS)) {
    fprintf(stderr,"# warning: value %u for flag 'group' exceeds limit %u. Using 0.\n",
	    val,~(~0<<(CGTYPE_BITS-GROUP_BITS)));
    val = 0;    
  } 
  *flag = 
    ( CGUser_GetActiveFlag(flag) << ACTIVE_BITS ) |
    ( CGUser_GetCGType(flag) << CGTYPE_BITS ) |
    ( val << GROUP_BITS ) |
    ( CGUser_GetMark(flag) << MARK_BITS) |
    ( CGUser_GetThermo(flag) << THERMO_BITS) |
    ( CGUser_GetPartner(flag) );

  return *flag;
}

INLINE uint32 CGUser_SetMark ( uint32 * flag, uint32 val )
{
  if (val >> (GROUP_BITS-MARK_BITS)) {
    fprintf(stderr,"# warning: value %u for flag 'mark' exceeds limit %u. Using 0.\n",
	    val,~(~0<<(GROUP_BITS-MARK_BITS)));
    val = 0;    
  } 
  *flag =
    ( CGUser_GetActiveFlag(flag) << ACTIVE_BITS ) |
    ( CGUser_GetCGType(flag) << CGTYPE_BITS ) |
    ( CGUser_GetGroup(flag) << GROUP_BITS ) |
    ( val << MARK_BITS) |
    ( CGUser_GetThermo(flag) << THERMO_BITS) |
    ( CGUser_GetPartner(flag) );

  return *flag;
}

INLINE uint32 CGUser_SetThermo ( uint32 * flag, uint32 val )
{
  if (val >> (MARK_BITS-THERMO_BITS)) {
    fprintf(stderr,"# warning: value %u for flag 'thermo' exceeds limit %u. Using 0.\n",
	    val,~(~0<<(MARK_BITS-THERMO_BITS)));
    val = 0;    
  }
  *flag =
    ( CGUser_GetActiveFlag(flag) << ACTIVE_BITS ) |
    ( CGUser_GetCGType(flag) << CGTYPE_BITS ) |
    ( CGUser_GetGroup(flag) << GROUP_BITS ) |
    ( CGUser_GetMark(flag) << MARK_BITS) |
    ( val << THERMO_BITS) |
    ( CGUser_GetPartner(flag) );

  return *flag;
}

INLINE uint32 CGUser_SetPartner (uint32 * flag, uint32 val )
{
  if (val >> (THERMO_BITS)) {
    fprintf(stderr,"# warning: value %u for flag 'partner' exceeds limit %u. Using 0.\n",
	    val,~(~0<<(THERMO_BITS)));
    val = 0;    
  } 
  *flag =
    ( CGUser_GetActiveFlag(flag) << ACTIVE_BITS ) |
    ( CGUser_GetCGType(flag) << CGTYPE_BITS ) |
    ( CGUser_GetGroup(flag) << GROUP_BITS ) |
    ( CGUser_GetMark(flag) << MARK_BITS) |
    ( CGUser_GetThermo(flag) << THERMO_BITS) |
    ( val );

  return *flag;
}

int main () {

  uint32 flag = 0;
  uint32 val;
  char code[3];

  char ln[255];

  CGUser_fprintf(stdout,&flag);

  while (fgets(ln,255,stdin)) {
    sscanf(ln,"%u %s",&val,code);
    if (!strcmp(code,"a")) CGUser_SetActiveFlag(&flag,val);
    if (!strcmp(code,"t")) CGUser_SetCGType(&flag,val);
    if (!strcmp(code,"g")) CGUser_SetGroup(&flag,val);
    if (!strcmp(code,"m")) CGUser_SetMark(&flag,val);
    if (!strcmp(code,"th")) CGUser_SetThermo(&flag,val);
    if (!strcmp(code,"p")) CGUser_SetPartner(&flag,val);
    if (!strcmp(code,"r")) memset(&flag,0,sizeof(uint32));
    CGUser_fprintf(stdout,&flag);
    
  }

}
