#include <stdio.h>
#include <string.h>

typedef unsigned int       uint32;
#define INLINE inline static
INLINE uint32 csl2yesnointarr ( uint32 * arr, char * str ) {
  char * p=str, *q=str;
  uint32 i=0;
  char ans[4];
  if (arr) { 
    while (*p) {
      if (*p==',') {
	*p=' ';sscanf(q,"%s",ans);q=p;
	if (ans[0]=='y'||ans[0]=='Y') arr[i]=0;
	else arr[i]=1;
	i++;
      }
      p++;
    }
    sscanf(q,"%s",ans);
    if (ans[0]=='y'||ans[0]=='Y') arr[i]=0;
    else arr[i]=1;
    i++;
  }
  return i;
}

int main () {

  uint32 arr[5];
  char str[255];
  int i,cnt;


  while (fscanf(stdin,"%s",str)) {
    cnt=csl2yesnointarr(arr,str);
    for (i=0;i<cnt;i++) fprintf(stdout,"%u ",arr[i]);
    fprintf(stdout,"\n");
  }

}
