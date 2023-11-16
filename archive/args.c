#include <stdio.h>

int main (int argc, char * argv[]) {
  int i;

  int n=strlen(argv[1]);
  char * p = argv[1];
  char * q;
  
  i=0;
  while ((int)(p-argv[1])<n) {
    while (isspace(*p)) p++;
    q = p+1;
    if ((int)*p==39) {
      while (*q&&(int)*q!=39) q++;
      if ((int)*q==39) q++;
      q[0]='\0';
    }
    else {
      while (*q && !isspace(*q)) q++;
      if (*q && isspace(*q)) q[0]='\0';
    }
    printf("%s ",p);
    p=q+1;
    i++;
  }
  //  printf("\n");
}
