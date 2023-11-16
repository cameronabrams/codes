/* Wichman - Hill pseudo random number generator in [0, 1) */
/* *x muss auf ein int-array mit 3 Eintraegen zeigen       */
/* Wiederholrate approx. 7*10e+12                          */
/* aus 'Random Processes in Physical Systems'              */
/* by Charles A. Whitney, Wiley & Sons Publ. 1990          */

float rng(int *x)
  {
   float t ;

   *(x++) = (171*(*x))%30269 ;
   *(x++) = (172*(*x))%30309 ;
   *(x  ) = (170*(*x))%30323 ;
   t  = ((float) *(x--))/30323 ;
   t += ((float) *(x--))/30309 ;
   t += ((float) *(x  ))/30269 ;
   return(t - (int) t) ;
  }
