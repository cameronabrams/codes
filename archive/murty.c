#include <stdio.h>
#include <math.h>
typedef float _fl_;
void MurtyAtwaterPlateauCutoff (_fl_ x, 
				_fl_ l_off, _fl_ l_on, _fl_ h_on, _fl_ h_off, 
				_fl_ * fcu, _fl_ * dfcu)
{
  if (x<l_off||x>h_off) {
    *fcu = *dfcu = 0.0;
  }
  else if (x<l_on) {
    _fl_ pd=M_PI/(l_on-l_off);
    _fl_ d=x-(l_on+l_off)/2.0;
    *fcu = 0.5 + 0.5625 * sin ( pd * d ) + 0.0625 * sin ( 3 * pd * d );
    *dfcu = 0.5625 * pd * cos ( pd * d ) + 0.1875 * pd * cos ( 3 * pd * d );
  }
  else if (x<h_on) {
    *fcu = 1.0;
    *dfcu = 0.0;
  }
  else if (x<h_off) {
    _fl_ pd=M_PI/(h_off-h_on);
    _fl_ d=x-(h_off+h_on)/2.0;
    *fcu = 0.5 - 0.5625 * sin ( pd * d ) - 0.0625 * sin ( 3 * pd * d );
    *dfcu = - 0.5625 * pd * cos ( pd * d ) - 0.1875 * pd * cos ( 3 * pd * d );
  }
}

int main ( int argc, char * argv[] )
{

  _fl_ cosTh, Th, fcu, dfcu, afcu, adfcu;
  _fl_ l_off = -0.9976, l_on = -0.98, h_on = 0.8480, h_off = 0.8829;
  _fl_ al_off = M_PI/180.0*28.0, al_on = M_PI/180.0*32.0, 
    ah_on = M_PI/180.0*172.0, ah_off = M_PI/180.0*176.0;
  
  for (cosTh=-1.0;cosTh<1.1;cosTh+=0.001) {
    //MurtyAtwaterPlateauCutoff(cosTh,l_off,l_on,h_on,h_off,&fcu,&dfcu);
    Th=acos(cosTh);
    MurtyAtwaterPlateauCutoff(Th,al_off,al_on,ah_on,ah_off,&afcu,&adfcu);
    fprintf(stdout,"%f %f %f\n",Th,afcu,adfcu);
  }
}
