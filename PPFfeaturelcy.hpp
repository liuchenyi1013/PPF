#include <iostream>
using namespace std;

  /** \brief A point structure for storing the Point Pair Feature (PPF) values
    * \ingroup common
    */


  struct PPFSignaturelcy
  {
    float f1, f2, f3, f4;
    float alpha_m;
  
    friend std::ostream& operator << (std::ostream& os, const PPFSignaturelcy& p);
  };

  ostream& operator<<(ostream& out,const PPFSignaturelcy& p)
  {
    out<<p.f1<<','<<p.f2<<','<<p.f3<<","<<p.f4<<","<<p.alpha_m;;
   
    return out;
  }

  