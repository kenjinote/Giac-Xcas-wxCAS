// -*- compile-command: "g++ -g pgcd.cc -lgiac -lgmp" -*-
#include <giac/giac.h>

using namespace std;
using namespace giac;

int main(){
  cout << "Tapez une fonction a integrer et une variable ";
  gen a,b;
  cin >> a >> b;
  cout << "La primitive est " << integrate(a,b,0) << endl;
  return 0;
}
