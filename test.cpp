#include <iostream>
#include <cmath>
#include <vector>
using namespace std;

int main() {
  vector< int > v;
  v.push_back(0);
  cout << "Vector v has following elements: " << endl;
  for (int i=0; i<v.size(); i++) {
    cout << v[i] << endl;
  }
  double x=-0.2, y=1.0;
  cout << "x/y " << x/y << endl;
  cout << "(floor)x/y " << floor(x/y) << endl;
  return 0;
}
