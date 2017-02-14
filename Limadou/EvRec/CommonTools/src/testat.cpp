#include "statq.hh"
#include "iostream"
#include <vector>

int main(int argc, char *argv[]) {

   std::vector<float> v{-2, -4, -10};
   statq stat(v);
   std::cout << stat.GetStdDev() << std::endl;
   std::cout << 2+2 << std::endl;

  return 0;
}
