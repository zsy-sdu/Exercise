#include "molecular.h"
#include <iostream>

int main(int argc, char *argv[]) {
  if (argc != 2) {
    std::cerr << "Usage: " << argv[0] << " filename" << std::endl;
    return 1;
  }
  /*
  const char *filename = argv[1];
  Molecular molecule(filename);
  std::cout << molecule << std::endl;
  */

  return 0;
}