#ifndef MOLECULAR_H
#define MOLECULAR_H

#include <iostream>
#include <string>
#include <vector>

struct Atom {
  std::string type; // element type
  std::size_t Z;    // atomic number
  double x, y, z;   // coordinates
};

class Molecular {
public:
  // Default constructor
  Molecular() = default;
  // Constructor with filename
  explicit Molecular(const std::string &filename);
  // Copy constructor
  Molecular(const Molecular &other);
  // Move constructor
  Molecular(Molecular &&other);
  // Destructor
  ~Molecular() = default;

  // Getters
  const std::vector<Atom> get_atoms() const;
  const int get_charge() const;
  const int get_spin() const;
  const std::string get_basis() const;

  // Setters
  void set_atoms(const std::vector<Atom> &atoms);
  void set_charge(const int charge);
  void set_spin(const int spin);
  void set_basis(const std::string &basis);

  // Overload << operator
  friend std::ostream &operator<<(std::ostream &os, const Molecular &molecule);

private:
  std::vector<Atom> atoms_; // atoms in the molecule
  int charge_;              // charge of the molecule
  int spin_;                // spin multiplicity of the molecule
  std::string basis_;       // basis set name
}; // class Molecular
#endif // MOLECULAR_H