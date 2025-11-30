#pragma once
#include <vector>
#ifndef GMAT_INTEGRAL_HPP_
#define GMAT_INTEGRAL_HPP_

#include <complex>
#include <iostream>
#include <omp.h>
#include <span>
#include <type_traits>

extern "C" {
#include <cint.h>
int cint1e_ovlp_cart(double *buf, int *shls, int *atm, int natm, int *bas,
                     int nbas, double *env);
int cint1e_nuc_cart(double *buf, int *shls, int *atm, int natm, int *bas,
                    int nbas, double *env);
int cint1e_kin_cart(double *buf, int *shls, int *atm, int natm, int *bas,
                    int nbas, double *env);
int cint2e_cart(double *buf, int *shls, int *atm, int natm, int *bas, int nbas,
                double *env, CINTOpt *opt);
int cint1e_ovlp_sph(double *buf, int *shls, int *atm, int natm, int *bas,
                    int nbas, double *env);
int cint1e_nuc_sph(double *buf, int *shls, int *atm, int natm, int *bas,
                   int nbas, double *env);
int cint1e_kin_sph(double *buf, int *shls, int *atm, int natm, int *bas,
                   int nbas, double *env);
int cint2e_sph(double *buf, int *shls, int *atm, int natm, int *bas, int nbas,
               double *env, CINTOpt *opt);
int cint1e_ovlp(double *buf, int *shls, int *atm, int natm, int *bas, int nbas,
                double *env);
int cint1e_nuc(double *buf, int *shls, int *atm, int natm, int *bas, int nbas,
               double *env);
int cint1e_kin(double *buf, int *shls, int *atm, int natm, int *bas, int nbas,
               double *env);
int cint2e(double *buf, int *shls, int *atm, int natm, int *bas, int nbas,
           double *env, CINTOpt *opt);
}

namespace core {
enum class FunctionType { SPH, CART, SPINOR }; // FunctionType

template <FunctionType FT>
  requires(FT == FunctionType::CART || FT == FunctionType::SPH ||
           FT == FunctionType::SPINOR)
struct Integral {

  using T = std::conditional_t<FT == FunctionType::SPINOR, std::complex<double>,
                               double>;
  // Libcint Data
  std::span<int> atm_;
  std::span<int> bas_;
  std::span<double> env_;
  CINTOpt *opt_ = nullptr; // managed by libcint

  int natm_ = 0;
  int nbas_ = 0;
  int nao_ = 0; // number of atomic orbitals

  Integral() = delete;
  Integral(int natm, int *atm, int nbas, int *bas, double *env, int env_size)
      : atm_(atm, natm * 6), // libcint: atm is [natm][6]
        bas_(bas, nbas * 8), // libcint: bas is [nbas][8]
        env_(env, env_size), natm_(natm), nbas_(nbas) {
    for (int i = 0; i < nbas_; ++i) {
      if constexpr (FT == FunctionType::CART) {
        nao_ += CINTcgto_cart(i, bas_.data());
      } else if constexpr (FT == FunctionType::SPH) {
        nao_ += CINTcgto_spheric(i, bas_.data());
      } else if constexpr (FT == FunctionType::SPINOR) {
        nao_ += CINTcgto_spinor(i, bas_.data());
      } // FunctionType
    }

    std::cout << "Integral: natm = " << natm_ << ", nbas = " << nbas_
              << ", env_size = " << env_size << ", nao = " << nao_ << std::endl;
  }

  ~Integral() {
    if (opt_ != nullptr) {
      CINTdel_optimizer(&opt_);
    }
  };

  T *ovlp() { // overlap integral in row_major format
    T *ovlp_mat = new T[nao_ * nao_];
    memset(ovlp_mat, 0, sizeof(T) * nao_ * nao_);
#pragma omp parallel for collapse(2)
    for (int i = 0; i < nbas_; ++i) {
      for (int j = i; j < nbas_; ++j) {
        int *shls = new int[2]{i, j};
        int di, dj;
        int x, y;
        if constexpr (FT == FunctionType::CART) {
          di = CINTcgto_cart(i, bas_.data());
          dj = CINTcgto_cart(j, bas_.data());
          x = CINTtot_cgto_cart(bas_.data(), i);
          y = CINTtot_cgto_cart(bas_.data(), j);
        } else if constexpr (FT == FunctionType::SPH) {
          di = CINTcgto_spheric(i, bas_.data());
          dj = CINTcgto_spheric(j, bas_.data());
          x = CINTtot_cgto_spheric(bas_.data(), i);
          y = CINTtot_cgto_spheric(bas_.data(), j);
        } else if constexpr (FT == FunctionType::SPINOR) {
          di = CINTcgto_spinor(i, bas_.data());
          dj = CINTcgto_spinor(j, bas_.data());
          x = CINTtot_cgto_spinor(bas_.data(), i);
          y = CINTtot_cgto_spinor(bas_.data(), j);
        } // FunctionType
        double *buf;
        if constexpr (FT == FunctionType::SPINOR) {
          buf = new double[2 * di * dj];
        } else {
          buf = new double[di * dj];
        }

        if constexpr (FT == FunctionType::CART) {
#if defined(_DEBUG_)
          if (!cint1e_ovlp_cart(buf, shls, atm_.data(), natm_, bas_.data(),
                                nbas_, env_.data()))
            throw std::runtime_error("cint1e_ovlp_cart failed");
#else
          cint1e_ovlp_cart(buf, shls, atm_.data(), natm_, bas_.data(), nbas_,
                           env_.data());
#endif // _DEBUG_
        } else if constexpr (FT == FunctionType::SPH) {
#if defined(_DEBUG_)
          if (!cint1e_ovlp_sph(buf, shls, atm_.data(), natm_, bas_.data(),
                               nbas_, env_.data()))
            throw std::runtime_error("cint1e_ovlp_sph failed");
#else
          cint1e_ovlp_sph(buf, shls, atm_.data(), natm_, bas_.data(), nbas_,
                          env_.data());
#endif // _DEBUG_
        } else if constexpr (FT == FunctionType::SPINOR) {
#if defined(_DEBUG_)
          if (!cint1e_ovlp(buf, shls, atm_.data(), natm_, bas_.data(), nbas_,
                           env_.data()))
            throw std::runtime_error("cint1e_ovlp failed");
#else
          cint1e_ovlp(buf, shls, atm_.data(), natm_, bas_.data(), nbas_,
                      env_.data());
#endif // _DEBUG_
        }
        // copy buf to ovlp_mat

        for (int ii = 0; ii < di; ++ii) {
          for (int jj = 0; jj < dj; ++jj) {
            T val;
            if constexpr (FT == FunctionType::SPINOR) {
              val = std::complex<double>(buf[2 * (ii + jj * di)],
                                         buf[2 * (ii + jj * di) + 1]);
            } else {
              val = buf[ii + jj * di];
            }
            ovlp_mat[(x + ii) * nao_ + (y + jj)] = val;
            if constexpr (std::is_same_v<T, std::complex<double>>) {
              ovlp_mat[(y + jj) * nao_ + (x + ii)] = std::conj(val);
            } else {
              ovlp_mat[(y + jj) * nao_ + (x + ii)] = val;
            }
          }
        }

        delete[] buf;
        delete[] shls;
      }
    }

    return ovlp_mat;
  }
  T *kin() { // kinetic integral in row_major format
    T *kin_mat = new T[nao_ * nao_];
    memset(kin_mat, 0, sizeof(T) * nao_ * nao_);
#pragma omp parallel for collapse(2)
    for (int i = 0; i < nbas_; ++i) {
      for (int j = i; j < nbas_; ++j) {
        int *shls = new int[2]{i, j};
        int di, dj;
        int x, y;
        if constexpr (FT == FunctionType::CART) {
          di = CINTcgto_cart(i, bas_.data());
          dj = CINTcgto_cart(j, bas_.data());
          x = CINTtot_cgto_cart(bas_.data(), i);
          y = CINTtot_cgto_cart(bas_.data(), j);
        } else if constexpr (FT == FunctionType::SPH) {
          di = CINTcgto_spheric(i, bas_.data());
          dj = CINTcgto_spheric(j, bas_.data());
          x = CINTtot_cgto_spheric(bas_.data(), i);
          y = CINTtot_cgto_spheric(bas_.data(), j);
        } else if constexpr (FT == FunctionType::SPINOR) {
          di = CINTcgto_spinor(i, bas_.data());
          dj = CINTcgto_spinor(j, bas_.data());
          x = CINTtot_cgto_spinor(bas_.data(), i);
          y = CINTtot_cgto_spinor(bas_.data(), j);
        } // FunctionType
        double *buf;
        if constexpr (FT == FunctionType::SPINOR) {
          buf = new double[2 * di * dj];
        } else {
          buf = new double[di * dj];
        }

        if constexpr (FT == FunctionType::CART) {
#if defined(_DEBUG_)
          if (!cint1e_kin_cart(buf, shls, atm_.data(), natm_, bas_.data(),
                               nbas_, env_.data()))
            throw std::runtime_error("cint1e_kin_cart failed");
#else
          cint1e_kin_cart(buf, shls, atm_.data(), natm_, bas_.data(), nbas_,
                          env_.data());
#endif // _DEBUG_
        } else if constexpr (FT == FunctionType::SPH) {
#if defined(_DEBUG_)
          if (!cint1e_kin_sph(buf, shls, atm_.data(), natm_, bas_.data(), nbas_,
                              env_.data()))
            throw std::runtime_error("cint1e_kin_sph failed");
#else
          cint1e_kin_sph(buf, shls, atm_.data(), natm_, bas_.data(), nbas_,
                         env_.data());
#endif // _DEBUG_
        } else if constexpr (FT == FunctionType::SPINOR) {
#if defined(_DEBUG_)
          if (!cint1e_kin(buf, shls, atm_.data(), natm_, bas_.data(), nbas_,
                          env_.data()))
            throw std::runtime_error("cint1e_kin failed");
#else
          cint1e_kin(buf, shls, atm_.data(), natm_, bas_.data(), nbas_,
                     env_.data());
#endif // _DEBUG_
        }
        // copy buf to kin_mat

        for (int ii = 0; ii < di; ++ii) {
          for (int jj = 0; jj < dj; ++jj) {
            T val;
            if constexpr (FT == FunctionType::SPINOR) {
              val = std::complex<double>(buf[2 * (ii + jj * di)],
                                         buf[2 * (ii + jj * di) + 1]);
            } else {
              val = buf[ii + jj * di];
            }
            kin_mat[(x + ii) * nao_ + (y + jj)] = val;
            if constexpr (std::is_same_v<T, std::complex<double>>) {
              kin_mat[(y + jj) * nao_ + (x + ii)] = std::conj(val);
            } else {
              kin_mat[(y + jj) * nao_ + (x + ii)] = val;
            }
          }
        }
        delete[] buf;
        delete[] shls;
      }
    }

    return kin_mat;
  }
  T *nuc() { // nuclear integral in row_major format
    T *nuc_mat = new T[nao_ * nao_];
    memset(nuc_mat, 0, sizeof(T) * nao_ * nao_);
#pragma omp parallel for collapse(2)
    for (int i = 0; i < nbas_; ++i) {
      for (int j = i; j < nbas_; ++j) {
        int *shls = new int[2]{i, j};
        int di, dj;
        int x, y;
        if constexpr (FT == FunctionType::CART) {
          di = CINTcgto_cart(i, bas_.data());
          dj = CINTcgto_cart(j, bas_.data());
          x = CINTtot_cgto_cart(bas_.data(), i);
          y = CINTtot_cgto_cart(bas_.data(), j);
        } else if constexpr (FT == FunctionType::SPH) {
          di = CINTcgto_spheric(i, bas_.data());
          dj = CINTcgto_spheric(j, bas_.data());
          x = CINTtot_cgto_spheric(bas_.data(), i);
          y = CINTtot_cgto_spheric(bas_.data(), j);
        } else if constexpr (FT == FunctionType::SPINOR) {
          di = CINTcgto_spinor(i, bas_.data());
          dj = CINTcgto_spinor(j, bas_.data());
          x = CINTtot_cgto_spinor(bas_.data(), i);
          y = CINTtot_cgto_spinor(bas_.data(), j);
        } // FunctionType
        double *buf;
        if constexpr (FT == FunctionType::SPINOR) {
          buf = new double[2 * di * dj];
        } else {
          buf = new double[di * dj];
        }

        if constexpr (FT == FunctionType::CART) {
#if defined(_DEBUG_)
          if (!cint1e_nuc_cart(buf, shls, atm_.data(), natm_, bas_.data(),
                               nbas_, env_.data()))
            throw std::runtime_error("cint1e_nuc_cart failed");
#else
          cint1e_nuc_cart(buf, shls, atm_.data(), natm_, bas_.data(), nbas_,
                          env_.data());
#endif // _DEBUG_
        } else if constexpr (FT == FunctionType::SPH) {
#if defined(_DEBUG_)
          if (!cint1e_nuc_sph(buf, shls, atm_.data(), natm_, bas_.data(), nbas_,
                              env_.data()))
            throw std::runtime_error("cint1e_nuc_sph failed");
#else
          cint1e_nuc_sph(buf, shls, atm_.data(), natm_, bas_.data(), nbas_,
                         env_.data());
#endif // _DEBUG_
        } else if constexpr (FT == FunctionType::SPINOR) {
#if defined(_DEBUG_)
          if (!cint1e_nuc(buf, shls, atm_.data(), natm_, bas_.data(), nbas_,
                          env_.data()))
            throw std::runtime_error("cint1e_nuc failed");
#else
          cint1e_nuc(buf, shls, atm_.data(), natm_, bas_.data(), nbas_,
                     env_.data());
#endif // _DEBUG_
        }
        // copy buf to nuc_mat

        for (int ii = 0; ii < di; ++ii) {
          for (int jj = 0; jj < dj; ++jj) {
            T val;
            if constexpr (FT == FunctionType::SPINOR) {
              val = std::complex<double>(buf[2 * (ii + jj * di)],
                                         buf[2 * (ii + jj * di) + 1]);
            } else {
              val = buf[ii + jj * di];
            }
            nuc_mat[(x + ii) * nao_ + (y + jj)] = val;
            if constexpr (std::is_same_v<T, std::complex<double>>) {
              nuc_mat[(y + jj) * nao_ + (x + ii)] = std::conj(val);
            } else {
              nuc_mat[(y + jj) * nao_ + (x + ii)] = val;
            }
          }
        }
        delete[] buf;
        delete[] shls;
      }
    }

    return nuc_mat;
  }
  T *int2e_s8() { // 8-fold symmetry two-electron integral in row_major format
    T *int2e_mat = new T[nao_ * nao_ * nao_ * nao_];
    memset(int2e_mat, 0, sizeof(T) * nao_ * nao_ * nao_ * nao_);
    if (opt_ != nullptr) {
      CINTdel_optimizer(&opt_);
      opt_ = nullptr;
    }
    if constexpr (FT == FunctionType::CART) {
      cint2e_cart_optimizer(&opt_, atm_.data(), natm_, bas_.data(), nbas_,
                            env_.data());
    } else if constexpr (FT == FunctionType::SPH) {
      cint2e_sph_optimizer(&opt_, atm_.data(), natm_, bas_.data(), nbas_,
                           env_.data());
    } else if constexpr (FT == FunctionType::SPINOR) {
      cint2e_optimizer(&opt_, atm_.data(), natm_, bas_.data(), nbas_,
                       env_.data());
    } // FunctionType

#pragma omp parallel for collapse(4)
    for (int i = 0; i < nbas_; ++i) {
      for (int j = 0; j <= i; ++j) {
        for (int k = 0; k <= i; ++k) {
          const auto l_max = (i == k) ? j : k;
          for (int l = 0; l <= l_max; ++l) {
            int *shls = new int[4]{i, j, k, l};
            int di, dj, dk, dl;
            int x, y, z, w;
            if constexpr (FT == FunctionType::CART) {
              di = CINTcgto_cart(i, bas_.data());
              dj = CINTcgto_cart(j, bas_.data());
              dk = CINTcgto_cart(k, bas_.data());
              dl = CINTcgto_cart(l, bas_.data());
              x = CINTtot_cgto_cart(bas_.data(), i);
              y = CINTtot_cgto_cart(bas_.data(), j);
              z = CINTtot_cgto_cart(bas_.data(), k);
              w = CINTtot_cgto_cart(bas_.data(), l);
            } else if constexpr (FT == FunctionType::SPH) {
              di = CINTcgto_spheric(i, bas_.data());
              dj = CINTcgto_spheric(j, bas_.data());
              dk = CINTcgto_spheric(k, bas_.data());
              dl = CINTcgto_spheric(l, bas_.data());
              x = CINTtot_cgto_spheric(bas_.data(), i);
              y = CINTtot_cgto_spheric(bas_.data(), j);
              z = CINTtot_cgto_spheric(bas_.data(), k);
              w = CINTtot_cgto_spheric(bas_.data(), l);
            } else if constexpr (FT == FunctionType::SPINOR) {
              di = CINTcgto_spinor(i, bas_.data());
              dj = CINTcgto_spinor(j, bas_.data());
              dk = CINTcgto_spinor(k, bas_.data());
              dl = CINTcgto_spinor(l, bas_.data());
              x = CINTtot_cgto_spinor(bas_.data(), i);
              y = CINTtot_cgto_spinor(bas_.data(), j);
              z = CINTtot_cgto_spinor(bas_.data(), k);
              w = CINTtot_cgto_spinor(bas_.data(), l);
            } // FunctionType
            double *buf;
            if constexpr (FT == FunctionType::SPINOR) {
              buf = new double[2 * di * dj * dk * dl];
            } else {
              buf = new double[di * dj * dk * dl];
            }

            if constexpr (FT == FunctionType::CART) {
#if defined(_DEBUG_)
              if (!cint2e_cart(buf, shls, atm_.data(), natm_, bas_.data(),
                               nbas_, env_.data(), opt_))
                throw std::runtime_error("cint2e_cart failed");
#else
              cint2e_cart(buf, shls, atm_.data(), natm_, bas_.data(), nbas_,
                          env_.data(), opt_));
#endif // _DEBUG_
            } else if constexpr (FT == FunctionType::SPH) {
#if defined(_DEBUG_)
              if (!cint2e_sph(buf, shls, atm_.data(), natm_, bas_.data(), nbas_,
                              env_.data(), opt_))
                throw std::runtime_error("cint2e_sph failed");
#else
              cint2e_sph(buf, shls, atm_.data(), natm_, bas_.data(), nbas_,
                              env_.data(), opt_));
#endif // _DEBUG_
            } else if constexpr (FT == FunctionType::SPINOR) {
#if defined(_DEBUG_)
              if (!cint2e(buf, shls, atm_.data(), natm_, bas_.data(), nbas_,
                          env_.data(), opt_))
                throw std::runtime_error("cint2e failed");
#else
              cint2e(buf, shls, atm_.data(), natm_, bas_.data(), nbas_,
                              env_.data(), opt_));
#endif // _DEBUG_
            }
            // copy buf
            for (int ii = 0; ii < di; ++ii) {
              auto I = x + ii;
              for (int jj = 0; jj < dj; ++jj) {
                auto J = y + jj;
                for (int kk = 0; kk < dk; ++kk) {
                  auto K = z + kk;
                  for (int ll = 0; ll < dl; ++ll) {
                    auto L = w + ll;

                    T val;
#define IJKL(I, J, K, L)                                                       \
  ((I) * nao_ * nao_ * nao_ + (J) * nao_ * nao_ + (K) * nao_ + (L))
                    if constexpr (FT == FunctionType::SPINOR) {
                      val = std::complex<double>(
                          buf[2 * (ll * di * dj * dk + kk * di * dj + jj * di +
                                   ii)],
                          buf[2 * (ll * di * dj * dk + kk * di * dj + jj * di +
                                   ii) +
                              1]);
                    } else {
                      val =
                          buf[ll * di * dj * dk + kk * di * dj + jj * di + ii];
                    }

                    int2e_mat[IJKL(I, J, K, L)] = val;
                    int2e_mat[IJKL(I, J, L, K)] = val;
                    int2e_mat[IJKL(J, I, K, L)] = val;
                    int2e_mat[IJKL(J, I, L, K)] = val;
                    int2e_mat[IJKL(K, L, I, J)] = val;
                    int2e_mat[IJKL(K, L, J, I)] = val;
                    int2e_mat[IJKL(L, K, I, J)] = val;
                    int2e_mat[IJKL(L, K, J, I)] = val;
                  }
                }
              }
            }
            delete[] buf;
            delete[] shls;
          }
        }
      }
    }
    CINTdel_optimizer(&opt_);
    opt_ = nullptr;
    return int2e_mat;
  }
  T *int2e_s4() { // 4-fold symmetry two-electron integral in row_major format
    T *int2e_mat = new T[nao_ * nao_ * nao_ * nao_];
    memset(int2e_mat, 0, sizeof(T) * nao_ * nao_ * nao_ * nao_);
    if (opt_ != nullptr) {
      CINTdel_optimizer(&opt_);
      opt_ = nullptr;
    }
    if constexpr (FT == FunctionType::CART) {
      cint2e_cart_optimizer(&opt_, atm_.data(), natm_, bas_.data(), nbas_,
                            env_.data());
    } else if constexpr (FT == FunctionType::SPH) {
      cint2e_sph_optimizer(&opt_, atm_.data(), natm_, bas_.data(), nbas_,
                           env_.data());
    } else if constexpr (FT == FunctionType::SPINOR) {
      cint2e_optimizer(&opt_, atm_.data(), natm_, bas_.data(), nbas_,
                       env_.data());
    } // FunctionType

#pragma omp parallel for collapse(4)
    for (int i = 0; i < nbas_; ++i) {
      for (int j = 0; j <= i; ++j) {
        for (int k = 0; k < nbas_; ++k) {
          for (int l = 0; l <= i; ++l) {
            int *shls = new int[4]{i, j, k, l};
            int di, dj, dk, dl;
            int x, y, z, w;
            if constexpr (FT == FunctionType::CART) {
              di = CINTcgto_cart(i, bas_.data());
              dj = CINTcgto_cart(j, bas_.data());
              dk = CINTcgto_cart(k, bas_.data());
              dl = CINTcgto_cart(l, bas_.data());
              x = CINTtot_cgto_cart(bas_.data(), i);
              y = CINTtot_cgto_cart(bas_.data(), j);
              z = CINTtot_cgto_cart(bas_.data(), k);
              w = CINTtot_cgto_cart(bas_.data(), l);
            } else if constexpr (FT == FunctionType::SPH) {
              di = CINTcgto_spheric(i, bas_.data());
              dj = CINTcgto_spheric(j, bas_.data());
              dk = CINTcgto_spheric(k, bas_.data());
              dl = CINTcgto_spheric(l, bas_.data());
              x = CINTtot_cgto_spheric(bas_.data(), i);
              y = CINTtot_cgto_spheric(bas_.data(), j);
              z = CINTtot_cgto_spheric(bas_.data(), k);
              w = CINTtot_cgto_spheric(bas_.data(), l);
            } else if constexpr (FT == FunctionType::SPINOR) {
              di = CINTcgto_spinor(i, bas_.data());
              dj = CINTcgto_spinor(j, bas_.data());
              dk = CINTcgto_spinor(k, bas_.data());
              dl = CINTcgto_spinor(l, bas_.data());
              x = CINTtot_cgto_spinor(bas_.data(), i);
              y = CINTtot_cgto_spinor(bas_.data(), j);
              z = CINTtot_cgto_spinor(bas_.data(), k);
              w = CINTtot_cgto_spinor(bas_.data(), l);
            } // FunctionType
            double *buf;
            if constexpr (FT == FunctionType::SPINOR) {
              buf = new double[2 * di * dj * dk * dl];
            } else {
              buf = new double[di * dj * dk * dl];
            }
            if constexpr (FT == FunctionType::CART) {
#if defined(_DEBUG_)
              if (!cint2e_cart(buf, shls, atm_.data(), natm_, bas_.data(),
                               nbas_, env_.data(), opt_))
                throw std::runtime_error("cint2e_cart failed");
#else
              cint2e_cart(buf, shls, atm_.data(), natm_, bas_.data(), nbas_,
                          env_.data(), opt_));
#endif // _DEBUG_
            } else if constexpr (FT == FunctionType::SPH) {
#if defined(_DEBUG_)
              if (!cint2e_sph(buf, shls, atm_.data(), natm_, bas_.data(), nbas_,
                              env_.data(), opt_))
                throw std::runtime_error("cint2e_sph failed");
#else
              cint2e_sph(buf, shls, atm_.data(), natm_, bas_.data(), nbas_,
                              env_.data(), opt_));
#endif // _DEBUG_
            } else if constexpr (FT == FunctionType::SPINOR) {
#if defined(_DEBUG_)
              if (!cint2e(buf, shls, atm_.data(), natm_, bas_.data(), nbas_,
                          env_.data(), opt_))
                throw std::runtime_error("cint2e_spinor failed");
#else
              cint2e(buf, shls, atm_.data(), natm_, bas_.data(), nbas_,
                              env_.data(), opt_));
#endif // _DEBUG_
            } // FunctionType
            // copy buf to ovlp_mat
            for (int ii = 0; ii < di; ++ii) {
              auto I = x + ii;
              for (int jj = 0; jj < dj; ++jj) {
                auto J = y + jj;
                for (int kk = 0; kk < dk; ++kk) {
                  auto K = z + kk;
                  for (int ll = 0; ll < dl; ++ll) {
                    auto L = w + ll;

#define IJKL(I, J, K, L)                                                       \
  ((I) * nao_ * nao_ * nao_ + (J) * nao_ * nao_ + (K) * nao_ + (L))
                    T val;
                    if constexpr (FT == FunctionType::SPINOR) {
                      val = std::complex<double>(
                          buf[2 * (ll * di * dj * dk + kk * di * dj + jj * di +
                                   ii)],
                          buf[2 * (ll * di * dj * dk + kk * di * dj + jj * di +
                                   ii) +
                              1]);

                    } else {
                      val =
                          buf[ll * di * dj * dk + kk * di * dj + jj * di + ii];
                    }
                    int2e_mat[IJKL(I, J, K, L)] = val;
                    int2e_mat[IJKL(K, L, I, J)] = val;
                    if constexpr (std::is_same_v<T, std::complex<double>>) {
                      int2e_mat[IJKL(J, I, L, K)] = std::conj(val);
                      int2e_mat[IJKL(L, K, J, I)] = std::conj(val);
                    } else {
                      int2e_mat[IJKL(J, I, L, K)] = val;
                      int2e_mat[IJKL(L, K, J, I)] = val;
                    }
                  }
                }
              }
            }
            delete[] buf;
            delete[] shls;
          }
        }
      }
    }
    CINTdel_optimizer(&opt_);
    opt_ = nullptr;
    return int2e_mat;
  }

  T *int2e_s1() { // 1-fold symmetry two-electron integral in row_major format
    T *int2e_mat = new T[nao_ * nao_ * nao_ * nao_];
    memset(int2e_mat, 0, sizeof(T) * nao_ * nao_ * nao_ * nao_);
    if (opt_ != nullptr) {
      CINTdel_optimizer(&opt_);
      opt_ = nullptr;
    }
    if constexpr (FT == FunctionType::CART) {
      cint2e_cart_optimizer(&opt_, atm_.data(), natm_, bas_.data(), nbas_,
                            env_.data());
    } else if constexpr (FT == FunctionType::SPH) {
      cint2e_sph_optimizer(&opt_, atm_.data(), natm_, bas_.data(), nbas_,
                           env_.data());
    } else if constexpr (FT == FunctionType::SPINOR) {
      cint2e_optimizer(&opt_, atm_.data(), natm_, bas_.data(), nbas_,
                       env_.data());
    } // FunctionType

#pragma omp parallel for collapse(4)
    for (int i = 0; i < nbas_; ++i) {
      for (int j = 0; j < nbas_; ++j) {
        for (int k = 0; k < nbas_; ++k) {
          for (int l = 0; l < nbas_; ++l) {
            int *shls = new int[4]{i, j, k, l};
            int di, dj, dk, dl;
            int x, y, z, w;
            if constexpr (FT == FunctionType::CART) {
              di = CINTcgto_cart(i, bas_.data());
              dj = CINTcgto_cart(j, bas_.data());
              dk = CINTcgto_cart(k, bas_.data());
              dl = CINTcgto_cart(l, bas_.data());
              x = CINTtot_cgto_cart(bas_.data(), i);
              y = CINTtot_cgto_cart(bas_.data(), j);
              z = CINTtot_cgto_cart(bas_.data(), k);
              w = CINTtot_cgto_cart(bas_.data(), l);
            } else if constexpr (FT == FunctionType::SPH) {
              di = CINTcgto_spheric(i, bas_.data());
              dj = CINTcgto_spheric(j, bas_.data());
              dk = CINTcgto_spheric(k, bas_.data());
              dl = CINTcgto_spheric(l, bas_.data());
              x = CINTtot_cgto_spheric(bas_.data(), i);
              y = CINTtot_cgto_spheric(bas_.data(), j);
              z = CINTtot_cgto_spheric(bas_.data(), k);
              w = CINTtot_cgto_spheric(bas_.data(), l);
            } else if constexpr (FT == FunctionType::SPINOR) {
              di = CINTcgto_spinor(i, bas_.data());
              dj = CINTcgto_spinor(j, bas_.data());
              dk = CINTcgto_spinor(k, bas_.data());
              dl = CINTcgto_spinor(l, bas_.data());
              x = CINTtot_cgto_spinor(bas_.data(), i);
              y = CINTtot_cgto_spinor(bas_.data(), j);
              z = CINTtot_cgto_spinor(bas_.data(), k);
              w = CINTtot_cgto_spinor(bas_.data(), l);
            } // FunctionType
            double *buf;
            if constexpr (FT == FunctionType::SPINOR) {
              buf = new double[2 * di * dj * dk * dl];
            } else {
              buf = new double[di * dj * dk * dl];
            }
            if constexpr (FT == FunctionType::CART) {
#if defined(_DEBUG_)
              if (!cint2e_cart(buf, shls, atm_.data(), natm_, bas_.data(),
                               nbas_, env_.data(), opt_))
                throw std::runtime_error("cint2e_cart failed");
#else
              cint2e_cart(buf, shls, atm_.data(), natm_, bas_.data(), nbas_,
                          env_.data(), opt_));
#endif // _DEBUG_
            } else if constexpr (FT == FunctionType::SPH) {
#if defined(_DEBUG_)
              if (!cint2e_sph(buf, shls, atm_.data(), natm_, bas_.data(), nbas_,
                              env_.data(), opt_))
                throw std::runtime_error("cint2e_sph failed");
#else
              cint2e_sph(buf, shls, atm_.data(), natm_, bas_.data(), nbas_,
                              env_.data(), opt_));
#endif // _DEBUG_
            } else if constexpr (FT == FunctionType::SPINOR) {
#if defined(_DEBUG_)
              if (!cint2e(buf, shls, atm_.data(), natm_, bas_.data(), nbas_,
                          env_.data(), opt_))
                throw std::runtime_error("cint2e_spinor failed");
#else
              cint2e(buf, shls, atm_.data(), natm_, bas_.data(), nbas_,
                              env_.data(), opt_));
#endif // _DEBUG_
            } // FunctionType
            // copy buf to ovlp_mat
            for (int ii = 0; ii < di; ++ii) {
              auto I = x + ii;
              for (int jj = 0; jj < dj; ++jj) {
                auto J = y + jj;
                for (int kk = 0; kk < dk; ++kk) {
                  auto K = z + kk;
                  for (int ll = 0; ll < dl; ++ll) {
                    auto L = w + ll;

#define IJKL(I, J, K, L)                                                       \
  ((I) * nao_ * nao_ * nao_ + (J) * nao_ * nao_ + (K) * nao_ + (L))
                    T val;
                    if constexpr (FT == FunctionType::SPINOR) {
                      val = std::complex<double>(
                          buf[2 * (ll * di * dj * dk + kk * di * dj + jj * di +
                                   ii)],
                          buf[2 * (ll * di * dj * dk + kk * di * dj + jj * di +
                                   ii) +
                              1]);

                    } else {
                      val =
                          buf[ll * di * dj * dk + kk * di * dj + jj * di + ii];
                    }
                    int2e_mat[IJKL(I, J, K, L)] = val;
                  }
                }
              }
            }
            delete[] buf;
            delete[] shls;
          }
        }
      }
    }
    CINTdel_optimizer(&opt_);
    opt_ = nullptr;
    return int2e_mat;
  }
}; // Integral
} // namespace core

#endif // GMAT_INTEGRAL_HPP_
