#pragma once
#include <vector>
#include <cmath>
#include <memory>
#include <random>
#include <iostream>
#include <algorithm>
#include <bitset>

#include "vectypes.h"


#include "gmxcomplex.h"

#include <immintrin.h>

constexpr real ONE_OVER_SQRTPI = M_2_SQRTPI/2.0;
constexpr real TWO_PI = M_PI*2.0;

constexpr double c_electronCharge = 1.602176634e-19; /* Exact definition, Coulomb NIST 2018 CODATA */
constexpr double c_angs           = 1e-10;
constexpr double c_kilo           = 1e3;

constexpr double c_avogadro  = 6.02214076e23;     /* 1/mol, Exact definition, NIST 2018 CODATA */

constexpr double c_epsilon0Si = 8.8541878128e-12; /* F/m,  NIST 2018 CODATA */
constexpr double c_epsilon0 =
        ((c_epsilon0Si * c_angs * c_kilo) / (c_electronCharge * c_electronCharge * c_avogadro));
constexpr double c_one4PiEps0 = (1.0 / (4.0 * M_PI * c_epsilon0));



//! particle stores coordinates and its atom_type
struct Particle{
    gmx::RVec r;
    size_t atom_type;
};


static inline gmx::RVec deltarPBC(const gmx::RVec r1, const gmx::RVec r2,
                               const gmx::RVec& half_box_length, const gmx::RVec& box_length) {
    auto delta_r = r2 - r1;
    for (int dim = 0; dim < 3; ++dim) {
        if (delta_r[dim] > half_box_length[dim])
            delta_r[dim] = -(box_length[dim] - delta_r[dim]);
        else if( delta_r[dim] < -half_box_length[dim]){
            delta_r[dim] = box_length[dim] + delta_r[dim] ;
        }
    }
    return delta_r;
}

static inline real distancePBC(const real x_coord1, const real x_coord2,
                         const real& half_box_length, const real& box_length){
    real result = std::abs(x_coord2 - x_coord1) ;
    if( result < half_box_length )
        return  result ;
    else
        return box_length - result;
}

inline real dist_sq_pbc(const gmx::RVec & r1, const gmx::RVec & r2, const gmx::RVec& half_box_length,
                        const gmx::RVec& box_length ){
    auto dx = std::fabs(r1[0] - r2[0]);
    auto dy = std::fabs(r1[1] - r2[1]);
    auto dz = std::fabs(r1[2] - r2[2]);

    dx = (dx > half_box_length[XX])*(box_length[XX] - dx) + (dx < half_box_length[XX])*dx;
    dy = (dy > half_box_length[YY])*(box_length[YY] - dy) + (dy < half_box_length[YY])*dy;
    dz = (dz > half_box_length[ZZ])*(box_length[ZZ] - dz) + (dz < half_box_length[ZZ])*dz;

    return dx*dx + dy*dy + dz*dz;
}

inline real sphereInBox(const real half_box_length, const real r) {
    if( r < half_box_length)
        return 4. / 3. * M_PI * r * r * r;
    else if( r < std::sqrt(2.0) * half_box_length) {
        return 4. / 3. * M_PI * r * r * r -
               M_PI / 4. * (2. * r - 2. * half_box_length) * (2. * r - 2. * half_box_length)
               * (4. * r + 2. * half_box_length);
    }
    return 4. / 3. * M_PI * r * r * r;
}

struct Coordinates{
    std::vector<real> x;
    std::vector<real> z;
    std::vector<real> y;
};
class TrajectoryFrame {

public:
    Coordinates trajectory;
    size_t n_atoms = 0;
    size_t step = 0;
    size_t number = 0;
    float timestamp = 0;
    real bounding_box[3];
    std::vector<real> box;

    TrajectoryFrame() = default;

    TrajectoryFrame(size_t n_atoms) : n_atoms(n_atoms), box(3) {
        trajectory.x.resize(n_atoms);
        trajectory.y.resize(n_atoms);
        trajectory.z.resize(n_atoms);
    }

    void reset(const size_t n_atoms) {
        this->n_atoms = n_atoms;
        trajectory.x.resize(n_atoms);
        trajectory.y.resize(n_atoms);
        trajectory.z.resize(n_atoms);

        bounding_box[0] = 0.;
        bounding_box[1] = 0.;
        bounding_box[2] = 0.;
    }
};


static real min(real a, real b){
    return std::min(a, b);
}

