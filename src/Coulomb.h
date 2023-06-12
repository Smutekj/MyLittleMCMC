

#ifndef GR_GCMC_CLION_COULOMB_H
#define GR_GCMC_CLION_COULOMB_H

#include "real.h"
#include "core.h"
#include "Grid.h"
#include "Systemx.h"

namespace Electrostatics{

    real calc_ewaldcoeff_q(real rc, real rtol)
    {
        real beta = 5, low, high;
        int  n, i = 0;

        do
        {
            i++;
            beta *= 2;
        } while (std::erfc(beta * rc) > rtol);

        /* Do a binary search with tolerance 2^-60 */
        n    = i + 60;
        low  = 0;
        high = beta;
        for (i = 0; i < n; i++)
        {
            beta = (low + high) / 2;
            if (std::erfc(beta * rc) > rtol)
            {
                low = beta;
            }
            else
            {
                high = beta;
            }
        }
        return beta;
    }

# if __WORDSIZE == 64
#  define __PRI64_PREFIX	"l"
#  define __PRIPTR_PREFIX	"l"
# else
    #  define __PRI64_PREFIX	"ll"
#  define __PRIPTR_PREFIX
# endif
# define PRId64		__PRI64_PREFIX "d"


    void* save_calloc(const char* name, const char* file, int line, size_t nelem, size_t elsize)
    {
        void* p = nullptr;

        if ((nelem == 0) || (elsize == 0))
        {
            p = nullptr;
        }
        else
        {
#ifdef PRINT_ALLOC_KB
            if (nelem * elsize >= PRINT_ALLOC_KB * 1024)
        {
            int rank = gmx_node_rank();
            printf("Allocating %.1f MB for %s (called from file %s, line %d on %d)\n",
                   nelem * elsize / 1048576.0,
                   name,
                   file,
                   line,
                   rank);
        }
#endif
#if GMX_BROKEN_CALLOC
            /* emulate calloc(3) with malloc/memset on machines with
           a broken calloc, e.g. in -lgmalloc on cray xt3. */
        if ((p = malloc((size_t)nelem * (size_t)elsize)) == NULL)
        {
            gmx_fatal(errno,
                      __FILE__,
                      __LINE__,
                      "Not enough memory. Failed to calloc %" PRId64 " elements of size %" PRId64
                      " for %s\n(called from file %s, line %d)",
                      (int64_t)nelem,
                      (int64_t)elsize,
                      name,
                      file,
                      line);
        }
        memset(p, 0, (size_t)(nelem * elsize));
#else
            if ((p = calloc(nelem, elsize)) == nullptr)
            {
                std::runtime_error("neco je v pici s allocaci asi");
            }
#endif
        }
        return p;
    }

    template<typename T>
    static inline void gmx_snew_impl(const char* name, const char* file, int line, T*& ptr, size_t nelem)
    {
        // TODO: Use std::is_pod_v when CUDA 11 is a requirement.
        static_assert(std::is_pod<T>::value, "snew() called on C++ type");
        // NOLINTNEXTLINE bugprone-sizeof-expression
        ptr = static_cast<T*>(save_calloc(name, file, line, nelem, sizeof(T)));
    }
#define snew(ptr, nelem) gmx_snew_impl(#ptr, __FILE__, __LINE__, (ptr), (nelem))

    struct gmx_ewald_tab_t
    {
        gmx_ewald_tab_t(const int nx_, const int ny_, const int nz_) :
                nx(nx_+1), ny(ny_+1), nz(nz_+1), kmax(std::max(nx, std::max(ny, nz))){}


        int nx;
        int ny;
        int nz;
        int kmax;

        std::vector<t_complex> tab_xy;
        std::vector<t_complex> tab_qxyz;
    };

    //! Calculates wave vectors.
    static void calc_lll(const rvec box, rvec lll)
    {
        lll[XX] = 2.0 * M_PI / box[XX];
        lll[YY] = 2.0 * M_PI / box[YY];
        lll[ZZ] = 2.0 * M_PI / box[ZZ];
    }

    using cvec = std::array<t_complex, DIM>;

    //! Make tables for the structure factor parts
    static void tabulateStructureFactors(int natom, std::vector<gmx::RVec> r_coords, int kmax, cvec** eir, const rvec lll)
    {
        int i, j, m;

        if (kmax < 1)
        {
            printf("Go away! kmax = %d\n", kmax);
            exit(1);
        }

        for (i = 0; (i < natom); i++)
        {
            for (m = 0; (m < 3); m++)
            {
                eir[0][i][m].re = 1;
                eir[0][i][m].im = 0;
            }

            for (m = 0; (m < 3); m++)
            {
                eir[1][i][m].re = std::cos(r_coords[i][m] * lll[m]);
                eir[1][i][m].im = std::sin(r_coords[i][m] * lll[m]);
            }
            for (j = 2; (j < kmax); j++)
            {
                for (m = 0; (m < 3); m++)
                {
                    eir[j][i][m] = cmul(eir[j - 1][i][m], eir[1][i][m]);
                }
            }
        }
    }

    real directPart(const std::vector<CellVector<gmx::RVec>>& cell_list, const Systemx& system,
                    const real* charges,  const real box_length, const real beta, const real cutoff, const real epsilon_rel){

        real half_box_length = box_length/2.;
        real total_energy = 0.0;
        real scaleDirect = c_one4PiEps0 / epsilon_rel;


        for (int cell_index_n = 0; cell_index_n<system.cell_to_world.size(); ++cell_index_n) {
            const auto &r_cell_n = cell_list[cell_index_n].r;
            const auto& neighbouring_cells = system.grid->cell2nearest_neighbours[cell_index_n];
            const auto n_particles_in_cell_n = cell_list[cell_index_n].atoms_in_cell;

            for (const auto cell_index_m: neighbouring_cells) {
                const auto &r_cell_m = cell_list[cell_index_m].r;
                const auto n_particles_in_cell_m = cell_list[cell_index_m].atoms_in_cell;

                for (int i = 0; i < n_particles_in_cell_n; ++i) {
                    const auto world_index_i = system.cell_to_world[cell_index_n][i];

#pragma clang loop vectorize(enable)
                    for (int j = 0; j < n_particles_in_cell_m; ++j) {
                        const auto world_index_j = system.cell_to_world[cell_index_m][j];

                        const auto dx = distancePBC(r_cell_n[i][XX], r_cell_m[j][XX],
                                                    half_box_length, box_length) + box_length * (world_index_j == world_index_i);
                        const auto dy = distancePBC(r_cell_n[i][YY], r_cell_m[j][YY],
                                                    half_box_length, box_length);
                        const auto dz = distancePBC(r_cell_n[i][ZZ], r_cell_m[j][ZZ],
                                                    half_box_length, box_length);
                        const auto dr = std::sqrt(dx * dx + dy * dy + dz * dz);

                        total_energy += charges[world_index_i]*charges[world_index_j] *
                                        std::erfc(beta * dr)/dr ;//* (dr <= cutoff);
//                        std::cout << charges[world_index_i]*charges[world_index_j] * std::erfc(beta * dr)/dr << "\n";
                    }
                }
            }
        }
        return 0.5*total_energy*scaleDirect;
    }

     real calc_norm2_of_reciprocal_density(const std::vector<gmx::RVec>& r_coords, const real* charges, const gmx::RVec k){
        real re_reciprocal_density = 0;
        real im_reciprocal_density = 0;
         for(int i = 0; i < r_coords.size(); ++i){
             re_reciprocal_density += charges[i]*std::cos(gmx::dot(r_coords[i], k));
             im_reciprocal_density += charges[i]*std::sin(gmx::dot(r_coords[i], k));
         }
         return re_reciprocal_density*re_reciprocal_density + im_reciprocal_density*im_reciprocal_density;
    };

    real reciprocalPart(const Systemx* system,
                        const real* charges, const real box_size, const real beta, const real k_cutoff, const real epsilon_rel = 1){

        real total_energy(0.);
        const real factor = -1.0 / (4 * beta * beta);
        real scaleRecip = 4.0 * M_PI / (box_size*box_size*box_size) * c_one4PiEps0 / epsilon_rel;


        gmx::RVec box = {box_size, box_size, box_size};
        gmx::RVec lll;
        calc_lll(box, lll);
        const auto n_atoms = system->world_to_rcoord.size();

        const auto volume = box_size*box_size*box_size;

        const auto nx_min = - static_cast<int>(std::floor(box_size*k_cutoff/TWO_PI));
        const auto nx_max = static_cast<int>(std::floor(box_size*k_cutoff/TWO_PI));

        const auto ny_min = - static_cast<int>(std::floor(box_size*k_cutoff/TWO_PI));
        const auto ny_max = static_cast<int>(std::floor(box_size*k_cutoff/TWO_PI));

        const auto nz_min = - static_cast<int>(std::floor(box_size*k_cutoff/TWO_PI));
        const auto nz_max = static_cast<int>(std::floor(box_size*k_cutoff/TWO_PI));

        const auto& r_coords = system->world_to_rcoord;

        for( int nx = nx_min; nx<=nx_max; ++nx){
            for( int ny = ny_min;  ny<=ny_max; ++ny){
               for( int nz = nz_min; nz<=nz_max; ++nz){
                   const auto kx = nx*TWO_PI/box_size;
                   const auto ky = ny*TWO_PI/box_size;
                   const auto kz = nz*TWO_PI/box_size;
                   const auto k2 = kx*kx + ky*ky + kz*kz;

                   if(not(nx==0 and ny==0 and nz==0)) {
                       const auto rho_of_k_2 = calc_norm2_of_reciprocal_density(r_coords, charges, {kx, ky, kz});
                       total_energy += std::exp(k2 * factor) / k2 * rho_of_k_2;
                   }
               }
            }
        }
        return 0.5*total_energy*scaleRecip;
    }

    real reciprocalPartGmx(const Systemx* system,
                        const real* charges, const real box_size, const real beta, const real k_cutoff, const real epsilon_rel = 1) {

        gmx::RVec box = {box_size, box_size, box_size};
        gmx::RVec lll;
        cvec** eir;
        int ix, iy, iz;

        real total_energy(0.);
        const real factor = -1.0 / (4 * beta * beta);
        real scaleRecip = 4.0 * M_PI / (box_size*box_size*box_size) * c_one4PiEps0 / epsilon_rel;

        const int n_atoms = system->world_to_rcoord.size();

        const auto volume = box_size*box_size*box_size;

        const auto nx_max = static_cast<int>(std::floor(box_size*k_cutoff/TWO_PI)) ;

        const auto ny_max = static_cast<int>(std::floor(box_size*k_cutoff/TWO_PI));

        const auto nz_max = static_cast<int>(std::floor(box_size*k_cutoff/TWO_PI));

        auto et = std::make_unique<gmx_ewald_tab_t>(nx_max, ny_max, nz_max);

        const auto& r_coords = system->world_to_rcoord;
        int n;


        snew(eir, et->kmax);
        for (n = 0; n < et->kmax; n++) {
            snew(eir[n], n_atoms);
        }
        et->tab_xy.resize(n_atoms);
        et->tab_qxyz.resize(n_atoms);

        calc_lll(box, lll);

        tabulateStructureFactors(n_atoms, r_coords,  et->kmax, eir, lll);

        int lowiy        = 0;
        int lowiz        = 1;
        real energy_AB = 0;
        real m2, mx, my, mz, ak, cs, ss;
        for (ix = 0; ix < et->nx; ix++)
        {
            mx = ix * lll[XX];
            for (iy = lowiy; iy < et->ny; iy++)
            {
                my = iy * lll[YY];
                if (iy >= 0)
                {
                    for (n = 0; n < n_atoms; n++)
                    {
                        et->tab_xy[n] = cmul(eir[ix][n][XX], eir[iy][n][YY]);
                    }
                }
                else
                {
                    for (n = 0; n < n_atoms; n++)
                    {
                        et->tab_xy[n] = cmul(eir[ix][n][XX], conjugate(eir[-iy][n][YY]));
                    }
                }
                for (iz = lowiz; iz < et->nz; iz++)
                {
                    mz  = iz * lll[ZZ];
                    m2  = mx * mx + my * my + mz * mz;
                    ak  = std::exp(m2 * factor) / m2;
                    if (iz >= 0)
                    {
                        for (n = 0; n < n_atoms; n++)
                        {
                            et->tab_qxyz[n] = rcmul(charges[n], cmul(et->tab_xy[n], eir[iz][n][ZZ]));
                        }
                    }
                    else
                    {
                        for (n = 0; n < n_atoms; n++)
                        {
                            et->tab_qxyz[n] =
                                    rcmul(charges[n], cmul(et->tab_xy[n], conjugate(eir[-iz][n][ZZ])));
                        }
                    }

                    cs = ss = 0;
                    for (n = 0; n < n_atoms; n++)
                    {
                        cs += et->tab_qxyz[n].re;
                        ss += et->tab_qxyz[n].im;
                    }
                    energy_AB += ak * (cs * cs + ss * ss);

                    lowiz = 1 - et->nz;
                }
                lowiy = 1 - et->ny;
            }
        }

        energy_AB *= scaleRecip;
        return energy_AB;
    }

real selfPart(const real* charges, size_t n_particles, const real beta, const real epsilon_rel = 1.0){
        real total_energy(0.);
        for(int i = 0; i < n_particles; ++i){
            total_energy += - beta * ONE_OVER_SQRTPI * charges[i]*charges[i];
        }
        return  total_energy * c_one4PiEps0/epsilon_rel;
    }

    real dipolePart(const std::vector<gmx::RVec>& r_coords, const real* charges, size_t n_particles,
                    const real box_size, const real epsilon_rel=1.0){
        gmx::RVec a = {0, 0, 0};
        real volume = box_size*box_size*box_size;
        for(int i = 0; i < n_particles; ++i){
            a += charges[i]*r_coords[i];
        }
        return  TWO_PI/((1.+2.*epsilon_rel) * volume) * gmx::norm2(a);
    }

    real EwaldEnergy(const std::vector<CellVector<gmx::RVec>>& cell_list, const Systemx* system,
                     const real* charges, real beta , const real box_size, real r_cutoff, const real k_cutoff, const real epsilon_rel){

        beta = calc_ewaldcoeff_q(r_cutoff, 1e-5);
        const auto n_atoms = system->world_to_rcoord.size();

        real Q2 = 0;
        for(int i = 0; i< n_atoms; ++i){
            Q2 += charges[i]*charges[i];
        }
        const auto k_max = static_cast<int>(std::ceil(box_size*k_cutoff/TWO_PI)) ;
        const real dUr = c_one4PiEps0/epsilon_rel* Q2 * (r_cutoff/(2.*std::pow(box_size,3))) * std::pow((beta*r_cutoff), -2) * std::exp(-std::pow((beta*r_cutoff), 2));
        const real dUk = c_one4PiEps0/epsilon_rel* Q2 * beta/M_PI/M_PI * std::pow(k_max, -3./2.) * std::exp(-std::pow(M_PI*k_max/(beta*box_size), 2));
        std::cout << "beta = " << beta << " dUr = " <<  dUr << " dUk = " << dUk << "\n";

        const auto reciprocal_part = reciprocalPart(system, charges, box_size, beta, k_cutoff, epsilon_rel);
        const auto reciprocal_part_gmx = reciprocalPartGmx(system, charges, box_size, beta, k_cutoff, epsilon_rel);
        const auto direct_part = directPart(cell_list, *system, charges, box_size, beta, r_cutoff, epsilon_rel);
        const auto self_part = selfPart(charges, n_atoms, beta, epsilon_rel);
        const auto dipole_part = dipolePart(system->world_to_rcoord, charges, n_atoms, box_size, epsilon_rel);

        return self_part + direct_part + dipole_part + reciprocal_part;
    }

    real BasicEnergyLj(const Systemx* system,
                       const real* As, const real* Bs,const real box_size){
        real total_energy(0.);
        const real half_box_size = box_size/2.f;

        const auto& r_coords = system->world_to_rcoord;
        const auto n_atoms = r_coords.size();
        for (int i = 0; i < n_atoms; ++i) {
            for (int j = i + 1; j < n_atoms; ++j) {
                const auto dx = distancePBC(r_coords[i][XX], r_coords[j][XX],
                                            half_box_size, box_size);
                const auto dy = distancePBC(r_coords[i][YY], r_coords[j][YY],
                                            half_box_size, box_size);
                const auto dz = distancePBC(r_coords[i][ZZ], r_coords[j][ZZ],
                                            half_box_size, box_size);
                const auto dr = std::sqrt(dx * dx + dy * dy + dz * dz);
                total_energy +=  As[i] * std::pow(dr, -12) - Bs[i] * std::pow(dr, -6);
            }
        }
        return total_energy;
    }

    real BasicEnergy(const Systemx* system,
                     const real* charges,const real box_size, const real epsilon_rel){


        real total_energy(0.);
        const real half_box_size = box_size/2.f;

        const auto& r_coords = system->world_to_rcoord;
        const auto n_atoms = r_coords.size();

        const int nx_max = 0;
        const int ny_max = 0;
        const int nz_max = 0;
        const int nx_min = -0;
        const int ny_min = -0;
        const int nz_min = -0;


        for( int nx = nx_min; nx<=nx_max; ++nx){
            for( int ny = ny_min;  ny<=ny_max; ++ny){
                for( int nz = nz_min; nz<=nz_max; ++nz) {

                    if(nx == 0 and ny == 0 and nz == 0) {
                        for (int i = 0; i < n_atoms; ++i) {
                            for (int j = 0; j < n_atoms; ++j) {
                                if(i==j){
                                    continue;
                                }
                                const real dr = std::sqrt(dist_sq_pbc(r_coords[j], r_coords[i], half_box_size, box_size));
                                total_energy += charges[i] * charges[j] / (dr);
                            }
                        }
                    } else{
                        for (int i = 0; i < n_atoms; ++i) {
                            for (int j = 0; j < n_atoms; ++j) {

                                const auto dx = std::abs(r_coords[i][XX] - r_coords[j][XX] + nx*box_size);
                                const auto dy = std::abs(r_coords[i][YY] - r_coords[j][YY] + ny*box_size);
                                const auto dz = std::abs(r_coords[i][ZZ] - r_coords[j][ZZ] + nz*box_size);
                                const auto dr = std::sqrt(dx * dx + dy * dy + dz * dz);
                                total_energy += charges[i] * charges[j] / (dr);
                            }
                        }
                    }
                }
            }
        }

        return total_energy*0.5*c_one4PiEps0/epsilon_rel;
    }

    real BasicEnergyChange(const Systemx* system, const size_t changed_index, const gmx::RVec r_new,
                           const real* charges, const real beta, const real box_size){

        real changed_energy(0.);
        const real half_box_size = box_size/2.0;

        const auto& r_coords = system->world_to_rcoord;
        const auto n_atoms = r_coords.size();

        for( int i = 0; i < n_atoms; ++i) {
            if (i != changed_index) {
                const real r_ij_new = std::sqrt(dist_sq_pbc(r_new, r_coords[i], half_box_size, box_size));
                const real r_ij_old = std::sqrt(dist_sq_pbc(r_coords[changed_index], r_coords[i], half_box_size, box_size));
                changed_energy += charges[i] * charges[changed_index] * (1./ r_ij_new - 1./r_ij_old);
            }
        }
        return changed_energy;
    }



}

#endif //GR_GCMC_CLION_COULOMB_H

