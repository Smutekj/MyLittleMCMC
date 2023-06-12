//
// Created by smutekj on 05.06.23.
//

#include "NeighbourSearcherImplementation.h"

#include "core.h"
#include "Systemx.h"
#include "Grid.h"

real closestDistanceSq1D(
        const real dx, const real box_size,
        const size_t nx, const real x){

    const real x_cell_left = nx*dx;
    const real x_cell_right = (nx+1)*dx;
    if(x < x_cell_right and x >= x_cell_left){
        return 0;
    }
    auto dx_closest = (x_cell_left - x)*(x<=x_cell_left) +  (x - x_cell_right)*(x>x_cell_left);
    const bool is_further_than_half_box_size =  dx_closest > box_size/2.0;
    dx_closest += (box_size - 2*dx_closest - dx)*(is_further_than_half_box_size);

    return dx_closest*dx_closest;
}

real closestDistanceSq1D( const real dx, const real box_size,
                          const int nx_total,
                          const int nxi, const int nxj){

    auto dnx = std::abs(nxi - nxj);
    dnx = std::min(nx_total - dnx, dnx);
    auto dx_closest =  dx * static_cast<real>(dnx - (1) + (nxi==nxj) );
    const bool is_further_than_half_box_size =  dx_closest >= box_size/2.0;
    dx_closest = dx_closest * (!is_further_than_half_box_size) + (box_size - dx_closest - dx)*(is_further_than_half_box_size);

    return dx_closest*dx_closest;
}

real closestDistanceToCellSq(const Grid& grid,
                             const size_t cell_index, const gmx::RVec r_coord){

    const auto nxi = grid.cellCoordX(cell_index);
    const auto nyi = grid.cellCoordY(cell_index);
    const auto nzi = grid.cellCoordZ(cell_index);

    const auto dxy      = grid.cell_sizexy;
    const auto dz       = grid.cell_size;
    const auto box_size = grid.box_size;

    return       closestDistanceSq1D( dxy, box_size[XX], nxi, r_coord[XX])
                 +  closestDistanceSq1D( dxy, box_size[YY], nyi, r_coord[YY])
                 +  closestDistanceSq1D( dz , box_size[ZZ], nzi, r_coord[ZZ]);
}

real closestDistanceSq(const SearchGrid& grid,
                       const size_t cell_index_i, const size_t cell_index_j){

    const auto nxi = grid.cellCoordX(cell_index_i);
    const auto nyi = grid.cellCoordY(cell_index_i);
    const auto nzi = grid.cellCoordZ(cell_index_i);
    const auto nxj = grid.cellCoordX(cell_index_j);
    const auto nyj = grid.cellCoordY(cell_index_j);
    const auto nzj = grid.cellCoordZ(cell_index_j);

    const auto dr      = grid.cell_size_;
    const auto box_size = grid.box_size_;


    return    closestDistanceSq1D(  dr[XX], box_size[XX], grid.n_cells_[XX], nxi, nxj)
              +  closestDistanceSq1D( dr[YY], box_size[YY], grid.n_cells_[YY], nyi, nyj)
              +  closestDistanceSq1D( dr[ZZ] , box_size[ZZ], grid.n_cells_[ZZ], nzi, nzj);
}


std::vector<size_t> findNearestCellIndices(const gmx::RVec r_coord, const Grid& grid, const real cutoff){
    std::vector<size_t> nearest_cell_indices;

    const real cutoff_sq = cutoff*cutoff;

    const size_t cell_index_center  = grid.coordToCell(r_coord);
    const int ix0 = static_cast<int>(grid.cellCoordX(r_coord));
    const int iy0 = static_cast<int>(grid.cellCoordY(r_coord));
    const int iz0 = static_cast<int>(grid.cellCoordZ(r_coord));

    int delta_ix_max = std::ceil(cutoff/grid.cell_sizexy);
    int delta_iy_max = std::ceil(cutoff/grid.cell_sizexy);
    const int delta_iz_max = std::ceil(cutoff/grid.cell_size);

    int delta_ix_min = -delta_ix_max;
    int delta_iy_min = -delta_iy_max;
    const int delta_iz_min = -delta_iz_max;

    const size_t nx = grid.nx;
    const size_t ny = grid.ny;
    const size_t nz = grid.nz;


    int delta_ix; int delta_iy; int delta_iz;
    for(delta_iz = delta_iz_min; delta_iz<=delta_iz_max; ++delta_iz) {
        delta_ix_max = std::ceil(cutoff/grid.cell_sizexy);
        delta_ix_min = -delta_ix_max;
        delta_iy_max = std::ceil(cutoff/grid.cell_sizexy);
        for(delta_iy = 0; delta_iy<=delta_iy_max; ++delta_iy) {

            for (delta_ix = 0; delta_ix <= delta_ix_max; ++delta_ix) {
                auto ix = delta_ix + ix0;
                auto iy = delta_iy + iy0;
                auto iz = delta_iz + iz0;
                ix += nx*(ix<0) - nx*(ix>=nx);
                iy += ny*(iy<0) - ny*(iy>=ny);
                iz += nz*(iz<0) - nz*(iz>=nz);
                assert(ix >= 0 and ix < nx and iy >= 0 and iy < ny and iz >= 0 and iz < nz);
                const auto cell_index_neighbour = grid.cellIndex(ix, iy, iz);
                const auto closest_dist_sq = closestDistanceToCellSq(grid, cell_index_neighbour, r_coord);
                if(closest_dist_sq <= cutoff_sq){
                    nearest_cell_indices.push_back(cell_index_neighbour);
                } else{
                    delta_ix_max = delta_ix;
                    break;
                }
            }
            if(delta_ix > 0) { //! look into left cells only if at least one cell was within cutoff
                for (delta_ix = -1; delta_ix >= delta_ix_min; --delta_ix) {
                    auto ix = delta_ix + ix0;
                    auto iy = delta_iy + iy0;
                    auto iz = delta_iz + iz0;
                    ix += nx * (ix < 0) - nx * (ix >= nx);
                    iy += ny * (iy < 0) - ny * (iy >= ny);
                    iz += nz * (iz < 0) - nz * (iz >= nz);
                    assert(ix >= 0 and ix < nx and iy >= 0 and iy < ny and iz >= 0 and iz < nz);
                    const auto cell_index_neighbour = grid.cellIndex(ix, iy, iz);
                    const auto closest_dist_sq = closestDistanceToCellSq(grid, cell_index_neighbour, r_coord);
                    if (closest_dist_sq <= cutoff_sq) {
                        nearest_cell_indices.push_back(cell_index_neighbour);
                    } else {
                        delta_ix_min = delta_ix;
                        break;
                    }
                }
            }
            else{ //! we move onto next delta_iy
                break;
            }
        }
    }

    for(delta_iz = delta_iz_min; delta_iz<=delta_iz_max; ++delta_iz) {
        delta_ix_max = std::ceil(cutoff/grid.cell_sizexy);
        delta_ix_min = -delta_ix_max;
        delta_iy_min = -std::ceil(cutoff/grid.cell_sizexy);
        for(delta_iy = -1; delta_iy >= delta_iy_min; --delta_iy) {

            for (delta_ix = 0; delta_ix <= delta_ix_max; ++delta_ix) {
                auto ix = delta_ix + ix0;
                auto iy = delta_iy + iy0;
                auto iz = delta_iz + iz0;
                ix += nx*(ix<0) - nx*(ix>=nx);
                iy += ny*(iy<0) - ny*(iy>=ny);
                iz += nz*(iz<0) - nz*(iz>=nz);
                assert(ix >= 0 and ix < nx and iy >= 0 and iy < ny and iz >= 0 and iz < nz);
                const auto cell_index_neighbour = grid.cellIndex(ix, iy, iz);
                const auto closest_dist_sq = closestDistanceToCellSq(grid, cell_index_neighbour, r_coord);
                if(closest_dist_sq <= cutoff_sq){
                    nearest_cell_indices.push_back(cell_index_neighbour);
                } else{
                    delta_ix_max = delta_ix;
                    break;
                }
            }
            if(delta_ix > 0) { //! look into left cells only if at least one cell was within cutoff
                for (delta_ix = -1; delta_ix >= delta_ix_min; --delta_ix) {
                    auto ix = delta_ix + ix0;
                    auto iy = delta_iy + iy0;
                    auto iz = delta_iz + iz0;
                    ix += nx * (ix < 0) - nx * (ix >= nx);
                    iy += ny * (iy < 0) - ny * (iy >= ny);
                    iz += nz * (iz < 0) - nz * (iz >= nz);
                    assert(ix >= 0 and ix < nx and iy >= 0 and iy < ny and iz >= 0 and iz < nz);
                    const auto cell_index_neighbour = grid.cellIndex(ix, iy, iz);
                    const auto closest_dist_sq = closestDistanceToCellSq(grid, cell_index_neighbour, r_coord);
                    if (closest_dist_sq <= cutoff_sq) {
                        nearest_cell_indices.push_back(cell_index_neighbour);
                    } else {
                        delta_ix_min = delta_ix;
                        break;
                    }
                }
            }
            else{
                break;
            }
        }
    }

    return nearest_cell_indices;
}


std::vector<size_t> findNearestCellIndices(const size_t cell_index, const SearchGrid& grid, const real cutoff){
    std::vector<size_t> nearest_cell_indices;

    const real cutoff_sq = cutoff*cutoff;

    const size_t cell_index_center  = cell_index;
    const int ix0 = static_cast<int>(grid.cellCoordX(cell_index_center));
    const int iy0 = static_cast<int>(grid.cellCoordY(cell_index_center));
    const int iz0 = static_cast<int>(grid.cellCoordZ(cell_index_center));

    int delta_ix_max = std::ceil(cutoff/grid.cell_size_[XX]);
    int delta_iy_max = std::ceil(cutoff/grid.cell_size_[YY]);
    const int delta_iz_max = std::ceil(cutoff/grid.cell_size_[ZZ]);

    int delta_ix_min = -delta_ix_max;
    int delta_iy_min = -delta_iy_max;
    const int delta_iz_min = -delta_iz_max;

    const auto nx = static_cast<int>(grid.n_cells_[XX]);
    const auto ny = static_cast<int>(grid.n_cells_[YY]);
    const auto nz = static_cast<int>(grid.n_cells_[ZZ]);


    int delta_ix; int delta_iy; int delta_iz;
    for(delta_iz = delta_iz_min; delta_iz<=delta_iz_max; ++delta_iz) {
        delta_ix_max = std::ceil(cutoff/grid.cell_size_[XX]);
        delta_ix_min = -delta_ix_max;
        delta_iy_max = std::ceil(cutoff/grid.cell_size_[YY]);
        for(delta_iy = 0; delta_iy<=delta_iy_max; ++delta_iy) {

            for (delta_ix = 0; delta_ix <= delta_ix_max; ++delta_ix) {
                int ix = delta_ix + ix0;
                int iy = delta_iy + iy0;
                int iz = delta_iz + iz0;
                ix += nx*(ix<0) - nx*(ix>=nx);
                iy += ny*(iy<0) - ny*(iy>=ny);
                iz += nz*(iz<0) - nz*(iz>=nz);

                assert(ix >= 0 and ix < nx and iy >= 0 and iy < ny and iz >= 0 and iz < nz);
                const auto cell_index_neighbour = grid.cellIndex(ix, iy, iz);
                const auto closest_dist_sq = closestDistanceSq(grid, cell_index_neighbour, cell_index_center);
                if(closest_dist_sq <= cutoff_sq){
                    nearest_cell_indices.push_back(cell_index_neighbour);
                } else{
                    delta_ix_max = delta_ix;
                    break;
                }
            }
            if(delta_ix > 0) { //! look into left cells only if at least one cell was within cutoff
                for (delta_ix = -1; delta_ix >= delta_ix_min; --delta_ix) {
                    auto ix = delta_ix + ix0;
                    auto iy = delta_iy + iy0;
                    auto iz = delta_iz + iz0;
                    ix += nx * (ix < 0) - nx * (ix >= nx);
                    iy += ny * (iy < 0) - ny * (iy >= ny);
                    iz += nz * (iz < 0) - nz * (iz >= nz);
                    assert(ix >= 0 and ix < nx and iy >= 0 and iy < ny and iz >= 0 and iz < nz);
                    const auto cell_index_neighbour = grid.cellIndex(ix, iy, iz);
                    const auto closest_dist_sq = closestDistanceSq(grid, cell_index_neighbour, cell_index_center);
                    if (closest_dist_sq <= cutoff_sq) {
                        nearest_cell_indices.push_back(cell_index_neighbour);
                    } else {
                        delta_ix_min = delta_ix;
                        break;
                    }
                }
            }
            else{ //! we move onto next delta_iy
                break;
            }
        }
    }

    for(delta_iz = delta_iz_min; delta_iz<=delta_iz_max; ++delta_iz) {
        delta_ix_max = std::ceil(cutoff/grid.cell_size_[XX]);
        delta_ix_min = -delta_ix_max;
        delta_iy_min = -std::ceil(cutoff/grid.cell_size_[YY]);
        for(delta_iy = -1; delta_iy >= delta_iy_min; --delta_iy) {

            for (delta_ix = 0; delta_ix <= delta_ix_max; ++delta_ix) {
                auto ix = delta_ix + ix0;
                auto iy = delta_iy + iy0;
                auto iz = delta_iz + iz0;
                ix += nx*(ix<0) - nx*(ix>=nx);
                iy += ny*(iy<0) - ny*(iy>=ny);
                iz += nz*(iz<0) - nz*(iz>=nz);
                assert(ix >= 0 and ix < nx and iy >= 0 and iy < ny and iz >= 0 and iz < nz);
                const auto cell_index_neighbour = grid.cellIndex(ix, iy, iz);
                const auto closest_dist_sq = closestDistanceSq(grid, cell_index_neighbour, cell_index_center);
                if(closest_dist_sq <= cutoff_sq){
                    nearest_cell_indices.push_back(cell_index_neighbour);
                } else{
                    delta_ix_max = delta_ix;
                    break;
                }
            }
            if(delta_ix > 0) { //! look into left cells only if at least one cell was within cutoff
                for (delta_ix = -1; delta_ix >= delta_ix_min; --delta_ix) {
                    auto ix = delta_ix + ix0;
                    auto iy = delta_iy + iy0;
                    auto iz = delta_iz + iz0;
                    ix += nx * (ix < 0) - nx * (ix >= nx);
                    iy += ny * (iy < 0) - ny * (iy >= ny);
                    iz += nz * (iz < 0) - nz * (iz >= nz);
                    assert(ix >= 0 and ix < nx and iy >= 0 and iy < ny and iz >= 0 and iz < nz);
                    const auto cell_index_neighbour = grid.cellIndex(ix, iy, iz);
                    const auto closest_dist_sq = closestDistanceSq(grid, cell_index_neighbour, cell_index_center);
                    if (closest_dist_sq <= cutoff_sq) {
                        nearest_cell_indices.push_back(cell_index_neighbour);
                    } else {
                        delta_ix_min = delta_ix;
                        break;
                    }
                }
            }
            else{
                break;
            }
        }
    }

    return nearest_cell_indices;
}


std::vector<size_t> findNearestCellIndicesStupid(const gmx::RVec r_coord, const Grid& grid, const real cutoff){

    std::vector<size_t> nearest_cell_indices;
    const real cutoff_sq = cutoff*cutoff;

    const int ix0 = static_cast<int>(grid.cellCoordX(r_coord));
    const int iy0 = static_cast<int>(grid.cellCoordY(r_coord));
    const int iz0 = static_cast<int>(grid.cellCoordZ(r_coord));

    const int delta_ix_max = std::ceil(cutoff/grid.cell_sizexy);
    const int delta_iy_max = std::ceil(cutoff/grid.cell_sizexy);
    const int delta_iz_max = std::ceil(cutoff/grid.cell_size);

    const int delta_ix_min = -delta_ix_max;
    const int delta_iy_min = -delta_iy_max;
    const int delta_iz_min = -delta_iz_max;

    const size_t nx = grid.nx;
    const size_t ny = grid.ny;
    const size_t nz = grid.nz;

    for( int dix = delta_ix_min; dix <= delta_ix_max; ++dix) {
        for( int diy = delta_iy_min; diy <= delta_iy_max; ++diy) {
            for (int diz = delta_iz_min; diz <= delta_iz_max; ++diz) {
                auto ix = dix + ix0;
                auto iy = diy + iy0;
                auto iz = diz + iz0;
                ix += nx*(ix<0) - nx*(ix>=nx);
                iy += ny*(iy<0) - ny*(iy>=ny);
                iz += nz*(iz<0) - nz*(iz>=nz);
                assert(ix >= 0 and ix < nx and iy >= 0 and iy < ny and iz >= 0 and iz < nz);
                const auto cell_index_neighbour = grid.cellIndex(ix, iy, iz);
                const auto closest_dist_sq = closestDistanceToCellSq(grid, cell_index_neighbour, r_coord);
                if(closest_dist_sq <= cutoff_sq){
                    nearest_cell_indices.push_back(cell_index_neighbour);
                }
            }
        }
    }
    return nearest_cell_indices;
}


std::vector<size_t> findNearestCellIndicesStupid(const size_t cell_index, const SearchGrid& grid, const real cutoff){

    std::vector<size_t> nearest_cell_indices;
    const real cutoff_sq = cutoff*cutoff;

    const auto cell_index_center = cell_index;
    const int ix0 = static_cast<int>(grid.cellCoordX(cell_index_center));
    const int iy0 = static_cast<int>(grid.cellCoordY(cell_index_center));
    const int iz0 = static_cast<int>(grid.cellCoordZ(cell_index_center));

    const int delta_ix_max = std::ceil(cutoff/grid.cell_size_[XX]);
    const int delta_iy_max = std::ceil(cutoff/grid.cell_size_[YY]);
    const int delta_iz_max = std::ceil(cutoff/grid.cell_size_[ZZ]);

    const int delta_ix_min = -delta_ix_max;
    const int delta_iy_min = -delta_iy_max;
    const int delta_iz_min = -delta_iz_max;

    const auto nx = static_cast<int>(grid.n_cells_[XX]);
    const auto ny = static_cast<int>(grid.n_cells_[YY]);
    const auto nz = static_cast<int>(grid.n_cells_[ZZ]);

    for( int dix = delta_ix_min; dix <= delta_ix_max; ++dix) {
        for( int diy = delta_iy_min; diy <= delta_iy_max; ++diy) {
            for (int diz = delta_iz_min; diz <= delta_iz_max; ++diz) {
                auto ix = dix + ix0;
                auto iy = diy + iy0;
                auto iz = diz + iz0;
                ix += nx*(ix<0) - nx*(ix>=nx);
                iy += ny*(iy<0) - ny*(iy>=ny);
                iz += nz*(iz<0) - nz*(iz>=nz);
                assert(ix >= 0 and ix < nx and iy >= 0 and iy < ny and iz >= 0 and iz < nz);
                const auto cell_index_neighbour = grid.cellIndex(ix, iy, iz);
                const auto closest_dist_sq = closestDistanceSq(grid, cell_index_neighbour, cell_index_center);
                if(closest_dist_sq <= cutoff_sq){
                    nearest_cell_indices.push_back(cell_index_neighbour);
                }
            }
        }
    }
    return nearest_cell_indices;
}