//
// Created by smutekj on 05.06.23.
//

#ifndef GR_GCMC_CLION_NEIGHBOURSEARCHERIMPLEMENTATION_H
#define GR_GCMC_CLION_NEIGHBOURSEARCHERIMPLEMENTATION_H

#include <vector>
#include "vectypes.h"

class Grid;
class SearchGrid;


real closestDistanceSq1D(
        real dx, real box_size,
        size_t nx, real x);
real closestDistanceSq1D( real dx, real box_size,
                          size_t nxi, size_t nxj);
real closestDistanceToCellSq(const Grid& grid,
                             size_t cell_index, gmx::RVec r_coord);
real closestDistanceSq(const SearchGrid& grid,
                       size_t cell_index_i, size_t cell_index_j);
std::vector<size_t> findNearestCellIndices(gmx::RVec r_coord, const Grid& grid, real cutoff);
std::vector<size_t> findNearestCellIndices(size_t cell_index, const SearchGrid& grid, real cutoff);
std::vector<size_t> findNearestCellIndicesStupid(const gmx::RVec r_coord, const Grid& grid, const real cutoff);
std::vector<size_t> findNearestCellIndicesStupid(const size_t cell_index, const SearchGrid& grid, const real cutoff);

class NeighbourSearcherImplementation {

public:
    NeighbourSearcherImplementation() = default;

protected:


};

class FullSearch : NeighbourSearcherImplementation {

};

#include <memory>
#include "Grid.h"

class GridSearch : NeighbourSearcherImplementation {

    std::shared_ptr<SearchGrid> search_grid_;

public:
    GridSearch(gmx::RVec box_size, gmx::RVec dr, const real cutoff, std::shared_ptr<SearchGrid> grid) : NeighbourSearcherImplementation()
    {
        search_grid_ = grid;
        updateNearestNeighbours(box_size, cutoff);
    }

    void updateNearestNeighbours(gmx::RVec box_size, const real cutoff){
        search_grid_->box_size_ = box_size;

        for(int i = 0; i<search_grid_->n_cells_[XX]; ++i){
            for(int j = 0; j<search_grid_->n_cells_[YY]; ++j) {
                for (int k = 0; k < search_grid_->n_cells_[ZZ]; ++k) {
                    const auto cell_index = search_grid_->cellIndex(i, j, k);
                    search_grid_->cell2nearest_neighbours.at(cell_index) = findNearestCellIndices(cell_index, *search_grid_, cutoff);
                }
            }
        }
    }

};



class KDTreeSearch : NeighbourSearcherImplementation {

};






#endif //GR_GCMC_CLION_NEIGHBOURSEARCHERIMPLEMENTATION_H
