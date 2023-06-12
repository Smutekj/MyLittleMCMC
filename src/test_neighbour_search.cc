#include <gtest/gtest.h>
#include "Grid.h"
#include "VerletList.h"
#include "NeighbourSearcherImplementation.h"


TEST(TestSmallestDistance, BasicAssertions) {

    real box_size = 50;
    real r_cutoff = 10.;


    real cell_size = 10.0;
    int nx = box_size / cell_size;

    Grid grid(nx, box_size);

    real x = 25;
    int cell_x = 0;
    auto closest_dist_sq = closestDistanceSq1D(cell_size, box_size, cell_x, x);
    EXPECT_EQ(closest_dist_sq , 225);
    x = 5;
    cell_x = 4;
    closest_dist_sq = closestDistanceSq1D(cell_size, box_size, cell_x, x);
    EXPECT_EQ(closest_dist_sq , 25);
    x = 45;
    cell_x = 3;
    closest_dist_sq = closestDistanceSq1D(cell_size, box_size, cell_x, x);
    EXPECT_EQ(closest_dist_sq , 25);
    x = 5;
    cell_x = 0;
    closest_dist_sq = closestDistanceSq1D(cell_size, box_size, cell_x, x);
    EXPECT_EQ(closest_dist_sq , 0);

    gmx::RVec r = {23, 22, 21};
    cell_x = grid.cellCoordX(r);
    auto cell_index = grid.cellIndex(cell_x-1, cell_x-1, cell_x);
    closest_dist_sq = closestDistanceToCellSq(grid, cell_index, r);
    EXPECT_EQ(closest_dist_sq, 3*3 + 2*2);


    box_size = 55.;
    SearchGrid sgrid({box_size, box_size, box_size}, {cell_size, cell_size, cell_size});
    x = 51;
    cell_x = 0;
    closest_dist_sq = closestDistanceSq1D(cell_size, box_size, cell_x, x);
    EXPECT_EQ(closest_dist_sq , 4*4);

}

TEST(TestNearestCellSearch, BasicAssertions) {

    real box_size = 50;
    real r_cutoff = 10.;

    real cell_size = 10.0;
    int nx = box_size / cell_size;

    Grid grid(nx, box_size);
    SearchGrid sgrid({box_size, box_size, box_size}, {cell_size, cell_size, cell_size});

    gmx::RVec r = {25, 25, 25};
    auto cell_index = grid.coordToCell(r);
    std::vector<size_t> nearest_cell_indices = findNearestCellIndices(r, grid, r_cutoff);
    std::vector<size_t> nearest_cell_indices2 = findNearestCellIndices(cell_index, sgrid, r_cutoff);

    std::vector<size_t> correct_indices = grid.cell2nearest_neighbours.at(cell_index);

    std::sort(nearest_cell_indices.begin(), nearest_cell_indices.end());
    std::sort(correct_indices.begin(), correct_indices.end());
    std::sort(nearest_cell_indices2.begin(), nearest_cell_indices2.end());
    for (int i = 0; i < nearest_cell_indices.size(); ++i) {
        EXPECT_EQ(correct_indices.at(i), nearest_cell_indices.at(i));
        EXPECT_EQ(nearest_cell_indices2.at(i), nearest_cell_indices.at(i));
    }

    r = {19.5, 15, 15};
    cell_index = grid.coordToCell(r);

    nearest_cell_indices = findNearestCellIndices(r, grid, r_cutoff);
    correct_indices = findNearestCellIndicesStupid(r, grid, r_cutoff);

    std::sort(nearest_cell_indices.begin(), nearest_cell_indices.end());
    std::sort(correct_indices.begin(), correct_indices.end());
    for (int i = 0; i < correct_indices.size(); ++i) {
        EXPECT_EQ(correct_indices.at(i), nearest_cell_indices.at(i));
    }

    Grid grid2(nx*10, box_size);
    nearest_cell_indices = findNearestCellIndices(r, grid2, r_cutoff);
    correct_indices = findNearestCellIndicesStupid(r, grid2, r_cutoff);

    std::sort(nearest_cell_indices.begin(), nearest_cell_indices.end());
    std::sort(correct_indices.begin(), correct_indices.end());

    for (int i = 0; i < correct_indices.size(); ++i) {
        EXPECT_EQ(correct_indices.at(i), nearest_cell_indices.at(i));
    }
}