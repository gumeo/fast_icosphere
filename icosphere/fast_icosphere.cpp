#include "fast_icosphere.h"

namespace {
/**
 * @brief Just fill the data for an icosahedron. Got from here:
 * https://observablehq.com/@mourner/fast-icosphere-mesh
 * Both inputs need to be empty, otherwise the face indices are not right.
 */
void BasicIcosahedron(std::vector<std::array<double, 3>>* vertexPositions,
                      std::vector<std::vector<size_t>>* faceIndices) {
  const double f = (1 + std::sqrt(5)) / 2;
  vertexPositions->push_back({-1, f, 0});
  vertexPositions->push_back({1, f, 0});
  vertexPositions->push_back({-1, -f, 0});
  vertexPositions->push_back({1, -f, 0});
  vertexPositions->push_back({0, -1, f});
  vertexPositions->push_back({0, 1, f});
  vertexPositions->push_back({0, -1, -f});
  vertexPositions->push_back({0, 1, -f});
  vertexPositions->push_back({f, 0, -1});
  vertexPositions->push_back({f, 0, 1});
  vertexPositions->push_back({-f, 0, -1});
  vertexPositions->push_back({-f, 0, 1});
  faceIndices->push_back({0, 11, 5});
  faceIndices->push_back({0, 5, 1});
  faceIndices->push_back({0, 1, 7});
  faceIndices->push_back({0, 7, 10});
  faceIndices->push_back({0, 10, 11});
  faceIndices->push_back({11, 10, 2});
  faceIndices->push_back({5, 11, 4});
  faceIndices->push_back({1, 5, 9});
  faceIndices->push_back({7, 1, 8});
  faceIndices->push_back({10, 7, 6});
  faceIndices->push_back({3, 9, 4});
  faceIndices->push_back({3, 4, 2});
  faceIndices->push_back({3, 2, 6});
  faceIndices->push_back({3, 6, 8});
  faceIndices->push_back({3, 8, 9});
  faceIndices->push_back({9, 8, 1});
  faceIndices->push_back({4, 9, 5});
  faceIndices->push_back({2, 4, 11});
  faceIndices->push_back({6, 2, 10});
  faceIndices->push_back({8, 6, 7});
}

void AddSingleEdge(
    const size_t first, const size_t second, const int num_additional,
    std::map<std::pair<size_t, size_t>, std::vector<std::pair<size_t, std::array<double, 3>>>>* edge_data,
    std::vector<std::array<double, 3>>* vertices) {
  if (edge_data->count(std::make_pair(first, second)) == 0) {
    std::pair<size_t, size_t> edge_key = std::make_pair(first, second);
    std::array<double, 3> first_pos = (*vertices)[first];
    std::array<double, 3> second_pos = (*vertices)[second];
    std::array<double, 3> jump_vec = {(second_pos[0] - first_pos[0]) / ((double)(num_additional + 1)),
                                      (second_pos[1] - first_pos[1]) / ((double)(num_additional + 1)),
                                      (second_pos[2] - first_pos[2]) / ((double)(num_additional + 1))};
    // Add first
    (*edge_data)[edge_key].push_back(std::make_pair(first, first_pos));
    // Add in between to vertices and map.
    size_t v_ind = vertices->size();
    for (size_t i = 0; i < (size_t)(num_additional); i++) {
      std::array<double, 3> pos = {first_pos[0] + jump_vec[0] * ((double)(i + 1)),
                                   first_pos[1] + jump_vec[1] * ((double)(i + 1)),
                                   first_pos[2] + jump_vec[2] * ((double)(i + 1))};
      vertices->push_back(pos);
      (*edge_data)[edge_key].push_back(std::make_pair(v_ind, pos));
      v_ind++;
    }
    // Add last, (already in vertices)
    (*edge_data)[edge_key].push_back(std::make_pair(second, second_pos));
  }
}

void AddIcoEdgePoints(
    const int num_additional, const std::vector<size_t>& inds,
    std::map<std::pair<size_t, size_t>, std::vector<std::pair<size_t, std::array<double, 3>>>>* edge_data,
    std::vector<std::array<double, 3>>* vertices) {
  std::vector<size_t> inds_s = inds;
  std::sort(inds_s.begin(), inds_s.end());
  AddSingleEdge(inds_s[0], inds_s[1], num_additional, edge_data, vertices);
  AddSingleEdge(inds_s[1], inds_s[2], num_additional, edge_data, vertices);
  AddSingleEdge(inds_s[0], inds_s[2], num_additional, edge_data, vertices);
}

void FillIcoFace(
    const std::vector<size_t>& inds,
    const std::map<std::pair<size_t, size_t>, std::vector<std::pair<size_t, std::array<double, 3>>>>& edge_data,
    std::vector<std::array<double, 3>>* vertices, std::vector<std::vector<size_t>>* faces) {
  // We fill each row of the tri
  std::vector<size_t> inds_s = inds;
  std::sort(inds_s.begin(), inds_s.end());
  size_t low = inds_s[0];
  size_t mid = inds_s[1];
  size_t high = inds_s[2];
  // Check if sort changes normal
  bool permutation_sign = false;
  if (inds[0] == low && inds[1] == mid) {
    permutation_sign = true;
  }
  if (inds[1] == low && inds[2] == mid) {
    permutation_sign = true;
  }
  if (inds[2] == low && inds[0] == mid) {
    permutation_sign = true;
  }
  const std::pair<size_t, size_t> bottom_key = std::make_pair(low, mid);
  const std::pair<size_t, size_t> row_key_left = std::make_pair(low, high);
  const std::pair<size_t, size_t> row_key_right = std::make_pair(mid, high);
  std::vector<std::pair<size_t, std::array<double, 3>>> bottom_row = edge_data.at(bottom_key);
  std::vector<std::pair<size_t, std::array<double, 3>>> left_side = edge_data.at(row_key_left);
  std::vector<std::pair<size_t, std::array<double, 3>>> right_side = edge_data.at(row_key_right);
  std::array<double, 3> jump = {bottom_row[1].second[0] - bottom_row[0].second[0],
                                bottom_row[1].second[1] - bottom_row[0].second[1],
                                bottom_row[1].second[2] - bottom_row[0].second[2]};
  for (size_t i = 1; i < left_side.size(); i++) {
    size_t v_ind = vertices->size();
    std::vector<std::pair<size_t, std::array<double, 3>>> new_bottom_row;
    new_bottom_row.push_back(left_side[i]);
    if (permutation_sign) {
      faces->push_back({bottom_row[0].first, left_side[i].first, bottom_row[1].first});
    } else {
      faces->push_back({bottom_row[0].first, bottom_row[1].first, left_side[i].first});
    }
    for (size_t j = 1; j < bottom_row.size() - 1; j++) {
      std::array<double, 3> pos;
      size_t new_tri_ind = v_ind;
      if (j == bottom_row.size() - 2) {
        // No need to create new point
        pos = right_side[i].second;
        new_tri_ind = right_side[i].first;
      } else {
        pos = {left_side[i].second[0] + (double)(j)*jump[0], left_side[i].second[1] + (double)(j)*jump[1],
               left_side[i].second[2] + (double)(j)*jump[2]};
        vertices->push_back(pos);
        v_ind++;
      }
      size_t ind_top_left = new_bottom_row[new_bottom_row.size() - 1].first;
      size_t ind_top_right = new_tri_ind;
      size_t ind_bottom_left = bottom_row[j].first;
      size_t ind_bottom_right = bottom_row[j + 1].first;
      if (permutation_sign) {
        faces->push_back({ind_top_left, ind_top_right, ind_bottom_left});
        faces->push_back({ind_bottom_left, ind_top_right, ind_bottom_right});
      } else {
        faces->push_back({ind_top_left, ind_bottom_left, ind_top_right});
        faces->push_back({ind_bottom_left, ind_bottom_right, ind_top_right});
      }
      new_bottom_row.push_back(std::make_pair(new_tri_ind, pos));
    }
    bottom_row = new_bottom_row;
  }
}

} // namespace

namespace icosphere {
std::pair<std::vector<std::array<double, 3>>, std::vector<std::vector<size_t>>> FastIcoSphere(const int num_additional,
                                                                                              const bool project) {
  // 1) Create basic icosahedron shape
  std::vector<std::array<double, 3>> vertexPositions;
  std::vector<std::vector<size_t>> faceIndices;
  std::vector<std::vector<size_t>> fineFaceIndices;
  BasicIcosahedron(&vertexPositions, &faceIndices);
  std::map<std::pair<size_t, size_t>, std::vector<std::pair<size_t, std::array<double, 3>>>> skeleton_points;
  for (size_t i = 0; i < faceIndices.size(); i++) {
    // 2) Add the "in between" vertices on the edges of the original, these hold the data to fill the rest of the
    // points.
    AddIcoEdgePoints(num_additional, faceIndices[i], &skeleton_points, &vertexPositions);
    // 3) Add the faces and rest of "internal" vertices, since the needed vertices have been defined
    FillIcoFace(faceIndices[i], skeleton_points, &vertexPositions, &fineFaceIndices);
  }

  // 4) Scale
  if (project) {
    for (size_t i = 0; i < vertexPositions.size(); i++) {
      std::array<double, 3> pos = vertexPositions[i];
      double pos_len = std::sqrt(pos[0] * pos[0] + pos[1] * pos[1] + pos[2] * pos[2]);
      pos[0] *= 1 / pos_len;
      pos[1] *= 1 / pos_len;
      pos[2] *= 1 / pos_len;
      vertexPositions[i] = pos;
    }
  }
  // 4) Flip normals out if that matters (reorder face inds).
  // TODO: Fix normals if needed.

  return std::make_pair(vertexPositions, fineFaceIndices);
}
} // namespace icosphere