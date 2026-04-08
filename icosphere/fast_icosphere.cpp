#include "fast_icosphere.h"

namespace {

struct EdgeKey {
  size_t first, second;
  bool operator==(const EdgeKey &o) const {
    return first == o.first && second == o.second;
  }
};

struct EdgeKeyHash {
  size_t operator()(const EdgeKey &k) const {
    // Szudzik pairing — injective for our index range
    return k.first >= k.second
               ? k.first * k.first + k.first + k.second
               : k.second * k.second + k.first;
  }
};

// Edge data stores the vertex index of the first endpoint, then intermediate
// vertices, then the last endpoint — exactly (num_additional + 2) entries.
using EdgeMap = std::unordered_map<EdgeKey, std::vector<size_t>, EdgeKeyHash>;

void BasicIcosahedron(std::vector<std::array<double, 3>> *vertexPositions,
                      std::vector<std::array<size_t, 3>> *faceIndices) {
  const double f = (1 + std::sqrt(5)) / 2;
  (*vertexPositions)[0] = {-1, f, 0};
  (*vertexPositions)[1] = {1, f, 0};
  (*vertexPositions)[2] = {-1, -f, 0};
  (*vertexPositions)[3] = {1, -f, 0};
  (*vertexPositions)[4] = {0, -1, f};
  (*vertexPositions)[5] = {0, 1, f};
  (*vertexPositions)[6] = {0, -1, -f};
  (*vertexPositions)[7] = {0, 1, -f};
  (*vertexPositions)[8] = {f, 0, -1};
  (*vertexPositions)[9] = {f, 0, 1};
  (*vertexPositions)[10] = {-f, 0, -1};
  (*vertexPositions)[11] = {-f, 0, 1};
  (*faceIndices)[0] = {0, 11, 5};
  (*faceIndices)[1] = {0, 5, 1};
  (*faceIndices)[2] = {0, 1, 7};
  (*faceIndices)[3] = {0, 7, 10};
  (*faceIndices)[4] = {0, 10, 11};
  (*faceIndices)[5] = {11, 10, 2};
  (*faceIndices)[6] = {5, 11, 4};
  (*faceIndices)[7] = {1, 5, 9};
  (*faceIndices)[8] = {7, 1, 8};
  (*faceIndices)[9] = {10, 7, 6};
  (*faceIndices)[10] = {3, 9, 4};
  (*faceIndices)[11] = {3, 4, 2};
  (*faceIndices)[12] = {3, 2, 6};
  (*faceIndices)[13] = {3, 6, 8};
  (*faceIndices)[14] = {3, 8, 9};
  (*faceIndices)[15] = {9, 8, 1};
  (*faceIndices)[16] = {4, 9, 5};
  (*faceIndices)[17] = {2, 4, 11};
  (*faceIndices)[18] = {6, 2, 10};
  (*faceIndices)[19] = {8, 6, 7};
}

inline void Sort3(size_t &a, size_t &b, size_t &c) {
  if (a > b) std::swap(a, b);
  if (b > c) std::swap(b, c);
  if (a > b) std::swap(a, b);
}

void AddSingleEdge(const size_t first, const size_t second,
                   const int num_additional, EdgeMap *edge_data,
                   std::vector<std::array<double, 3>> *vertices,
                   size_t *v_ind) {
  EdgeKey key{first, second};
  if (edge_data->count(key))
    return;

  const int total = num_additional + 2;
  std::vector<size_t> indices;
  indices.reserve(total);

  const auto &first_pos = (*vertices)[first];
  const auto &second_pos = (*vertices)[second];
  const double inv = 1.0 / (num_additional + 1);
  const double dx = (second_pos[0] - first_pos[0]) * inv;
  const double dy = (second_pos[1] - first_pos[1]) * inv;
  const double dz = (second_pos[2] - first_pos[2]) * inv;

  indices.push_back(first);
  for (int i = 0; i < num_additional; i++) {
    const double t = (double)(i + 1);
    (*vertices)[*v_ind] = {first_pos[0] + dx * t, first_pos[1] + dy * t,
                           first_pos[2] + dz * t};
    indices.push_back((*v_ind)++);
  }
  indices.push_back(second);

  edge_data->emplace(key, std::move(indices));
}

void AddIcoEdgePoints(const int num_additional, size_t i0, size_t i1,
                      size_t i2, EdgeMap *edge_data,
                      std::vector<std::array<double, 3>> *vertices,
                      size_t *v_ind) {
  Sort3(i0, i1, i2);
  AddSingleEdge(i0, i1, num_additional, edge_data, vertices, v_ind);
  AddSingleEdge(i1, i2, num_additional, edge_data, vertices, v_ind);
  AddSingleEdge(i0, i2, num_additional, edge_data, vertices, v_ind);
}

void FillIcoFace(const size_t orig0, const size_t orig1, const size_t orig2,
                 const EdgeMap &edge_data,
                 std::vector<std::array<double, 3>> *vertices, size_t *v_ind,
                 std::vector<std::array<size_t, 3>> *faces, size_t *face_ind) {
  size_t low = orig0, mid = orig1, high = orig2;
  Sort3(low, mid, high);

  // Determine if the sort is an even permutation of the original
  bool permutation_sign = false;
  if (orig0 == low && orig1 == mid)
    permutation_sign = true;
  if (orig1 == low && orig2 == mid)
    permutation_sign = true;
  if (orig2 == low && orig0 == mid)
    permutation_sign = true;

  const auto &bottom_inds = edge_data.at(EdgeKey{low, mid});
  const auto &left_inds = edge_data.at(EdgeKey{low, high});
  const auto &right_inds = edge_data.at(EdgeKey{mid, high});

  // Compute the jump vector from the bottom edge
  const auto &p0 = (*vertices)[bottom_inds[0]];
  const auto &p1 = (*vertices)[bottom_inds[1]];
  const double jx = p1[0] - p0[0];
  const double jy = p1[1] - p0[1];
  const double jz = p1[2] - p0[2];

  // We build rows bottom-to-top. Store only the index list for the current
  // bottom row — avoid storing positions redundantly.
  const size_t edge_len = bottom_inds.size();
  // Use a flat buffer for the current bottom row indices (avoids
  // vector-of-pair copies)
  std::vector<size_t> bot_row(bottom_inds.begin(), bottom_inds.end());
  std::vector<size_t> new_row;
  new_row.reserve(edge_len);

  for (size_t i = 1; i < left_inds.size(); i++) {
    new_row.clear();
    const size_t left_vi = left_inds[i];
    new_row.push_back(left_vi);

    // First triangle of the row
    if (permutation_sign) {
      (*faces)[(*face_ind)++] = {bot_row[0], left_vi, bot_row[1]};
    } else {
      (*faces)[(*face_ind)++] = {bot_row[0], bot_row[1], left_vi};
    }

    const auto &left_pos = (*vertices)[left_vi];
    for (size_t j = 1; j < bot_row.size() - 1; j++) {
      size_t new_vi;
      if (j != bot_row.size() - 2) {
        const double t = (double)j;
        (*vertices)[*v_ind] = {left_pos[0] + t * jx, left_pos[1] + t * jy,
                               left_pos[2] + t * jz};
        new_vi = (*v_ind)++;
      } else {
        new_vi = right_inds[i];
      }

      const size_t tl = new_row.back();
      const size_t tr = new_vi;
      const size_t bl = bot_row[j];
      const size_t br = bot_row[j + 1];
      if (permutation_sign) {
        (*faces)[(*face_ind)++] = {tl, tr, bl};
        (*faces)[(*face_ind)++] = {bl, tr, br};
      } else {
        (*faces)[(*face_ind)++] = {tl, bl, tr};
        (*faces)[(*face_ind)++] = {bl, br, tr};
      }
      new_row.push_back(new_vi);
    }
    std::swap(bot_row, new_row);
  }
}

size_t NumInternalVertices(int num_additional) {
  if (num_additional < 2) {
    return 0;
  }
  int max_num = num_additional - 1;
  return max_num * (max_num + 1) / 2;
}

} // namespace

namespace icosphere {
std::pair<std::vector<std::array<double, 3>>, std::vector<std::array<size_t, 3>>>
FastIcoSphere(const int num_additional, const bool project) {
  size_t num_vertices =
      12 + 30 * num_additional + 20 * NumInternalVertices(num_additional);
  std::vector<std::array<double, 3>> vertexPositions(num_vertices);
  std::vector<std::array<size_t, 3>> faceIndices(20);
  size_t num_faces = (size_t)(num_additional + 1) * (num_additional + 1) * 20;
  std::vector<std::array<size_t, 3>> fineFaceIndices(num_faces);
  size_t face_counter = 0;
  size_t vertex_counter = 12;
  BasicIcosahedron(&vertexPositions, &faceIndices);

  EdgeMap skeleton_points;
  skeleton_points.reserve(64);

  for (size_t i = 0; i < 20; i++) {
    const auto &fi = faceIndices[i];
    AddIcoEdgePoints(num_additional, fi[0], fi[1], fi[2], &skeleton_points,
                     &vertexPositions, &vertex_counter);
    FillIcoFace(fi[0], fi[1], fi[2], skeleton_points, &vertexPositions,
                &vertex_counter, &fineFaceIndices, &face_counter);
  }

  if (project) {
    for (size_t i = 0; i < vertexPositions.size(); i++) {
      auto &p = vertexPositions[i];
      const double inv_len =
          1.0 / std::sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]);
      p[0] *= inv_len;
      p[1] *= inv_len;
      p[2] *= inv_len;
    }
  }
  return std::make_pair(std::move(vertexPositions), std::move(fineFaceIndices));
}
} // namespace icosphere
