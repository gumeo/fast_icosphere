/*
Copyright 2020 Gudmundur Einarsson

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated
documentation files (the "Software"), to deal in the Software without restriction, including without limitation the
rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit
persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the
Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE
WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#ifndef FAST_ICOSPHERE_H_
#define FAST_ICOSPHERE_H_


#include <algorithm>
#include <cmath>
#include <map>
#include <tuple>
#include <vector>

namespace icosphere {
/**
 * @brief Generate an icosphere with detailed control for the number of faces fast.
 *
 * @param num_additional Number of additional vertices on the edges of an original icosahedron.
 * @param project if true, project to sphere, otherwise only subdivide.
 */
std::pair<std::vector<std::array<double, 3>>, std::vector<std::vector<size_t>>> FastIcoSphere(const int num_additional,
                                                                                              const bool project);
} // namespace icosphere
#endif