#pragma once

#include <algorithm>
#include <cstdint>
#include <sstream>
#include <string>
#include <vector>

#include <boost/graph/adjacency_list.hpp>
#include <boost/range/iterator_range.hpp>

namespace graph_hash {
    constexpr std::size_t WL_ITERATIONS = 4;
    constexpr std::uint64_t FNV_OFFSET = 1469598103934665603ULL;
    constexpr std::uint64_t FNV_PRIME = 1099511628211ULL;

    // https://en.wikipedia.org/wiki/Fowler%E2%80%93Noll%E2%80%93Vo_hash_function
    inline std::uint64_t fnv1a_mix(std::uint64_t h, std::uint64_t x) {
        for (int i = 0; i < 8; ++i) {
            h ^= static_cast<unsigned char>((x >> (8 * i)) & 0xffU);
            h *= FNV_PRIME;
        }
        return h;
    }

    inline std::uint64_t stable_hash_string(const std::string &s) {
        std::uint64_t h = FNV_OFFSET;
        for (unsigned char c : s) {
            h ^= c;
            h *= FNV_PRIME;
        }
        return h;
    }

    template <typename Graph>
    std::string vertex_label(const Graph &g, typename Graph::vertex_descriptor v) {
        std::ostringstream oss;
        oss << boost::get(boost::vertex_name, g, v);
        return oss.str();
    }

    template <typename Graph>
    std::vector<std::uint64_t> wl_vertex_colors(
        const Graph &g,
        std::size_t iterations = WL_ITERATIONS) {
        const std::size_t n = boost::num_vertices(g);
        auto index = boost::get(boost::vertex_index, g);

        std::vector<std::uint64_t> color(n);
        std::vector<std::uint64_t> next(n);

        for (auto v : boost::make_iterator_range(boost::vertices(g))) {
            color[index[v]] = stable_hash_string(vertex_label(g, v));
        }

        for (std::size_t it = 0; it < iterations; ++it) {
            for (auto v : boost::make_iterator_range(boost::vertices(g))) {
                std::vector<std::uint64_t> neigh;
                neigh.reserve(boost::out_degree(v, g));

                for (auto u : boost::make_iterator_range(boost::adjacent_vertices(v, g))) {
                    neigh.push_back(color[index[u]]);
                }
                std::sort(neigh.begin(), neigh.end());

                std::uint64_t h = FNV_OFFSET;
                h = fnv1a_mix(h, color[index[v]]);
                h = fnv1a_mix(h, static_cast<std::uint64_t>(neigh.size()));
                for (std::uint64_t c : neigh) {
                    h = fnv1a_mix(h, c);
                }
                next[index[v]] = h;
            }

            if (next == color) {
                break;
            }
            color.swap(next);
        }

        return color;
    }

    template <typename Graph>
    std::uint64_t wl_multiset_hash(const Graph &g, std::vector<std::uint64_t> colors) {
        std::sort(colors.begin(), colors.end());

        std::uint64_t h = FNV_OFFSET;
        h = fnv1a_mix(h, static_cast<std::uint64_t>(boost::num_vertices(g)));
        h = fnv1a_mix(h, static_cast<std::uint64_t>(boost::num_edges(g)));
        for (std::uint64_t c : colors) {
            h = fnv1a_mix(h, c);
        }
        return h;
    }

    template <typename Graph>
    std::uint64_t wl_hash(const Graph &g, std::size_t iterations = WL_ITERATIONS) {
        return wl_multiset_hash(g, wl_vertex_colors(g, iterations));
    }
}

// Convert a 64-bit value to a lowercase, zero-padded hex string.
// Always returns 16 hex characters.
inline std::string to_hex_string(std::uint64_t value)
{
    static constexpr char hex[] = "0123456789abcdef";

    std::string out(16, '0');
    for (int i = 15; i >= 0; --i) {
        out[i] = hex[value & 0xF];
        value >>= 4;
    }
    return out;
}
