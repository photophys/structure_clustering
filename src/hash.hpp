#pragma once

#include <string>
#include <cstdint>

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
