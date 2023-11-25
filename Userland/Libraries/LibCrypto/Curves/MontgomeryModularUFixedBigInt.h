/*
 * Copyright (c) 2023, Michiel Visser <opensource@webmichiel.nl>
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#pragma once

#include <AK/UFixedBigInt.h>

namespace Crypto::Curves {

template<size_t bit_size, AK::UFixedBigInt<bit_size> const& MODULUS>
class MontgomeryModularUFixedBigInt {
private:
    using StorageType = AK::UFixedBigInt<bit_size>;
    using StorageTypeX2 = AK::UFixedBigInt<bit_size * 2>;

    static constexpr StorageType calculate_modular_inverse_mod_r(StorageType value)
    {
        // Calculate the modular multiplicative inverse of value mod 2^bit_size using the extended euclidean algorithm
        using StorageTypeP1 = AK::UFixedBigInt<bit_size + 1>;

        StorageTypeP1 old_r = value;
        StorageTypeP1 r = static_cast<StorageTypeP1>(1u) << bit_size;
        StorageTypeP1 old_s = 1u;
        StorageTypeP1 s = 0u;

        while (!r.is_zero_constant_time()) {
            StorageTypeP1 r_save = r;
            StorageTypeP1 quotient = old_r.div_mod(r, r);
            old_r = r_save;

            StorageTypeP1 s_save = s;
            s = old_s - quotient * s;
            old_s = s_save;
        }

        return static_cast<StorageType>(old_s);
    }

    static constexpr StorageType calculate_r2_mod(StorageType modulus)
    {
        // Calculate the value of R^2 mod modulus, where R = 2^bit_size
        using StorageTypeX2P1 = AK::UFixedBigInt<bit_size * 2 + 1>;

        StorageTypeX2P1 r2 = static_cast<StorageTypeX2P1>(1u) << (2 * bit_size);
        StorageTypeX2P1 result = r2 % modulus;
        return static_cast<StorageType>(result);
    }

public:
    static constexpr StorageType REDUCE_MODULUS = StorageType { 0 } - MODULUS;
    static constexpr StorageType MODULUS_INVERSE_MOD_R = StorageType { 0 } - calculate_modular_inverse_mod_r(MODULUS);
    static constexpr StorageType R2_MOD_MODULUS = calculate_r2_mod(MODULUS);

    constexpr MontgomeryModularUFixedBigInt()
        : m_value(to_montgomery(0u))
    {
    }

    explicit constexpr MontgomeryModularUFixedBigInt(u32 value)
        : m_value(to_montgomery(value))
    {
    }

    explicit constexpr MontgomeryModularUFixedBigInt(StorageType const& value)
        : m_value(to_montgomery(value))
    {
    }

    explicit constexpr operator StorageType() const
    {
        return from_montgomery(m_value);
    }

    constexpr StorageType const& value() const
    {
        return m_value;
    }

    constexpr MontgomeryModularUFixedBigInt& operator+=(MontgomeryModularUFixedBigInt const& other)
    {
        m_value = modular_add(m_value, other.m_value);
        return *this;
    }

    constexpr MontgomeryModularUFixedBigInt operator+(MontgomeryModularUFixedBigInt const& other) const
    {
        MontgomeryModularUFixedBigInt result = *this;
        result += other;
        return result;
    }

    constexpr MontgomeryModularUFixedBigInt& operator-=(MontgomeryModularUFixedBigInt const& other)
    {
        m_value = modular_sub(m_value, other.m_value);
        return *this;
    }

    constexpr MontgomeryModularUFixedBigInt operator-(MontgomeryModularUFixedBigInt const& other) const
    {
        MontgomeryModularUFixedBigInt result = *this;
        result -= other;
        return result;
    }

    constexpr MontgomeryModularUFixedBigInt& operator*=(MontgomeryModularUFixedBigInt const& other)
    {
        m_value = modular_multiply(m_value, other.m_value);
        return *this;
    }

    constexpr MontgomeryModularUFixedBigInt operator*(MontgomeryModularUFixedBigInt const& other) const
    {
        MontgomeryModularUFixedBigInt result = *this;
        result *= other;
        return result;
    }

    constexpr MontgomeryModularUFixedBigInt squared() const
    {
        return *this * *this;
    }

    constexpr MontgomeryModularUFixedBigInt cubed() const
    {
        return *this * *this * *this;
    }

    constexpr MontgomeryModularUFixedBigInt multiply_by_adding(u32 multiplicand) const
    {
        MontgomeryModularUFixedBigInt result { 0u };
        for (u32 i = 0; i < multiplicand; ++i) {
            result += *this;
        }
        return result;
    }

    constexpr MontgomeryModularUFixedBigInt inverse() const
    {
        MontgomeryModularUFixedBigInt result;
        result.m_value = modular_inverse(m_value);
        return result;
    }

    constexpr bool is_zero_constant_time() const
    {
        return m_value.is_zero_constant_time();
    }

    constexpr bool is_equal_to_constant_time(MontgomeryModularUFixedBigInt const& other) const
    {
        return m_value.is_equal_to_constant_time(other.m_value);
    }

    static constexpr MontgomeryModularUFixedBigInt select(MontgomeryModularUFixedBigInt const& left, MontgomeryModularUFixedBigInt const& right, bool condition)
    {
        MontgomeryModularUFixedBigInt result;
        result.m_value = select(left.m_value, right.m_value, condition);
        return result;
    }

private:
    static constexpr StorageType to_montgomery(StorageType const& value)
    {
        return modular_multiply(value, R2_MOD_MODULUS);
    }

    static constexpr StorageType from_montgomery(StorageType const& value)
    {
        return modular_multiply(value, 1u);
    }

    static constexpr StorageType select(StorageType const& left, StorageType const& right, bool condition)
    {
        // If condition = 0 return left else right
        StorageType mask = static_cast<StorageType>(condition) - 1;
        AK::taint_for_optimizer(mask);

        return (left & mask) | (right & ~mask);
    }

    static constexpr StorageType modular_add(StorageType const& left, StorageType const& right, bool carry_in = false)
    {
        bool carry = carry_in;
        StorageType output = left.addc(right, carry);

        // If there is a carry, subtract p by adding 2^KEY_BIT_SIZE - p
        StorageType addend = select(0u, REDUCE_MODULUS, carry);
        carry = false;
        output = output.addc(addend, carry);

        // If there is still a carry, subtract p by adding 2^KEY_BIT_SIZE - p
        addend = select(0u, REDUCE_MODULUS, carry);
        return output + addend;
    }

    static constexpr StorageType modular_sub(StorageType const& left, StorageType const& right)
    {
        bool borrow = false;
        StorageType output = left.subc(right, borrow);

        // If there is a borrow, add p by subtracting 2^KEY_BIT_SIZE - p
        StorageType sub = select(0u, REDUCE_MODULUS, borrow);
        borrow = false;
        output = output.subc(sub, borrow);

        // If there is still a borrow, add p by subtracting 2^KEY_BIT_SIZE - p
        sub = select(0u, REDUCE_MODULUS, borrow);
        return output - sub;
    }

    static constexpr StorageType modular_multiply(StorageType const& left, StorageType const& right)
    {
        // Modular multiplication using the Montgomery method: https://en.wikipedia.org/wiki/Montgomery_modular_multiplication
        // This requires that the inputs to this function are in Montgomery form.

        // T = left * right
        StorageTypeX2 mult = left.wide_multiply(right);
        StorageType mult_mod_r = static_cast<StorageType>(mult);

        // m = ((T mod R) * curve_p')
        StorageType m = mult_mod_r * MODULUS_INVERSE_MOD_R;

        // mp = (m mod R) * curve_p
        StorageTypeX2 mp = m.wide_multiply(MODULUS);

        // t = (T + mp)
        bool carry = false;
        mult_mod_r.addc(static_cast<StorageType>(mp), carry);

        // output = t / R
        StorageType mult_high = static_cast<StorageType>(mult >> bit_size);
        StorageType mp_high = static_cast<StorageType>(mp >> bit_size);
        return modular_add(mult_high, mp_high, carry);
    }

    static constexpr StorageType modular_inverse(StorageType const& value)
    {
        // Modular inverse modulo the curve prime can be computed using Fermat's little theorem: a^(p-2) mod p = a^-1 mod p.
        // Calculating a^(p-2) mod p can be done using the square-and-multiply exponentiation method, as p-2 is constant.
        StorageType base = value;
        StorageType result = to_montgomery(1u);
        StorageType prime_minus_2 = MODULUS - 2u;

        for (size_t i = 0; i < bit_size; i++) {
            if ((prime_minus_2 & 1u) == 1u) {
                result = modular_multiply(result, base);
            }
            base = modular_multiply(base, base);
            prime_minus_2 >>= 1u;
        }

        return result;
    }

    StorageType m_value;
};

}
