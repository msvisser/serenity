/*
 * Copyright (c) 2023, Michiel Visser <opensource@webmichiel.nl>
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#pragma once

#include <AK/ByteBuffer.h>
#include <AK/Endian.h>
#include <AK/MemoryStream.h>
#include <AK/Random.h>
#include <AK/StdLibExtras.h>
#include <AK/StringView.h>
#include <AK/UFixedBigInt.h>
#include <AK/UFixedBigIntDivision.h>
#include <LibCrypto/ASN1/DER.h>
#include <LibCrypto/Curves/EllipticCurve.h>

#include <LibCrypto/Curves/MontgomeryModularUFixedBigInt.h>

namespace Crypto::Curves {

struct SECPxxxr1CurveParameters {
    StringView prime;
    StringView a;
    StringView b;
    StringView order;
    StringView generator_point;
};

template<size_t bit_size, SECPxxxr1CurveParameters const& CURVE_PARAMETERS>
class SECPxxxr1 : public EllipticCurve {
private:
    using StorageType = AK::UFixedBigInt<bit_size>;
    using StorageTypeX2 = AK::UFixedBigInt<bit_size * 2>;

    // Curve parameters
    static constexpr size_t KEY_BIT_SIZE = bit_size;
    static constexpr size_t KEY_BYTE_SIZE = KEY_BIT_SIZE / 8;
    static constexpr size_t POINT_BYTE_SIZE = 1 + 2 * KEY_BYTE_SIZE;

    static constexpr StorageType make_unsigned_fixed_big_int_from_string(StringView str)
    {
        StorageType result { 0 };
        for (auto c : str) {
            if (c == '_')
                continue;

            result <<= 4;
            result |= parse_ascii_hex_digit(c);
        }
        return result;
    }

    static constexpr Array<u8, POINT_BYTE_SIZE> make_generator_point_bytes(StringView generator_point)
    {
        Array<u8, POINT_BYTE_SIZE> buf_array { 0 };

        auto it = generator_point.begin();
        for (size_t i = 0; i < POINT_BYTE_SIZE; i++) {
            if (it == CURVE_PARAMETERS.generator_point.end())
                break;

            while (*it == '_') {
                it++;
            }

            buf_array[i] = parse_ascii_hex_digit(*it) * 16;
            it++;
            if (it == CURVE_PARAMETERS.generator_point.end())
                break;

            buf_array[i] += parse_ascii_hex_digit(*it);
            it++;
        }

        return buf_array;
    }

    static constexpr StorageType PRIME = make_unsigned_fixed_big_int_from_string(CURVE_PARAMETERS.prime);
    static constexpr StorageType A = make_unsigned_fixed_big_int_from_string(CURVE_PARAMETERS.a);
    static constexpr StorageType B = make_unsigned_fixed_big_int_from_string(CURVE_PARAMETERS.b);
    static constexpr StorageType ORDER = make_unsigned_fixed_big_int_from_string(CURVE_PARAMETERS.order);
    static constexpr Array<u8, POINT_BYTE_SIZE> GENERATOR_POINT = make_generator_point_bytes(CURVE_PARAMETERS.generator_point);

    // Check that the generator point starts with 0x04
    static_assert(GENERATOR_POINT[0] == 0x04);

    // Verify that A = -3 mod p, which is required for some optimizations
    static_assert(A == PRIME - 3);

    using StorageTypeModPrime = MontgomeryModularUFixedBigInt<bit_size, PRIME>;
    using StorageTypeModOrder = MontgomeryModularUFixedBigInt<bit_size, ORDER>;

    struct JacobianPoint {
        StorageTypeModPrime x;
        StorageTypeModPrime y;
        StorageTypeModPrime z;
    };

public:
    size_t key_size() override { return POINT_BYTE_SIZE; }

    ErrorOr<ByteBuffer> generate_private_key() override
    {
        auto buffer = TRY(ByteBuffer::create_uninitialized(KEY_BYTE_SIZE));
        fill_with_random(buffer);
        return buffer;
    }

    ErrorOr<ByteBuffer> generate_public_key(ReadonlyBytes a) override
    {
        return compute_coordinate(a, GENERATOR_POINT);
    }

    ErrorOr<ByteBuffer> compute_coordinate(ReadonlyBytes scalar_bytes, ReadonlyBytes point_bytes) override
    {
        AK::FixedMemoryStream scalar_stream { scalar_bytes };
        AK::FixedMemoryStream point_stream { point_bytes };

        StorageType scalar = TRY(scalar_stream.read_value<BigEndian<StorageType>>());
        JacobianPoint point = TRY(read_uncompressed_point(point_stream));
        JacobianPoint result = TRY(compute_coordinate_internal(scalar, point));

        // Export the values into an output buffer
        auto buf = TRY(ByteBuffer::create_uninitialized(POINT_BYTE_SIZE));
        AK::FixedMemoryStream buf_stream { buf.bytes() };
        TRY(buf_stream.write_value<u8>(0x04));
        TRY(buf_stream.write_value<BigEndian<StorageType>>(static_cast<StorageType>(result.x)));
        TRY(buf_stream.write_value<BigEndian<StorageType>>(static_cast<StorageType>(result.y)));
        return buf;
    }

    ErrorOr<ByteBuffer> derive_premaster_key(ReadonlyBytes shared_point) override
    {
        VERIFY(shared_point.size() == POINT_BYTE_SIZE);
        VERIFY(shared_point[0] == 0x04);

        ByteBuffer premaster_key = TRY(ByteBuffer::create_uninitialized(KEY_BYTE_SIZE));
        premaster_key.overwrite(0, shared_point.data() + 1, KEY_BYTE_SIZE);
        return premaster_key;
    }

    ErrorOr<bool> verify(ReadonlyBytes hash, ReadonlyBytes pubkey, ReadonlyBytes signature)
    {
        Crypto::ASN1::Decoder asn1_decoder(signature);
        TRY(asn1_decoder.enter());

        auto r_bigint = TRY(asn1_decoder.read<Crypto::UnsignedBigInteger>(Crypto::ASN1::Class::Universal, Crypto::ASN1::Kind::Integer));
        auto s_bigint = TRY(asn1_decoder.read<Crypto::UnsignedBigInteger>(Crypto::ASN1::Class::Universal, Crypto::ASN1::Kind::Integer));

        StorageType r = 0u;
        StorageType s = 0u;
        for (size_t i = 0; i < (KEY_BIT_SIZE / 32); i++) {
            StorageType rr = r_bigint.words()[i];
            StorageType ss = s_bigint.words()[i];
            r |= (rr << (i * 32));
            s |= (ss << (i * 32));
        }

        StorageTypeModOrder r_mo = static_cast<StorageTypeModOrder>(r);
        StorageTypeModOrder s_mo = static_cast<StorageTypeModOrder>(s);

        // z is the hash
        AK::FixedMemoryStream hash_stream { hash };
        StorageType z = TRY(hash_stream.read_value<BigEndian<StorageType>>());
        StorageTypeModOrder z_mo = static_cast<StorageTypeModOrder>(z);

        AK::FixedMemoryStream pubkey_stream { pubkey };
        JacobianPoint pubkey_point = TRY(read_uncompressed_point(pubkey_stream));

        StorageTypeModOrder s_inv = s_mo.inverse();

        StorageTypeModOrder u1_mo = z_mo * s_inv;
        StorageTypeModOrder u2_mo = r_mo * s_inv;

        StorageType u1 = static_cast<StorageType>(u1_mo);
        StorageType u2 = static_cast<StorageType>(u2_mo);

        JacobianPoint point1 = TRY(generate_public_key_internal(u1));
        JacobianPoint point2 = TRY(compute_coordinate_internal(u2, pubkey_point));
        JacobianPoint result = point_add(point1, point2);

        // Convert from Jacobian coordinates back to Affine coordinates
        convert_jacobian_to_affine(result);

        // Make sure the resulting point is on the curve
        VERIFY(is_point_on_curve(result));

        return r.is_equal_to_constant_time(static_cast<StorageType>(result.x));
    }

private:
    ErrorOr<JacobianPoint> generate_public_key_internal(StorageType a)
    {
        AK::FixedMemoryStream generator_point_stream { GENERATOR_POINT };
        JacobianPoint point = TRY(read_uncompressed_point(generator_point_stream));
        return compute_coordinate_internal(a, point);
    }

    ErrorOr<JacobianPoint> compute_coordinate_internal(StorageType scalar, JacobianPoint point)
    {
        // FIXME: This will slightly bias the distribution of client secrets
        scalar = static_cast<StorageType>(static_cast<StorageTypeModOrder>(scalar));
        if (scalar.is_zero_constant_time())
            return Error::from_string_literal("SECPxxxr1: scalar is zero");

        // Check that the point is on the curve
        if (!is_point_on_curve(point))
            return Error::from_string_literal("SECPxxxr1: point is not on the curve");

        JacobianPoint result {};
        JacobianPoint temp_result {};

        // Calculate the scalar times point multiplication in constant time
        for (size_t i = 0; i < KEY_BIT_SIZE; i++) {
            temp_result = point_add(result, point);

            auto condition = (scalar & 1u) == 1u;
            result.x = StorageTypeModPrime::select(result.x, temp_result.x, condition);
            result.y = StorageTypeModPrime::select(result.y, temp_result.y, condition);
            result.z = StorageTypeModPrime::select(result.z, temp_result.z, condition);

            point = point_double(point);
            scalar >>= 1u;
        }

        // Convert from Jacobian coordinates back to Affine coordinates
        convert_jacobian_to_affine(result);

        // Make sure the resulting point is on the curve
        VERIFY(is_point_on_curve(result));

        return result;
    }

    ErrorOr<JacobianPoint> read_uncompressed_point(Stream& stream)
    {
        // Make sure the point is uncompressed
        if (TRY(stream.read_value<u8>()) != 0x04)
            return Error::from_string_literal("SECPxxxr1: point is not uncompressed format");

        JacobianPoint point {
            static_cast<StorageTypeModPrime>(TRY(stream.read_value<BigEndian<StorageType>>())),
            static_cast<StorageTypeModPrime>(TRY(stream.read_value<BigEndian<StorageType>>())),
            static_cast<StorageTypeModPrime>(1u),
        };

        return point;
    }

    JacobianPoint point_double(JacobianPoint const& point)
    {
        // Based on "Point Doubling" from http://point-at-infinity.org/ecc/Prime_Curve_Jacobian_Coordinates.html

        // if (Y == 0)
        //   return POINT_AT_INFINITY
        if (point.y.is_zero_constant_time()) {
            VERIFY_NOT_REACHED();
        }

        // Y2 = Y^2
        StorageTypeModPrime y2 = point.y.squared();

        // S = 4*X*Y2
        StorageTypeModPrime s = (point.x * y2).multiply_by_adding(4u);

        // M = 3*X^2 + a*Z^4 = 3*(X + Z^2)*(X - Z^2)
        // This specific equation from https://github.com/earlephilhower/bearssl-esp8266/blob/6105635531027f5b298aa656d44be2289b2d434f/src/ec/ec_p256_m64.c#L811-L816
        // This simplification only works because a = -3 mod p
        StorageTypeModPrime m = ((point.x + point.z.squared()) * (point.x - point.z.squared())).multiply_by_adding(3u);

        // X' = M^2 - 2*S
        StorageTypeModPrime xp = m.squared() - s.multiply_by_adding(2u);

        // Y' = M*(S - X') - 8*Y2^2
        StorageTypeModPrime yp = m * (s - xp) - y2.squared().multiply_by_adding(8u);

        // Z' = 2*Y*Z
        StorageTypeModPrime zp = (point.y * point.z).multiply_by_adding(2u);

        // return (X', Y', Z')
        return JacobianPoint { xp, yp, zp };
    }

    JacobianPoint point_add(JacobianPoint const& point_a, JacobianPoint const& point_b)
    {
        // Based on "Point Addition" from  http://point-at-infinity.org/ecc/Prime_Curve_Jacobian_Coordinates.html
        if (point_a.x.is_zero_constant_time() && point_a.y.is_zero_constant_time() && point_a.z.is_zero_constant_time()) {
            return point_b;
        }

        // U1 = X1*Z2^2
        StorageTypeModPrime u1 = point_a.x * point_b.z.squared();
        // S1 = Y1*Z2^3
        StorageTypeModPrime s1 = point_a.y * point_b.z.cubed();

        // U2 = X2*Z1^2
        StorageTypeModPrime u2 = point_b.x * point_a.z.squared();
        // S2 = Y2*Z1^3
        StorageTypeModPrime s2 = point_b.y * point_a.z.cubed();

        // if (U1 == U2)
        //   if (S1 != S2)
        //     return POINT_AT_INFINITY
        //   else
        //     return POINT_DOUBLE(X1, Y1, Z1)
        if (u1.is_equal_to_constant_time(u2)) {
            if (s1.is_equal_to_constant_time(s2)) {
                return point_double(point_a);
            } else {
                VERIFY_NOT_REACHED();
            }
        }

        // H = U2 - U1
        StorageTypeModPrime h = u2 - u1;
        StorageTypeModPrime h2 = h.squared();
        StorageTypeModPrime h3 = h.cubed();

        // R = S2 - S1
        StorageTypeModPrime r = s2 - s1;

        // X3 = R^2 - H^3 - 2*U1*H^2
        StorageTypeModPrime x3 = r.squared() - h3 - (u1 * h2).multiply_by_adding(2u);

        // Y3 = R*(U1*H^2 - X3) - S1*H^3
        StorageTypeModPrime y3 = r * ((u1 * h2) - x3) - (s1 * h3);

        // Z3 = H*Z1*Z2
        StorageTypeModPrime z3 = h * point_a.z * point_b.z;

        // return (X3, Y3, Z3)
        return JacobianPoint { x3, y3, z3 };
    }

    void convert_jacobian_to_affine(JacobianPoint& point)
    {
        // X' = X/Z^2
        point.x *= point.z.squared().inverse();
        // Y' = Y/Z^3
        point.y *= point.z.cubed().inverse();
        // Z' = 1
        point.z = StorageTypeModPrime(1u);
    }

    bool is_point_on_curve(JacobianPoint const& point)
    {
        // This check requires the point to be in Montgomery form, with Z=1
        // Calulcate Y^2 - X^3 - a*X - b = Y^2 - X^3 + 3*X - b
        StorageTypeModPrime result = point.y.squared();
        result -= point.x.cubed();
        result += point.x.multiply_by_adding(3u);
        result -= static_cast<StorageTypeModPrime>(B);

        return result.is_zero_constant_time() && point.z.is_equal_to_constant_time(static_cast<StorageTypeModPrime>(1u));
    }
};

// SECP256r1 curve
static constexpr SECPxxxr1CurveParameters SECP256r1_CURVE_PARAMETERS {
    .prime = "FFFFFFFF_00000001_00000000_00000000_00000000_FFFFFFFF_FFFFFFFF_FFFFFFFF"sv,
    .a = "FFFFFFFF_00000001_00000000_00000000_00000000_FFFFFFFF_FFFFFFFF_FFFFFFFC"sv,
    .b = "5AC635D8_AA3A93E7_B3EBBD55_769886BC_651D06B0_CC53B0F6_3BCE3C3E_27D2604B"sv,
    .order = "FFFFFFFF_00000000_FFFFFFFF_FFFFFFFF_BCE6FAAD_A7179E84_F3B9CAC2_FC632551"sv,
    .generator_point = "04_6B17D1F2_E12C4247_F8BCE6E5_63A440F2_77037D81_2DEB33A0_F4A13945_D898C296_4FE342E2_FE1A7F9B_8EE7EB4A_7C0F9E16_2BCE3357_6B315ECE_CBB64068_37BF51F5"sv,
};
using SECP256r1 = SECPxxxr1<256, SECP256r1_CURVE_PARAMETERS>;

// SECP384r1 curve
static constexpr SECPxxxr1CurveParameters SECP384r1_CURVE_PARAMETERS {
    .prime = "FFFFFFFF_FFFFFFFF_FFFFFFFF_FFFFFFFF_FFFFFFFF_FFFFFFFF_FFFFFFFF_FFFFFFFE_FFFFFFFF_00000000_00000000_FFFFFFFF"sv,
    .a = "FFFFFFFF_FFFFFFFF_FFFFFFFF_FFFFFFFF_FFFFFFFF_FFFFFFFF_FFFFFFFF_FFFFFFFE_FFFFFFFF_00000000_00000000_FFFFFFFC"sv,
    .b = "B3312FA7_E23EE7E4_988E056B_E3F82D19_181D9C6E_FE814112_0314088F_5013875A_C656398D_8A2ED19D_2A85C8ED_D3EC2AEF"sv,
    .order = "FFFFFFFF_FFFFFFFF_FFFFFFFF_FFFFFFFF_FFFFFFFF_FFFFFFFF_C7634D81_F4372DDF_581A0DB2_48B0A77A_ECEC196A_CCC52973"sv,
    .generator_point = "04_AA87CA22_BE8B0537_8EB1C71E_F320AD74_6E1D3B62_8BA79B98_59F741E0_82542A38_5502F25D_BF55296C_3A545E38_72760AB7_3617DE4A_96262C6F_5D9E98BF_9292DC29_F8F41DBD_289A147C_E9DA3113_B5F0B8C0_0A60B1CE_1D7E819D_7A431D7C_90EA0E5F"sv,
};
using SECP384r1 = SECPxxxr1<384, SECP384r1_CURVE_PARAMETERS>;

}
