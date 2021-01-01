#pragma once

// Underwater ray tracer
// Copyright 2020 Brede Johnsen
//
// This implementation is based on
// Ray Trace Modeling of Underwater Sound Propagation
// and
// PlaneRay: An acoustic underwater propagation model based 
// on ray tracingand plane-wave reflection coefficients
// both papers by Jens M. Hovem
// http://dx.doi.org/10.5772/55935

#include <algorithm>
#include <cassert>
#include <map>
#include <type_traits>
#include <vector>

namespace urt
{

/// @brief Linear interpolation
/// Available in std namespace since C++20
template <typename T, std::enable_if_t<std::is_floating_point_v<T>, bool> = true>
constexpr T lerp(T a, T b, T t)
{
    return a + t*(b-a);
}

// Exiting ray
template <typename T, std::enable_if_t<std::is_floating_point_v<T>, bool> = true>
struct basic_trace
{
    T dt;          // layer travel time [s]
    T dr;          // horizontal offset [m]
    T cos_th1;     // cosine of grazing angle [rad]
    bool turns;    // ray turns
};


template <typename T, std::enable_if_t<std::is_floating_point_v<T>, bool> = true>
struct basic_layer_boundary
{
    T z;           // depth [m]
    T c;           // sound speed [m/s]
};

template <typename T, std::enable_if_t<std::is_floating_point_v<T>, bool> = true>
struct basic_trace_pos
{
    T z;           // depth [m]
    T r;           // horizontal offset [m]
};

/// @brief Ray path as a collection of positions (z, r)
/// @tparam T 
template <typename T, std::enable_if_t<std::is_floating_point_v<T>, bool> = true>
using basic_ray_path = std::vector<basic_trace_pos<T>>;


/// @brief Sound speed c(z)
/// Eq. (16)
/// Within the layer the sound speed profile is approximated as linear
/// @param z0 depth [m]
/// @param c0 sound speed [m/s]
/// @param g gradient [1/s]
/// @param z depth [m]
/// @return sound speed
template <typename T, std::enable_if_t<std::is_floating_point_v<T>, bool> = true>
constexpr T sound_speed(T z0, T c0, T g, T z)
{
    return c0 + g * (z - z0);
}

/// @brief Sound speed gradient g (Eq. 4)
/// The positive or negative sign of the gradient determines whether
/// the sign of R is negative or positive, and thereby determines
/// if the ray path curves downward or upward.
/// @param z0 
/// @param c0 
/// @param z1 
/// @param c1 
/// @return gradient
template <typename T, std::enable_if_t<std::is_floating_point_v<T>, bool> = true>
constexpr T sound_speed_gradient(T z0, T c0, T z1, T c1)
{
    return (c1 - c0) / (z1 - z0);
}

/// @brief Ray parameter
/// (Greek letter Xi)
/// Eq. (1)
/// @param theta gracing angle of ray at given depth (Greek theta)
/// @param c sound speed
/// @return ray parameter
template <typename T, std::enable_if_t<std::is_floating_point_v<T>, bool> = true>
inline T ray_parameter(T theta, T c)
{
    return std::cos(theta) / c;
}

/// @brief Ray’s radius of curvature at given depth R(z)
/// Eq. (3)
/// @param xi ray parameter (Greek letter Xi)
/// @return radius
template <typename T, std::enable_if_t<std::is_floating_point_v<T>, bool> = true>
constexpr T ray_curvature_radius(T xi, T g)
{
    return -1.0 / (xi * g);
}

/// @brief Reflected gracing angle
/// Eq. 5
/// @param theta Gracing angle
/// @return angle out
template <typename T, std::enable_if_t<std::is_floating_point_v<T>, bool> = true>
constexpr T reflect(T theta)
{
    const T a = 0.0;       // angle of layer [rad]
    return theta + 2 * a;
}


/// @brief Trace ray through layer with linear sound speed
/// @tparam T 
template <typename T, std::enable_if_t<std::is_floating_point_v<T>, bool> = true>
struct basic_layer_tracer
{
    const T _z0, _z1;       // depth of layer boundaries
    const T _c0, _c1;       // sound speed at layer boundaries

    // Ray always enters at boundary 0 (implementation doesn't care what's up or down)
    basic_layer_tracer(basic_layer_boundary<T> b0, basic_layer_boundary<T> b1) :
        _z0(b0.z), _z1(b1.z), _c0(b0.c), _c1(b1.c)
    {
        assert(std::fabs(_z0-_z1) > 1.0e-8);
    }

    /// @brief Trace ray through layer
    /// Ray always enters at boundary 0
    /// @param cos_th0 Cosine of gracing angle to boundary layer
    /// @return ray travel
    basic_trace<T> trace(T cos_theta) const
    {
        assert(cos_theta <= 1);
        T g = (_c1 - _c0) / (_z1 - _z0);    // Sound speed gradient
        T cos_th0 = cos_theta;
        T xi = cos_th0 / _c0;               // Ray parameter Eq. (1)
        T cos_th1 = xi * _c1;               // From Eq. 1
        bool turns = false;                 // Does ray turn inside layer?
        if (cos_th1 >= 1.0)
        {
            turns = true;
            cos_th1 = cos_theta;  // if ray turns inside layer, it exits at same angle as it entered
        }
        T sin_th0 = std::sqrt(1 - cos_th0 * cos_th0);
        T sin_th1 = std::sqrt(1 - cos_th1 * cos_th1);

        // vertical ray
        if (std::fabs(cos_th0) < 1.0e-7)
        {
            if (std::fabs(g) < 1.0e-8)
            {
                T dz = std::fabs(_z1 - _z0);
                T dt = dz / _c1;
                T dr = 0;
                bool turns = false;
                assert(dt > 1.0e-8);
                return { dt, dr, cos_theta, turns };
            }
            else
            {
                T dt_log_arg = _c1 / _c0;
                if (dt_log_arg < 1)
                    dt_log_arg = 1 / dt_log_arg;      // swap direction
                assert(dt_log_arg >= 1);
                T dt = 1 / std::fabs(g) * std::log(dt_log_arg);
                T dr = 0;
                bool turns = false;
                assert(dt > 1.0e-8);
                return { dt, dr, cos_theta, turns };
            }
        }

        // constant sound speed
        if (std::fabs(g) < 1.0e-8)
        {
            // ray as a straight line
            // Eq. 11 & 12 in https://ffi-publikasjoner.archive.knowledgearc.net/bitstream/handle/20.500.12242/2128/08-00610.pdf
            T dz = std::fabs(_z1 - _z0);
            T dt = dz / (_c1 * sin_th1);
            T dr = dz * cos_th1 / sin_th1;
            bool turns = false;
            //assert(dr >= 0);
            assert(dt > 1.0e-8);
            return { dt, dr, cos_theta, turns };
        }

        // Ray going up
        if (_z0 > _z1)
        {
            sin_th0 = -sin_th0;
            sin_th1 = -sin_th1;
        }

        //const bool turns = cos_th1 >= 1.0;
        if (!turns)
        {
            // Travel time equation (Eq.18) only works for positive gradients.
            // If gradient is negative calculate travel time in opposite direction.
            T dt_log_arg = (_c1 / _c0) * (1 + sin_th0) / (1 + sin_th1);
            if (g < 0)
                dt_log_arg = 1 / dt_log_arg;                    // swap direction
            T dt = 1 / std::fabs(g) * std::log(dt_log_arg);     // ln(1/x) == -ln(x)
            assert(dt > 1.0e-8);
            T dr = 1 / (xi * g) * (sin_th0 - sin_th1);
            return { dt, dr, cos_th1, turns };
        }
        else // ray turns
        {
            T dt_log_arg = (1 + sin_th0) / cos_th0;
            if (g < 0)
                dt_log_arg = 1 / dt_log_arg;                    // swap direction
            T dt = 2 / std::fabs(g) * std::log(dt_log_arg);
            assert(dt > 1.0e-8);
            T dr = 2 / (xi * g) * sin_th0;
            return { dt, dr, cos_th1, turns };
        }
    }
};

// Sound speed profile
// key: depth [m]
// value: sound speed [m/s]
template <typename T, std::enable_if_t<std::is_floating_point_v<T>, bool> = true>
using sound_speed_profile = std::map<T, T>;


template <typename T, std::enable_if_t<std::is_floating_point_v<T>, bool> = true>
inline T sound_speed(const sound_speed_profile<T>& svp, T z)
{
    if (svp.empty())
        throw std::invalid_argument("Empty sound speed profile");

    if (z <= svp.begin()->first)   return svp.begin()->second;
    if (z >= svp.rbegin()->first)  return svp.rbegin()->second;

    auto last = svp.upper_bound(z);
    auto first = std::prev(last);
    auto z0 = first->first;
    auto z1 = last->first;
    auto c0 = first->second;
    auto c1 = last->second;
    auto t = (z - z0) / (z1 - z0);      // [0.0 - 1.0]
    return lerp(c0, c1, t);
}

/// @brief Find next depth that is a multiple of |dz|
/// @param z Start depth
/// @param dz Layer thickness (z1-z0).  Negative if ray direction is up.
/// @return Next depth
template <typename T, std::enable_if_t<std::is_floating_point_v<T>, bool> = true>
inline T next_depth(T z, T dz)
{
    if (dz > 0) // round up to multiple of dz
    {
        return std::ceil(std::nextafter(z / dz, std::numeric_limits<T>::max())) * dz;
    }
    else        // round down to multiple of |dz|
    {
        return std::floor(std::nextafter(z / -dz, std::numeric_limits<T>::lowest())) * -dz;
    }
}


/// @brief Ray tracer
/// @tparam T 
template <typename T, std::enable_if_t<std::is_floating_point_v<T>, bool> = true>
class ray_tracer
{
public:

    explicit ray_tracer(const sound_speed_profile<T>& svp) :
        _svp(svp)
    {
    }

    // Trace ray using sound speed profile
    // @param z Start depth for trace [m]
    // @param theta gracing angle [rad]
    // @param duration Time limit for trace [sec]
    // @return Ray path
    basic_ray_path<T> trace(T z, T theta, T duration, T dz = 50.0)
    {
        if (_svp.empty())
            throw std::out_of_range("Sound speed profile is empty");
        dz = std::copysign(dz, std::sin(theta));
        basic_ray_path<T> path;
        path.reserve(10000);
        basic_trace_pos<T> pos{ z, 0 };
        path.push_back(pos);

        auto z0 = z;
        auto c0 = sound_speed(_svp, z0);
        auto cos_th0 = std::cos(theta);
        while (duration > 0)
        {
            auto z1 = next_depth(z0,  dz);      // = z0 + dz;
            auto c1 = sound_speed(_svp, z1);

            basic_layer_tracer<T> layer_tracer(
                basic_layer_boundary<T>{ z0, c0 },
                basic_layer_boundary<T>{ z1, c1 });

            auto [dt, dr, cos_th1, turns] = layer_tracer.trace(cos_th0);
            if (dt <= 0) break;
            //assert(dr >= 0.0);
            assert(dt >= 0.0);

            if (turns)
            {
                dz = -dz;
                z1 = z0;
                c1 = c0;
            }

            pos.r += dr;
            pos.z = z1;
            path.push_back(pos);

            duration -= dt;
            cos_th0 = cos_th1;
            z0 = z1;
            c0 = c1;
        }
        return path;
    }

private:

    const sound_speed_profile<T>& _svp;
};


} // end namespace urt
