#pragma once

// Underwater ray tracer
// Copyright 2020 Brede Johnsen
//
// This implementation is based on
// Ray Trace Modeling of Underwater Sound Propagation
// by Jens M. Hovem
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
T lerp(T a, T b, T t)
{
    return a + t*(b-a);
}

// Exiting ray
struct trace_t
{
    double dt;          // layer travel time [s]
    double dr;          // Horizontal offset [m]
    double cos_th1;     // cosine of grazing angle [rad]
    bool turns;         // ray turns
};


struct layer_boundary
{
    double z;           // next_depth [m]
    double c;           // sound speed [m/s]
};

struct trace_pos
{
    double z;           // next_depth [m]
    double r;           // Horizontal offset [m]
};

using ray_path_t = std::vector<trace_pos>;


/// @brief Sound speed c(z)
/// Eq. (16)
/// Within the layer the sound speed profile is approximated as linear
/// @param z0 next_depth [m]
/// @param c0 sound speed [m/s]
/// @param g gradient [1/s]
/// @param z next_depth [m]
/// @return sound speed
constexpr double sound_speed(double z0, double c0, double g, double z)
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
constexpr double sound_speed_gradient(double z0, double c0, double z1, double c1)
{
    return (c1 - c0) / (z1 - z0);
}

/////////////////////////////////////
//

/// @brief Ray parameter
/// (Greek letter Xi)
/// Eq. (1)
/// @param theta gracing angle of ray at given next_depth (Greek theta)
/// @param c sound speed
/// @return ray parameter
inline double ray_parameter(double theta, double c)
{
    return std::cos(theta) / c;
}

/// @brief Ray’s radius of curvature at given next_depth R(z)
/// Eq. (3)
/// @param xi ray parameter (Greek letter Xi)
/// @return radius
constexpr double ray_curvature_radius(double xi, double g)
{
    return -1.0 / (xi * g);
}

/// @brief Reflected gracing angle
/// Eq. 5
/// @param theta Gracing angle
/// @return angle out
constexpr double reflect(double theta)
{
    const double a = 0.0;       // angle of layer [rad]
    return theta + 2 * a;
}

//
///////////////////////////////////////

// Trace ray through layer with linear sound speed
struct linear_layer_tracer
{
    const double _z0, _z1;       // next_depth of layer boundaries
    const double _c0, _c1;       // sound speed at layer boundaries

    // Ray always enters at boundary 0 (implementation doesn't care what's up or down)
    linear_layer_tracer(layer_boundary b0, layer_boundary b1) :
        _z0(b0.z), _z1(b1.z), _c0(b0.c), _c1(b1.c)
    {
        assert(std::fabs(_z0-_z1) > 1.0e-8);
    }

    /// @brief Trace ray through layer
    /// Ray always enters at boundary 0
    /// @param cos_th0 Cosine of gracing angle to boundary layer
    /// @return 
    trace_t trace(double cos_theta) const
    {
        double g = (_c1 - _c0) / (_z1 - _z0);    // Sound speed gradient
        double cos_th0 = cos_theta;
        double sin_th0 = std::sqrt(1.0 - cos_th0 * cos_th0);
        double xi = cos_th0 / _c0;               // Ray parameter Eq. (1)
        double cos_th1 = xi * _c1;               // From Eq. 1
        double sin_th1 = std::sqrt(1.0 - cos_th1 * cos_th1);

        // constant sound speed
        if (std::fabs(g) < 1.0e-8)
        {
            // ray as a straight line
            // Eq. 11 & 12 in https://ffi-publikasjoner.archive.knowledgearc.net/bitstream/handle/20.500.12242/2128/08-00610.pdf?sequence=1&isAllowed=y
            double dz = std::fabs(_z1 - _z0);
            double dt = dz / (_c1 * sin_th1);
            double dr = dz * cos_th1 / sin_th1;
            bool turns = false;
            assert(dr > 0);
            assert(dt > 1.0e-8);
            return { dt, dr, cos_theta, turns };
        }

        // Ray going up
        if (_z0 > _z1)
        {
            sin_th1 = -sin_th1;
            sin_th0 = -sin_th0;
        }

        const bool turns = cos_th1 >= 1.0;
        if (!turns)
        {
            // Travel time equation only works for positive gradients (dt_log_arg > 1).
            // If gradient is negative calculate travel time in opposite direction.
            double dt_log_arg = (_c1 / _c0) * (1.0 + sin_th0) / (1.0 + sin_th1);
            if (g < 0)
                dt_log_arg = 1 / dt_log_arg;      // swap direction
            assert(dt_log_arg >= 1);
            double dt = 1.0 / std::fabs(g) * std::log(dt_log_arg);
            assert(dt > 1.0e-8);
            double dr = 1.0 / (xi * g) * (sin_th0 - sin_th1);
            assert(dr > 0);
            return { dt, dr, cos_th1, turns };
        }
        else // ray turns
        {
            double dt = 2.0 / std::fabs(g) * std::log((1.0 + std::fabs(sin_th0)) / (cos_th0));
            assert(dt > 1.0e-8);
            double dr = 2.0 / (xi * g) * sin_th0;
            assert(dr > 0);
            return { dt, dr, cos_theta, turns };
        }
    }


#if 0
    /// @brief Range increments (Horizontal delta)
    /// Eq. (17) and (19)
    /// @param xi ray parameter (Greek letter Xi)
    /// @return horizontal increment
    double horizontal_delta(double xi, double sin_th0, double sin_th1) const
    {
        //auto xic0 = xi * xi * _c0 * _c0;
        //auto xic1 = xi * xi * _c1 * _c1;
        auto cos_th1 = std::sqrt(1.0 - sin_th1 * sin_th1);
        if (cos_th1 < 1.0)
        {
            //auto sin_th0 = std::sqrt(1.0 - xic0);
            //auto sin_th1 = std::sqrt(1.0 - xic1);
            auto r0 = ray_curvature_radius(xi);
            return r0 * (sin_th0 - sin_th1);
        }
        else // ray turns
        {
            //auto sin_th0 = std::sqrt(1.0 - xic0);
            auto r0 = ray_curvature_radius(xi);
            return 2.0 * r0 * sin_th0;
        }
    }
#endif

#if 0
    /// @brief Layer travel time
    /// Eq. (18) and (20)
    /// @param xi ray parameter
    /// @return duration
    double travel_time_delta(double xi, double sin_th0, double sin_th1) const
    {
        //const auto xic0 = xi * xi * _c0 * _c0;
        //const auto xic1 = xi * xi * _c1 * _c1;
        //auto cos_th0 = std::sqrt(1.0 - sin_th0*sin_th0);
        //auto cos_th1 = std::sqrt(1.0 - sin_th1*sin_th1);
        //double cos_th1 = xi * _c1;                  // Eq. 1
        if (cos_th1 < 1.0)
        {
            //auto sin_th0 = std::sqrt(1.0 - xic0);
            //auto sin_th1 = std::sqrt(1.0 - xic1);
            auto ln_expr = (cos_th1 / cos_th0) * (1.0 + sin_th0) / (1.0 + sin_th1);
            // assert(ln_expr >= 1.0);
            return 1.0 / /*std::fabs*/(g) * std::log(ln_expr);
        }
        else // ray turns
        {
            //auto sin_th0 = std::sqrt(1.0 - xic0);
            double ln_expr = (1.0 + sin_th0) / (cos_th0);
            return 2.0 * std::fabs(g) * std::log(ln_expr);
        }
    }
#endif
};

// Sound speed profile
// key: next_depth [m]
// value: sound speed [m/s]
using sound_speed_profile = std::map<double, double>;


inline double sound_speed(const sound_speed_profile& svp, double z)
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
inline double next_depth(double z, double dz)
{
    constexpr double e = 1.0e-9;
    if (dz > 0) // round up to multiple of dz
    {
        return std::ceil((z + e) / dz) * dz;
    }
    else        // round down to multiple of |dz|
    {
        dz = std::fabs(dz);
        return std::floor((z - e) / dz) * dz;
    }
}

class ray_tracer
{
public:

    explicit ray_tracer(const sound_speed_profile& svp) :
        _svp(svp)
    {
    }

    // Trace ray using sound speed profile to apply refraction
    // @param z Start next_depth for trace [m]
    // @param theta gracing angle [rad]
    // @param duration Time limit for trace [sec]
    // @return Ray path
    ray_path_t trace(double z, double theta, double duration, double dz = 50.0)
    {
        if (_svp.empty())
            throw std::out_of_range("Sound speed profile is empty");

        ray_path_t path;
        path.reserve(10000);
        trace_pos pos{ z, 0 };
        path.push_back(pos);

        auto z0 = z;
        auto c0 = sound_speed(_svp, z0);
        auto cos_th0 = std::cos(theta);
        while (duration > 0)
        {
            auto z1 = next_depth(z0,  dz);      // = z0 + dz;
            auto c1 = sound_speed(_svp, z1);

            linear_layer_tracer layer_tracer(
                layer_boundary{ z0, c0 },
                layer_boundary{ z1, c1 });

            auto [dt, dr, cos_th1, turns] = layer_tracer.trace(cos_th0);
            if (dt <= 0) break;
            assert(dr > 0.0);
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
        //path.pop_back();  // TODO invalid r
        return path;
    }

private:

    const sound_speed_profile& _svp;
};


} // end namespace urt
