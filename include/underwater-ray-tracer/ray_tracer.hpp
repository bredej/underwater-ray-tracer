#pragma once

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

struct trace_t
{
    double r;           // Horizontal offset [m]
    double theta;       // grazing angle [rad]
    double dt;          // layer travel time [s]
};


struct layer_boundary
{
    double z;           // depth [m]
    double c;           // sound speed [m/s]
};

struct trace_pos
{
    double z;           // depth [m]
    double r;           // Horizontal offset [m]
};

using ray_path_t = std::vector<trace_pos>;


/// @brief Sound speed c(z)
/// Eq. (16)
/// Within the layer the sound speed profile is approximated as linear
/// @param z depth
/// @return sound speed

/// @brief Sound speed c(z)
/// @param z0 depth [m]
/// @param c0 sound speed [m/s]
/// @param g gradient [1/s]
/// @param z depth [m]
/// @return 
constexpr double sound_speed(double z0, double c0, double g, double z)
{
    return c0 + g * (z - z0);
}

/// @brief Sound speed gradient g
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


// Trace ray through uniform layer
struct uniform_layer_tracer
{
    const double _z0, _z1;       // depth of layer boundaries
    const double _c0, _c1;       // sound speed at layer boundaries
    const double _g;             // sound speed gradient


    uniform_layer_tracer(layer_boundary b0, layer_boundary b1) :
        _z0(b0.z), _z1(b1.z), _c0(b0.c), _c1(b1.c),
        _g(sound_speed_gradient())
    {
        assert(std::fabs(_z0-_z1) > 1.0e-8);
        assert(std::fabs(_g) > 1.0e-8);
    }

    /// @brief Trace ray through layer
    /// @param theta Gracing angle to boundary layer
    /// @return 
    trace_t trace(double theta) const
    {
        double xi = ray_parameter(theta, _c0);
        double r = horizontal_delta(xi);
        double dt = travel_time_delta(xi);
        double theta_out = std::copysign(std::acos(xi * _c1), _g);
        return { r, theta_out, dt };
    }

    /// @brief Ray parameter
    /// (Greek letter Xi)
    /// Eq. (1)
    /// @param theta gracing angle of ray at given depth (Greek theta)
    /// @param c sound speed
    /// @return ray parameter
    double ray_parameter(double theta, double c) const
    {
        return std::cos(theta) / c;
    }

    /// @brief Sound speed gradient g(z)
    /// Eq. (4)
    /// The positive or negative sign of the gradient determines whether
    /// the sign of R is negative or positive, and thereby determines
    /// if the ray path curves downward or upward.
    /// @return gradient
    constexpr double sound_speed_gradient() const
    {
        return (_c1 - _c0) / (_z1 - _z0);
    }

    /// @brief Reflected gracing angle
    /// Eq. 5
    /// @param theta Gracing angle
    /// @return angle out
    constexpr double reflect(double theta) const
    {
        const double a = 0.0;       // angle of layer [rad]
        return theta + 2 * a;
    }

    /// @brief Sound speed c(z)
    /// Eq. (16)
    /// Within the layer the sound speed profile is approximated as linear
    /// @param z depth
    /// @return sound speed
    constexpr double sound_speed(double z) const
    {
        return _c0 + _g * (z - _z0);
    }

    /// @brief Ray’s radius of curvature at given depth R(z)
    /// Eq. (3)
    /// @param xi ray parameter (Greek letter Xi)
    /// @return radius
    constexpr double ray_curvature_radius(double xi) const
    {
        return 1.0 / (xi * _g);
    }

    /// @brief Range increments (Horizontal delta)
    /// Eq. (17) and (19)
    /// @param xi ray parameter (Greek letter Xi)
    /// @return horizontal increment
    double horizontal_delta(double xi) const
    {
        auto xic0 = xi * xi * _c0 * _c0;
        auto xic1 = xi * xi * _c1 * _c1;
        auto r0 = ray_curvature_radius(xi);
        if (xic1 < 1.0)
        {
            auto sin_th0 = std::sqrt(1.0 - xic0);
            auto sin_th1 = std::sqrt(1.0 - xic1);
            return r0 * (sin_th0 - sin_th1);
        }
        else // ray turns
        {
            auto sin_th0 = std::sqrt(1.0 - xic0);
            return 2.0 * r0 * sin_th0;
        }
    }

    /// @brief Layer travel time
    /// Eq. (18) and (20)
    /// @param xi ray parameter
    /// @return duration
    double travel_time_delta(double xi) const
    {
        const auto xic0 = xi * xi * _c0 * _c0;
        const auto xic1 = xi * xi * _c1 * _c1;
        if (xic1 < 1.0)
        {
            auto sin_th0 = std::sqrt(1.0 - xic0);
            auto sin_th1 = std::sqrt(1.0 - xic1);
            //auto ln_expr = _c1 / _c0 * (1.0 + sin_th0) / (1.0 + sin_th1);
            // assert(ln_expr >= 1.0);
            return 1.0 / std::fabs(_g) * std::log(_c1/_c0  * (1.0 + sin_th0) / (1.0 + sin_th1));
        }
        else // ray turns
        {
            auto sin_th0 = std::sqrt(1.0 - xic0);
            double nom = 1 + sin_th0;
            double denom = xi * _c0;
            return 2.0 * std::fabs(_g) * std::log(nom / denom);
        }
    }
};

// Sound speed profile
// key: depth [m]
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

class ray_tracer
{
public:

    explicit ray_tracer(const sound_speed_profile& svp) :
        _svp(svp)
    {
    }

    // Trace ray using sound speed profile to apply refraction
    // @param z Start depth for trace [m]
    // @param theta gracing angle [rad]
    // @param duration Time limit for trace [sec]
    // @return Ray path
    ray_path_t trace(double z, double theta, double duration)
    {
        if (_svp.empty())
            throw std::out_of_range("Sound speed profile is empty");

        ray_path_t path;
        path.reserve(10000);
        trace_pos pos{ z, 0 };
        path.push_back(pos);

        double dz = 10.0;
        auto z0 = z;
        auto c0 = sound_speed(_svp, z0);
        while (duration > 0)
        {
            auto z1 = z0 + dz;
            auto c1 = sound_speed(_svp, z1);

            uniform_layer_tracer layer_tracer(
                layer_boundary{ z0, c0 },
                layer_boundary{ z1, c1 });

            auto [r, th, dt] = layer_tracer.trace(theta);
            pos.r += r;
            pos.z = z1;     // or z0
            path.push_back(pos);

            if (std::signbit(th) != std::signbit(theta))
                dz = -dz;

            theta = th;
            duration -= dt;
            z0 = z1;
            c0 = c1;
        }
        return path;
    }

private:

    const sound_speed_profile& _svp;
};


} // end namespace urt
