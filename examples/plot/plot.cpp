#define _USE_MATH_DEFINES
#include <cmath>
#include <underwater-ray-tracer/ray_tracer.hpp>

#define GNUPLOT_DEPRECATE_WARN
#include <gnuplot-iostream/gnuplot-iostream.h>


using namespace gnuplotio;

using curve_t = std::vector<std::pair<double, double>>;
using curves_t = std::vector<curve_t>;



int main()
{
    using namespace urt;
    using float_t = float;

    sound_speed_profile<float_t> profile;
    profile.insert({ float_t(0.0),    float_t(1500.0) });
    profile.insert({ float_t(750.0),  float_t(1462.5) });
    profile.insert({ float_t(4000.0), float_t(1517.75) });

    ray_tracer tracer(profile);

    Gnuplot gp;
    gp << "set xrange [-5:60000]\nset yrange [4500 : 0] reverse\n";
    gp << "plot '-' with lines title 'Ray path'\n";

    float_t z = 100;
    float_t duration = float_t(60000.0 / 1480.0);

    curves_t curves;
    curves.reserve(100);

    for (float_t ang = 1.0; ang<12.1; ang+=1.0)
    {
        float_t theta = ang * float_t(M_PI / 180.0);        // grazing angle
        auto path = tracer.trace(z, theta, duration, 5.0);

        std::vector<std::pair<double, double>> curve;
        curve.reserve(path.size());
        for (auto& pos : path)
        {
            curve.emplace_back(pos.r, pos.z);
        }

        curves.emplace_back(std::move(curve));
    }
    gp.send2d(curves);
}
