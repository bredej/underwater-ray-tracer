#define _USE_MATH_DEFINES
#include <cmath>
#include <underwater-ray-tracer/ray_tracer.hpp>

#define GNUPLOT_DEPRECATE_WARN
#include <gnuplot-iostream/gnuplot-iostream.h>


using namespace gnuplotio;


int main()
{
    using namespace urt;

    sound_speed_profile profile;
    profile.insert({ 0.0, 1500.0 });
    profile.insert({ 750.0, 1462.5 });
    profile.insert({ 4000.0, 1517.75 });

    ray_tracer tracer(profile);

    Gnuplot gp;
    gp << "set xrange [-5:150]\nset yrange [300 : 0] reverse\n";
    gp << "plot '-' with lines title 'Ray path'\n";

    double theta = 0.00000001 * M_PI / 180.0;        // grazing angle
    double z = 100.0;
    double duration = 60000.0 / 1480.0;

    {
        auto path = tracer.trace(z, theta, duration);

        std::vector<std::pair<double, double>> curve;
        curve.reserve(path.size());
        for (auto& pos : path)
        {
            curve.emplace_back(pos.r, pos.z);
        }

        gp.send1d(curve);
    }
}
