#define _USE_MATH_DEFINES
#include <cmath>
#include "gplot++.h"
#include <underwater-ray-tracer/ray_tracer.hpp>

int main()
{
    using namespace urt;
    using float_t = float;

    sound_speed_profile<float_t> profile;
    profile.insert({ float_t(0.0),    float_t(1500.0) });
    profile.insert({ float_t(750.0),  float_t(1462.5) });
    profile.insert({ float_t(4000.0), float_t(1517.75) });

    ray_tracer tracer(profile);

    using line_t = std::pair<std::vector<double>, std::vector<double>>;
    std::vector<line_t> lines;
    for (float_t ang = 1.0; ang<12.1; ang+=1.0)
    {
        const float_t z = 100;
        const float_t duration = float_t(60000.0 / 1480.0);
        const float_t theta = ang * float_t(M_PI / 180.0);        // grazing angle
        auto path = tracer.trace(z, theta, duration, 5.0);

        // Make plot dataset
        std::vector<double> x, y;
        x.reserve(path.size());
        y.reserve(path.size());
        for (auto& pos : path)
        {
            x.push_back((double)pos.r);
            y.push_back((double)pos.z);
        }

        lines.push_back({std::move(x), std::move(y)});
    }

    Gnuplot plt{};
	plt.set_yrange(6000, -100);
    for (const auto& line : lines)
        plt.plot(line.first, line.second);
    plt.show();
}
