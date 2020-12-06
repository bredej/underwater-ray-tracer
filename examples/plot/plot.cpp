#include <underwater-ray-tracer/ray_tracer.hpp>

#define GNUPLOT_DEPRECATE_WARN
#include <gnuplot-iostream/gnuplot-iostream.h>


using namespace gnuplotio;


int main()
{
    using namespace urt;

    sound_speed_profile profile;
    profile.insert(0, 1400);
    profile.insert(20, 1500);
    profile.insert(60, 1300);
    profile.insert(200, 1200);
    profile.insert(10000, 1200);

    ray_tracer tracer(profile);

    Gnuplot gp;
    gp << "set xrange [-5:150]\nset yrange [300 : 0] reverse\n";
    gp << "plot '-' with lines title 'Ray path'\n";

    //for (auto ang = 10.0_deg; ang <= 85.1_deg; ang += 5.0_deg)
    auto ang = 20.0_deg;        // grazing angle
    {
        vector3d start_pos{ 0, 0, 10 };
        //vector3d start_dir = rotate_vector(z_axis, quaternion(-x_axis, ang));
        double duration = 150.0 / 1480.0;   // ca. 150m
        auto path = tracer.trace(start_pos, ang, duration);

        std::vector<std::pair<double, double>> curve;
        curve.reserve(path.size());
        for (auto& pos : path)
        {
            curve.emplace_back(pos.y(), pos.z());
        }

        gp.send1d(curve);
    }
}
