#include <gtest/gtest.h>
#include <underwater-ray-tracer/ray_tracer.hpp>

using namespace urt;

TEST(urt, sound_speed)
{
    sound_speed_profile svp
    {
        { 0.0, 1400.0 },
        { 100.0 , 1500.0 },
        { 1000 , 1600.0 }
    };

    EXPECT_DOUBLE_EQ(1400.0, sound_speed(svp, -1.0));
    EXPECT_DOUBLE_EQ(1400.0, sound_speed(svp, 0.0));
    EXPECT_DOUBLE_EQ(1450.0, sound_speed(svp, 50.0));
    EXPECT_DOUBLE_EQ(1500.0, sound_speed(svp, 100.0));
    EXPECT_DOUBLE_EQ(1550.0, sound_speed(svp, 550.0));
    EXPECT_DOUBLE_EQ(1600.0, sound_speed(svp, 1000.0));
    EXPECT_DOUBLE_EQ(1600.0, sound_speed(svp, 110000.0));
}

TEST(urt, sound_speed_4000m)
{
    double z0 = 750.0;
    double c0 = 1462.5;
    double z = 4000;
    double g = 0.017;
    double expected = 1517.75;
    double actual = c0 + g * (z - z0);
    EXPECT_DOUBLE_EQ(expected, actual);
}

TEST(uniform_layer_tracer, sound_speed_gradient)
{
    layer_boundary l0, l1;
    l0.z = 0.0;
    l0.c = 1500.0;
    l1.z = 750.0;
    l1.c = 1462.5;
    uniform_layer_tracer layer_tracer(l0, l1);
    EXPECT_DOUBLE_EQ(-0.05, layer_tracer.sound_speed_gradient());
}
