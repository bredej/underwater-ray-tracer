#include <gtest/gtest.h>
#include <cmath>
#include <limits>
#include <underwater-ray-tracer/ray_tracer.hpp>

using namespace urt;

TEST(urt, sound_speed)
{
    sound_speed_profile<double> svp
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

TEST(urt, next_depth)
{
    EXPECT_FLOAT_EQ(300.0f, next_depth(0.0f, 300.0f));
    EXPECT_FLOAT_EQ(300.0f, next_depth(100.0f, 300.0f));
    EXPECT_FLOAT_EQ(300.0f, next_depth(200.0f, 300.0f));
    EXPECT_FLOAT_EQ(600.0f, next_depth(300.0f, 300.0f));
    EXPECT_FLOAT_EQ(600.0f, next_depth(400.0f, 300.0f));

    EXPECT_FLOAT_EQ(-300.0f, next_depth(0.0f, -300.0f));      // Should we throw if z < 0?
    EXPECT_FLOAT_EQ(0.0f, next_depth(100.0f, -300.0f));
    EXPECT_FLOAT_EQ(0.0f, next_depth(200.0f, -300.0f));
    EXPECT_FLOAT_EQ(0.0f, next_depth(300.0f, -300.0f));
    EXPECT_FLOAT_EQ(300.0f, next_depth(400.0f, -300.0f));

    EXPECT_DOUBLE_EQ(300.0, next_depth(0.0, 300.0));
    EXPECT_DOUBLE_EQ(300.0, next_depth(100.0, 300.0));
    EXPECT_DOUBLE_EQ(300.0, next_depth(200.0, 300.0));
    EXPECT_DOUBLE_EQ(600.0, next_depth(300.0, 300.0));
    EXPECT_DOUBLE_EQ(600.0, next_depth(400.0, 300.0));

    EXPECT_DOUBLE_EQ(-300.0, next_depth(0.0, -300.0));      // Should we throw if z < 0?
    EXPECT_DOUBLE_EQ(0.0, next_depth(100.0, -300.0));
    EXPECT_DOUBLE_EQ(0.0, next_depth(200.0, -300.0));
    EXPECT_DOUBLE_EQ(0.0, next_depth(300.0, -300.0));
    EXPECT_DOUBLE_EQ(300.0, next_depth(400.0, -300.0));
}
