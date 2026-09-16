#include "ERF_PBLDiagnosticRecorder.H"

#include <gtest/gtest.h>

TEST(PBLDiagnosticRecorder, SiteFilenamesAndDefaultIntervals)
{
    EXPECT_EQ(PBLDiagnosticRecorder::make_filename("pbl_diag_", 0, false),
              "pbl_diag_00_summary.txt");
    EXPECT_EQ(PBLDiagnosticRecorder::make_filename("pbl_diag_", 12, true),
              "pbl_diag_12_profile.txt");
    PBLDiagnosticRecorder recorder;
    EXPECT_TRUE(recorder.summary_step(0));
    EXPECT_TRUE(recorder.summary_step(1));
    EXPECT_TRUE(recorder.profile_step(0));
    EXPECT_FALSE(recorder.profile_step(1));
}

TEST(PBLDiagnosticRecorder, UsesInstrumentMovementConvention)
{
    InstrumentSimUtil::MoveSchedule move;
    move.start_time = {10.0, 30.0};
    move.stop_time = {20.0, 40.0};
    move.speed_x = {2.0, -1.0};
    move.speed_y = {0.5, 3.0};
    EXPECT_DOUBLE_EQ(move.offset_x(5.0), 0.0);
    EXPECT_DOUBLE_EQ(move.offset_x(25.0), 20.0);
    EXPECT_DOUBLE_EQ(move.offset_x(35.0), 15.0);
    EXPECT_DOUBLE_EQ(move.offset_y(35.0), 20.0);
}
