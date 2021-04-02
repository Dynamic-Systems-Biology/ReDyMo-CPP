#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include <map>
#include <memory>

#include "../include/util.hpp"

TEST(UtilTest, TestValidOpenclErrors)
{
    std::map<int, char *> opencl_errors;
    opencl_errors.insert(std::pair<int, char *>(0,"CL_SUCCESS"));
    opencl_errors.insert(std::pair<int, char *>(-1,"CL_DEVICE_NOT_FOUND"));
    opencl_errors.insert(std::pair<int, char *>(-2,"CL_DEVICE_NOT_AVAILABLE"));
    opencl_errors.insert(std::pair<int, char *>(-3,"CL_COMPILER_NOT_AVAILABLE"));
    opencl_errors.insert(std::pair<int, char *>(-4,"CL_MEM_OBJECT_ALLOCATION_FAILURE"));
    opencl_errors.insert(std::pair<int, char *>(-5,"CL_OUT_OF_RESOURCES"));
    opencl_errors.insert(std::pair<int, char *>(-6,"CL_OUT_OF_HOST_MEMORY"));
    opencl_errors.insert(std::pair<int, char *>(-7,"CL_PROFILING_INFO_NOT_AVAILABLE"));
    opencl_errors.insert(std::pair<int, char *>(-8,"CL_MEM_COPY_OVERLAP"));
    opencl_errors.insert(std::pair<int, char *>(-9,"CL_IMAGE_FORMAT_MISMATCH"));
    opencl_errors.insert(std::pair<int, char *>(-10,"CL_IMAGE_FORMAT_NOT_SUPPORTED"));
    opencl_errors.insert(std::pair<int, char *>(-11,"CL_BUILD_PROGRAM_FAILURE"));
    opencl_errors.insert(std::pair<int, char *>(-12,"CL_MAP_FAILURE"));
    opencl_errors.insert(std::pair<int, char *>(-13,"CL_MISALIGNED_SUB_BUFFER_OFFSET"));
    opencl_errors.insert(std::pair<int, char *>(-14,"CL_EXEC_STATUS_ERROR_FOR_EVENTS_IN_WAIT_LIST"));
    opencl_errors.insert(std::pair<int, char *>(-15,"CL_COMPILE_PROGRAM_FAILURE"));
    opencl_errors.insert(std::pair<int, char *>(-16,"CL_LINKER_NOT_AVAILABLE"));
    opencl_errors.insert(std::pair<int, char *>(-17,"CL_LINK_PROGRAM_FAILURE"));
    opencl_errors.insert(std::pair<int, char *>(-18,"CL_DEVICE_PARTITION_FAILED"));
    opencl_errors.insert(std::pair<int, char *>(-19,"CL_KERNEL_ARG_INFO_NOT_AVAILABLE"));

    opencl_errors.insert(std::pair<int, char *>(-30,"CL_INVALID_VALUE"));
    opencl_errors.insert(std::pair<int, char *>(-31,"CL_INVALID_DEVICE_TYPE"));
    opencl_errors.insert(std::pair<int, char *>(-32,"CL_INVALID_PLATFORM"));
    opencl_errors.insert(std::pair<int, char *>(-33,"CL_INVALID_DEVICE"));
    opencl_errors.insert(std::pair<int, char *>(-34,"CL_INVALID_CONTEXT"));
    opencl_errors.insert(std::pair<int, char *>(-35,"CL_INVALID_QUEUE_PROPERTIES"));
    opencl_errors.insert(std::pair<int, char *>(-36,"CL_INVALID_COMMAND_QUEUE"));
    opencl_errors.insert(std::pair<int, char *>(-37,"CL_INVALID_HOST_PTR"));
    opencl_errors.insert(std::pair<int, char *>(-38,"CL_INVALID_MEM_OBJECT"));
    opencl_errors.insert(std::pair<int, char *>(-39,"CL_INVALID_IMAGE_FORMAT_DESCRIPTOR"));
    opencl_errors.insert(std::pair<int, char *>(-40,"CL_INVALID_IMAGE_SIZE"));
    opencl_errors.insert(std::pair<int, char *>(-41,"CL_INVALID_SAMPLER"));
    opencl_errors.insert(std::pair<int, char *>(-42,"CL_INVALID_BINARY"));
    opencl_errors.insert(std::pair<int, char *>(-43,"CL_INVALID_BUILD_OPTIONS"));
    opencl_errors.insert(std::pair<int, char *>(-44,"CL_INVALID_PROGRAM"));
    opencl_errors.insert(std::pair<int, char *>(-45,"CL_INVALID_PROGRAM_EXECUTABLE"));
    opencl_errors.insert(std::pair<int, char *>(-46,"CL_INVALID_KERNEL_NAME"));
    opencl_errors.insert(std::pair<int, char *>(-47,"CL_INVALID_KERNEL_DEFINITION"));
    opencl_errors.insert(std::pair<int, char *>(-48,"CL_INVALID_KERNEL"));
    opencl_errors.insert(std::pair<int, char *>(-49,"CL_INVALID_ARG_INDEX"));
    opencl_errors.insert(std::pair<int, char *>(-50,"CL_INVALID_ARG_VALUE"));
    opencl_errors.insert(std::pair<int, char *>(-51,"CL_INVALID_ARG_SIZE"));
    opencl_errors.insert(std::pair<int, char *>(-52,"CL_INVALID_KERNEL_ARGS"));
    opencl_errors.insert(std::pair<int, char *>(-53,"CL_INVALID_WORK_DIMENSION"));
    opencl_errors.insert(std::pair<int, char *>(-54,"CL_INVALID_WORK_GROUP_SIZE"));
    opencl_errors.insert(std::pair<int, char *>(-55,"CL_INVALID_WORK_ITEM_SIZE"));
    opencl_errors.insert(std::pair<int, char *>(-56,"CL_INVALID_GLOBAL_OFFSET"));
    opencl_errors.insert(std::pair<int, char *>(-57,"CL_INVALID_EVENT_WAIT_LIST"));
    opencl_errors.insert(std::pair<int, char *>(-58,"CL_INVALID_EVENT"));
    opencl_errors.insert(std::pair<int, char *>(-59,"CL_INVALID_OPERATION"));
    opencl_errors.insert(std::pair<int, char *>(-60,"CL_INVALID_GL_OBJECT"));
    opencl_errors.insert(std::pair<int, char *>(-61,"CL_INVALID_BUFFER_SIZE"));
    opencl_errors.insert(std::pair<int, char *>(-62,"CL_INVALID_MIP_LEVEL"));
    opencl_errors.insert(std::pair<int, char *>(-63,"CL_INVALID_GLOBAL_WORK_SIZE"));
    opencl_errors.insert(std::pair<int, char *>(-64,"CL_INVALID_PROPERTY"));
    opencl_errors.insert(std::pair<int, char *>(-65,"CL_INVALID_IMAGE_DESCRIPTOR"));
    opencl_errors.insert(std::pair<int, char *>(-66,"CL_INVALID_COMPILER_OPTIONS"));
    opencl_errors.insert(std::pair<int, char *>(-67,"CL_INVALID_LINKER_OPTIONS"));
    opencl_errors.insert(std::pair<int, char *>(-68,"CL_INVALID_DEVICE_PARTITION_COUNT"));

    opencl_errors.insert(std::pair<int, char *>(-1000,"CL_INVALID_GL_SHAREGROUP_REFERENCE_KHR"));
    opencl_errors.insert(std::pair<int, char *>(-1001,"CL_PLATFORM_NOT_FOUND_KHR"));
    opencl_errors.insert(std::pair<int, char *>(-1002,"CL_INVALID_D3D10_DEVICE_KHR"));
    opencl_errors.insert(std::pair<int, char *>(-1003,"CL_INVALID_D3D10_RESOURCE_KHR"));
    opencl_errors.insert(std::pair<int, char *>(-1004,"CL_D3D10_RESOURCE_ALREADY_ACQUIRED_KHR"));
    opencl_errors.insert(std::pair<int, char *>(-1005,"CL_D3D10_RESOURCE_NOT_ACQUIRED_KHR"));

    for (const auto &opencl_error : opencl_errors){
        ASSERT_EQ (getCLErrorString(opencl_error.first), opencl_error.second);
    }
}

TEST(UtilTest, TestInvalidOpenclErrors)
{
        ASSERT_EQ (getCLErrorString(1), "Unknown OpenCL error");
}

int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
