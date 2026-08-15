cmake_minimum_required(VERSION 3.10)
project(decom_h)

# 设置 C++ 标准
set(CMAKE_CXX_STANDARD 14)

# 查找 pybind11 和 igraph 库
find_package(pybind11 REQUIRED)
find_package(igraph REQUIRED)

# pybind11 封装 C++ + 算法实现 C
set(SRC_FILES decom_h.cpp decom_h.c)

# 平台特定编译选项
if (MSVC)  # Windows 编译器
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /O2 /Ot /GL")
    set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} /O2 /Ot /GL")
    set(CMAKE_EXE_LINKER_FLAGS_RELEASE "${CMAKE_EXE_LINKER_FLAGS_RELEASE} /LTCG")
elseif(UNIX)  # Linux 或 macOS
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -Ofast -flto")
    set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -O3 -Ofast -flto")
    set(CMAKE_EXE_LINKER_FLAGS_RELEASE "${CMAKE_EXE_LINKER_FLAGS_RELEASE} -flto")
endif()

# 创建 Python 扩展模块
add_library(decom_h MODULE ${SRC_FILES})

# 0416 连通分量版本扩展模块（用于性能对比）
# add_library(decom_h_0416 MODULE decom_h_0416.cpp "decom_h_0416连通分量.c")

# 链接 pybind11 和 igraph 库
target_link_libraries(decom_h PRIVATE
    pybind11::module
    igraph::igraph
    igraph::bliss
    igraph::prpack
)

# target_link_libraries(decom_h_0416 PRIVATE
#     pybind11::module
#     igraph::igraph
#     igraph::bliss
#     igraph::prpack
# )

# 指定 include 目录
target_include_directories(decom_h PRIVATE
    ${PYBIND11_INCLUDE_DIRS}
    ${IGRAPH_INCLUDE_DIRS}
    .
    ../include
)

# target_include_directories(decom_h_0416 PRIVATE
#     ${PYBIND11_INCLUDE_DIRS}
#     ${IGRAPH_INCLUDE_DIRS}
#     .
#     ../include
# )

# 根据操作系统设置不同的输出文件后缀
if (WIN32)  # Windows 系统
    set_target_properties(decom_h PROPERTIES
        PREFIX ""
        SUFFIX ".pyd"
        LIBRARY_OUTPUT_DIRECTORY_RELEASE "${CMAKE_BINARY_DIR}"
        LIBRARY_OUTPUT_DIRECTORY_DEBUG "${CMAKE_BINARY_DIR}"
        LIBRARY_OUTPUT_DIRECTORY_RELWITHDEBINFO "${CMAKE_BINARY_DIR}"
        LIBRARY_OUTPUT_DIRECTORY_MINSIZEREL "${CMAKE_BINARY_DIR}"
    )
    # set_target_properties(decom_h_0416 PROPERTIES
    #     PREFIX ""
    #     SUFFIX ".pyd"
    #     LIBRARY_OUTPUT_DIRECTORY_RELEASE "${CMAKE_BINARY_DIR}"
    #     LIBRARY_OUTPUT_DIRECTORY_DEBUG "${CMAKE_BINARY_DIR}"
    #     LIBRARY_OUTPUT_DIRECTORY_RELWITHDEBINFO "${CMAKE_BINARY_DIR}"
    #     LIBRARY_OUTPUT_DIRECTORY_MINSIZEREL "${CMAKE_BINARY_DIR}"
    # )
elseif(UNIX)  # Linux 系统
    set_target_properties(decom_h PROPERTIES
        PREFIX "lib"
        SUFFIX ".so"
        LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}"
    )
    # set_target_properties(decom_h_0416 PROPERTIES
    #     PREFIX "lib"
    #     SUFFIX ".so"
    #     LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}"
    # )
endif()

