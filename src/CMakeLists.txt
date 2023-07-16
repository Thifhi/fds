add_executable(ADS main.cpp)
target_include_directories(ADS PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
set_target_properties(${PROJECT_NAME} PROPERTIES OUTPUT_NAME "ads_program")