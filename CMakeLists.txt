cmake_minimum_required(VERSION 3.8)
project(dual_arm_tp_cooperative_control)

if(NOT CMAKE_CXX_STANDARD)
  set(CMAKE_CXX_STANDARD 17)
endif()


if(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  add_compile_options(-Wall -Wextra -Wpedantic)
endif()

# find dependencies
find_package(ament_cmake REQUIRED)
find_package(rclcpp REQUIRED)
find_package(controller_interface REQUIRED)
find_package(pluginlib REQUIRED)
find_package(rclcpp_lifecycle REQUIRED)
find_package(hardware_interface REQUIRED)
find_package(franka_semantic_components REQUIRED)
find_package(generate_parameter_library)
find_package(franka_msgs REQUIRED)
find_package(Eigen3 REQUIRED)

foreach(Dependency IN ITEMS ${THIS_PACKAGE_INCLUDE_DEPENDS})
  find_package(${Dependency} REQUIRED)
endforeach()
# uncomment the following section in order to fill in
# further dependencies manually.
# find_package(<dependency> REQUIRED)

add_library(dual_arm_tp_cooperative_control SHARED src/panda_arms.cpp)
target_include_directories(dual_arm_tp_cooperative_control 
			   PUBLIC 
			   include 
			   ${EIGEN3_INCLUDE_DIRS})
ament_target_dependencies(
			dual_arm_tp_cooperative_control
                        controller_interface
        		hardware_interface
        		pluginlib
        		rclcpp
        		rclcpp_lifecycle
        		franka_semantic_components
        		${THIS_PACKAGE_INCLUDE_DEPENDS})


install(TARGETS dual_arm_tp_cooperative_control DESTINATION lib)
install(DIRECTORY include/ DESTINATION include)

if(BUILD_TESTING)
  find_package(ament_lint_auto REQUIRED)
  set(CPP_DIRECTORIES src include)
  # the following line skips the linter which checks for copyrights
  # comment the line when a copyright and license is added to all source files
  set(ament_cmake_copyright_FOUND TRUE)
  # the following line skips cpplint (only works in a git repo)
  # comment the line when this package is in a git repo and when
  # a copyright and license is added to all source files
  set(ament_cmake_cpplint_FOUND TRUE)
  ament_lint_auto_find_test_dependencies()
endif()

ament_export_include_directories(
        include
)
ament_export_libraries(
        dual_arm_tp_cooperative_control
)
ament_export_dependencies(
        ${THIS_PACKAGE_INCLUDE_DEPENDS}
)

ament_package()
