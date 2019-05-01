########### next target ###############

file(GLOB header *.h)
file(GLOB sources *.cpp)

add_executable(HybridizeTest main.cpp ${sources} ${header})

target_include_directories(HybridizeTest PUBLIC)

target_link_libraries(HybridizeTest pz ${OpenCV_LIBS})


