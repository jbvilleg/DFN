########### next target ###############


add_executable(HybridizeTest main.cpp )



target_include_directories(HybridizeTest PUBLIC)

target_link_libraries(HybridizeTest pz ${OpenCV_LIBS})


