add_library(_stVenant SHARED bem.cpp fem.cpp)
target_link_libraries(_stVenant ${LAPLACE_LIB_2D})

install(TARGETS _stVenant LIBRARY DESTINATION lib)

add_library(stVenant SHARED stVenantPy.cpp)
set_target_properties(stVenant PROPERTIES PREFIX "")
target_link_libraries(stVenant _stVenant)
set_target_properties(stVenant PROPERTIES INSTALL_RPATH ${CMAKE_INSTALL_PREFIX}/lib)
install(TARGETS stVenant DESTINATION ${python_libs})