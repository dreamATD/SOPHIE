set(BINARY ${CMAKE_PROJECT_NAME}_tst)


add_executable(${BINARY} TestUtils.cpp Utils.cpp Context.cpp Poly.cpp Encoder.cpp Evaluator.cpp Encryptor.cpp)

target_link_libraries(${BINARY} ${CMAKE_PROJECT_NAME}_lib)