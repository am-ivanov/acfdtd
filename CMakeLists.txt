# Define an executable
add_executable(acfdtd
    src/acfdtd.cpp
    src/config.cpp)

# Define the libraries this project depends upon
target_link_libraries(acfdtd PRIVATE RGrid)
target_link_libraries(acfdtd PRIVATE kutils)
