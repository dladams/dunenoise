# dunenoise/Fcl/CMakeLists.txt

# Install fcl files.
install(
  DIRECTORY .
  DESTINATION fcl
  FILES_MATCHING PATTERN "*.fcl"
)

add_subdirectory("Iceberg")
add_subdirectory("Protodune")
