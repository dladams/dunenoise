# dunenoise/Fcl/Protodune/CMakeLists.txt

# Install fcl files.
install(
  DIRECTORY .
  DESTINATION fcl
  FILES_MATCHING PATTERN "*.fcl"
)
