file(REMOVE_RECURSE
  "../../lib/libOPENFHEcore.pdb"
  "../../lib/libOPENFHEcore.so"
  "../../lib/libOPENFHEcore.so.1"
  "../../lib/libOPENFHEcore.so.1.2.1"
)

# Per-language clean rules from dependency scanning.
foreach(lang C CXX)
  include(CMakeFiles/OPENFHEcore.dir/cmake_clean_${lang}.cmake OPTIONAL)
endforeach()
