macro(fetchduneanaobj DUNE_ANAOBJ_BRANCH)

  include(FetchContent)
  FetchContent_Declare(
    duneanaobj
    GIT_REPOSITORY https://github.com/DUNE/duneanaobj.git
    GIT_TAG        ${DUNE_ANAOBJ_BRANCH}
  )


  # Check if population has already been performed
  FetchContent_GetProperties(duneanaobj)
  if(NOT duneanaobj_POPULATED)
    cmake_policy(SET CMP0169 OLD)
    # Fetch the content using previously declared details
    FetchContent_Populate(duneanaobj)

  endif()

  include_directories(${duneanaobj_SOURCE_DIR})

  ROOT_GENERATE_DICTIONARY(StandardRecordDict
    ${duneanaobj_SOURCE_DIR}/duneanaobj/StandardRecord/classes.h
    LINKDEF
      ${duneanaobj_SOURCE_DIR}/duneanaobj/StandardRecord/classes_def.xml)

  add_library(libStandardRecordDict SHARED ${CMAKE_CURRENT_BINARY_DIR}/StandardRecordDict.cxx)
  target_link_libraries(libStandardRecordDict PUBLIC ROOT::RIO)
  set_target_properties(libStandardRecordDict PROPERTIES OUTPUT_NAME "StandardRecordDict")

  target_include_directories(libStandardRecordDict PUBLIC
    $<BUILD_INTERFACE:${duneanaobj_SOURCE_DIR}>
    $<INSTALL_INTERFACE:include>)
  target_include_directories(libStandardRecordDict PRIVATE
    $<BUILD_INTERFACE:${duneanaobj_SOURCE_DIR}/duneanaobj/StandardRecord> #root puts this in the dictionary
    )

  file(GLOB SR_IMPL_FILES "${duneanaobj_SOURCE_DIR}/duneanaobj/StandardRecord/*.cxx")
  add_library(duneanaobj_StandardRecord SHARED ${SR_IMPL_FILES})
  target_link_libraries(duneanaobj_StandardRecord PUBLIC libStandardRecordDict ROOT::MathCore ROOT::Physics)

  target_include_directories(duneanaobj_StandardRecord PUBLIC
    $<BUILD_INTERFACE:${duneanaobj_SOURCE_DIR}>
    $<INSTALL_INTERFACE:include>)

  set_target_properties(duneanaobj_StandardRecord PROPERTIES EXPORT_NAME StandardRecord)

  file(GLOB SR_HEADER_FILES "${duneanaobj_SOURCE_DIR}/duneanaobj/StandardRecord/*.h")

  install(FILES ${SR_HEADER_FILES} DESTINATION include/duneanaobj/StandardRecord)
  install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/libStandardRecordDict_rdict.pcm
    ${CMAKE_CURRENT_BINARY_DIR}/libStandardRecordDict.rootmap
    DESTINATION lib)

  install(TARGETS duneanaobj_StandardRecord libStandardRecordDict
    EXPORT duneanaobj-targets
    LIBRARY DESTINATION lib)

  install(EXPORT duneanaobj-targets
    FILE duneanaobjTargets.cmake
    NAMESPACE duneanaobj::
    DESTINATION ${CMAKE_INSTALL_PREFIX}/
  )

  add_library(duneanaobj::StandardRecord ALIAS duneanaobj_StandardRecord)
  add_library(duneanaobj::all ALIAS duneanaobj_StandardRecord)

endmacro()
