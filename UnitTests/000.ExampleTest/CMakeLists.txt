get_filename_component(ProjectId ${CMAKE_CURRENT_SOURCE_DIR} NAME)
#message(${ProjectId})
string(REPLACE " " "_" ProjectId ${ProjectId}) # In case there are blank spaces
project(${ProjectId})

include(../TestsCommon.cmake)
my_add_test(${ProjectId})
