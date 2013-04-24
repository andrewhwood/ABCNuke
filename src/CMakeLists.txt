#set PACKAGE_ROOT = "/tools/package"
#set NUKE_INCLUDE_DIR = "${PACKAGE_ROOT}/nuke/6.3v8/include"
#set Boost_INCLUDE_DIRS = "${PACKAGE_ROOT}/boost/1.47.0/include"
#set ALEMBIC_INCLUDE_DIR = "${PACKAGE_ROOT}/alembic/1.0.5.beta9/include"
#set OPENEXR_INCLUDE_PATHS =
#"${PACKAGE_ROOT}/ilmbase/1.0.2PIC/include:${PACKAGE_ROOT}/openexr/1.7.0/include"
#set ALEMBIC_ILMBASE_INCLUDE_DIRECTORY = "${PACKAGE_ROOT}/ilmbase/1.0.2PIC/include"
#set NUKE_LIBRARY_DIR = "${PACKAGE_ROOT}/nuke/6.3v8"
#set ALEMBIC_LIBRARY_DIR = "${PACKAGE_ROOT}/alembic/1.0.5.beta9/lib"


include_directories ( ${NUKE_INCLUDE_DIR}
					  ${Boost_INCLUDE_DIRS}
					  ${ALEMBIC_INCLUDE_DIR}	  
					  ${OPENEXR_INCLUDE_PATHS}
					  ${ALEMBIC_ILMBASE_INCLUDE_DIRECTORY}
					 )

link_directories    ( ${NUKE_LIBRARY_DIR}
				      ${ALEMBIC_LIBRARY_DIR}
				      #${OPENEXR_LIB_PATH}
				     )

add_library 		( ABCReadGeo SHARED
			  ABCNuke_ArchiveHelper.cpp
			  ABCNuke_Interpolation.cpp
			  ABCNuke_MatrixHelper.cpp
			  ABCNuke_GeoHelper.cpp
		          ABCReadGeo.cpp
				   	 )


set_target_properties ( ABCReadGeo
			PROPERTIES
			PREFIX ""
		        OUTPUT_NAME "ABCReadGeo"
			COMPILE_FLAGS "-g -c -Wall -DUSE_GLEW -fPIC -O3"
			  		   )

target_link_libraries ( ABCReadGeo 
		        ${DDIMAGE_LIBRARY}
		        GLEW
			Iex
			Half
			Imath
                        pthread
			${ALEMBIC_LIBRARIES}
			${HDF5_HL_LIBRARIES}
			${HDF5_LIBRARIES}
                      )


install ( TARGETS ABCReadGeo
          DESTINATION ${CMAKE_BINARY_DIR}/plugins/${NUKE_MAJOR_VERSION}.${NUKE_MINOR_VERSION})
