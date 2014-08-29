SET( STrenDEditor_SRCS
	STrenDEditor.cpp
)

SET( STrenDEditor_HDRS
	STrenDEditor.h
)

# and finally this will run moc:
QT4_WRAP_CPP( STrenDEditor_MOC_SRCS ${STrenDEditor_HDRS} )


ADD_EXECUTABLE(STrenDEditor main.cpp
             ${STrenDEditor_SRCS}
             ${STrenDEditor_HDRS} 
             ${STrenDEditor_MOC_SRCS} 
)

TARGET_LINK_LIBRARIES( STrenDEditor ${QT_LIBRARIES} ftkCommon STrenD ClusClus_lib)

#MAKE A GUI EXECUTABLE, NOT A CONSOLE APPLICATION
#SET_TARGET_PROPERTIES( STrenDEditor PROPERTIES WIN32_EXECUTABLE 1)
#OR MAKE A CONSOLE APPLICATION
SET_TARGET_PROPERTIES( STrenDEditor PROPERTIES WIN32_EXECUTABLE 0)
