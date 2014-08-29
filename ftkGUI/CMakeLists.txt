SET( GUI_SRCS
    PlotWindow.cpp
    ScatterView.cpp
    TableWindow.cpp
    HistoWindow.cpp
    ObjectSelection.cpp
	StatisticsToolbar.cpp
	GraphWindowForNewSelection.cpp
	GraphWindow.cpp
	SelectiveClustering.cpp
	QvtkTableView.cpp
    vtk3DView.cpp
)
SET( NONQT_GUI_SRCS
	SelectionUtilities.cpp
)

SET( GUI_HDRS
    PlotWindow.h
    ScatterView.h
    TableWindow.h
    HistoWindow.h
    ObjectSelection.h
	StatisticsToolbar.h
	GraphWindowForNewSelection.h
	GraphWindow.h
	SelectiveClustering.h
	QvtkTableView.h
    vtk3DView.h
)
SET( NONQT_GUI_HDRS
	SelectionUtilities.h
)

SET( GUI_LIBS
  ${QT_LIBRARIES} ftkCommon vnl)

SET( GUI_UIS
)

SET( GUI_RCS
)

# this command will generate rules that will run rcc on all files from SAMPLE_RCS
# in result SAMPLE_RC_SRCS variable will contain paths to files produced by rcc
QT4_ADD_RESOURCES( GUI_RC_SRCS ${GUI_RCS} )

# this will run uic on .ui files:
QT4_WRAP_UI( GUI_UI_HDRS ${GUI_UIS} )

# and finally this will run moc:
QT4_WRAP_CPP( GUI_MOC_SRCS ${GUI_HDRS} )

# This must be built static because libKPLS is only available as a static
# library.
ADD_LIBRARY( ftkGUI STATIC
  ${GUI_SRCS}
  ${GUI_MOC_SRCS}
  ${GUI_HDRS}
  ${GUI_UI_HDRS}
  ${GUI_RC_SRCS}
  ${NONQT_GUI_SRCS}
  ${NONQT_GUI_HDRS}
  )
# For linking into shared libraries.
IF( CMAKE_COMPILER_IS_GNUCC )
  SET_TARGET_PROPERTIES( ftkGUI PROPERTIES COMPILE_FLAGS -fPIC )
ENDIF()

TARGET_LINK_LIBRARIES( ftkGUI ${GUI_LIBS} )
