SET( FTKCOMMON_SRCS
  ftkUtils.cpp
)

SET( FTKCOMMON_HDRS
  ftkUtils.h
)

# GnT -> Graphs and Tree Classes
# Contains RAG and Merge Tree Codes
ADD_LIBRARY( ftkCommon
            ${FTKCOMMON_SRCS} 
            ${FTKCOMMON_HDRS} )
            
TARGET_LINK_LIBRARIES( ftkCommon TinyXML ${VTK_LIBRARIES} ${ITK_LIBRARIES} )