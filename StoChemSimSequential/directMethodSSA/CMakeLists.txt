# sequential direct method CMakeLists
add_library(directMethod STATIC 
			dependencyGraph.cpp
			dependencyGraph.h
			directMethodSSA.cpp
			directMethodSSA.h
			reactionTree.cpp
			reactionTree.h)