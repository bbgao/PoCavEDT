#=================================================================
#
#	Author: Sebastian Daberdaku
#	Project: PoCavEDT	v1.0
#
#	This makefile searches SRC_DIR and its subdirectories
#	recursively for source files and builds them.
#
#	Object and dependency files can be placed in separate 
#	directories.
#
#=================================================================
CC := g++ -m64
MKDIR_P := @mkdir -p 

GPPOPT := -c -std=c++11 -MMD -O3 -funroll-loops 
LDOPT := -static -std=c++11 -O3 

BOOST_LIBS := -lboost_program_options -lboost_regex -lpthread -lboost_system -lboost_filesystem
LIBS := $(BOOST_LIBS)

BIN_DIR := bin
OBJ_DIR := bin/obj
DEP_DIR := bin/dep
OUT_DIR := bin/output
SRC_DIR := src

# Defines: just add the defines to this variable
DEFS := -D _GLIBCXX_PARALLEL -D NDEBUG #-D NO_OUTPUT_TEST #-D RANGECHECK_TEST #-D TIGHT_PACKING

# Make does not offer a recursive wildcard function, so here's one:
rwildcard = $(wildcard $1$2) $(foreach d,$(wildcard $1*),$(call rwildcard,$d/,$2))

# Recursively find all source files in SRC_DIR
SOURCES := $(call rwildcard,$(SRC_DIR)/,*.cpp)
OBJECTS := $(SOURCES:$(SRC_DIR)/%.cpp=$(OBJ_DIR)/%.o)
DEPENDENCIES := $(SOURCES:$(SRC_DIR)/%.cpp=$(DEP_DIR)/%.d)
EXECUTABLE := $(BIN_DIR)/PoCavEDT

.PHONY: all clean	

all: $(EXECUTABLE) 
	
$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDOPT) $^ -o $@ $(LIBS) 
	$(MKDIR_P) $(OUT_DIR)
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(MKDIR_P) $(@D) $(@D:$(OBJ_DIR)%=$(DEP_DIR)%)
	$(CC) $(GPPOPT) $< -o $@ -MF $(@:$(OBJ_DIR)%.o=$(DEP_DIR)%.d) $(DEFS) $(LIBS)
	
clean: 
	rm -rf $(EXECUTABLE) $(OBJ_DIR) $(DEP_DIR)

-include $(DEPENDENCIES)
	