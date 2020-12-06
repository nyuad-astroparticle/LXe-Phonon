# Define directories
SRC_DIR := src
OBJ_DIR := obj
BIN_DIR := bin

# Define executable
EXE := $(BIN_DIR)/particle_sound

# get the source files
SRC := $(wildcard $(SRC_DIR)/*.cpp)

# get the object fiels from source files
OBJ := $(SRC:$(SRC_DIR)/%.cpp=$(OBJ_DIR)/%.o)

# For flags
CPPFLAGS := -Iinclude -MMD -MP # -I is a preprocessor flag, not a compiler flag
CFLAGS   := -Wall -g             # some warnings about bad code
LDFLAGS  := -Llib              # -L is a linker flag
LDLIBS   := -lm                # Left empty if no libs are needed
CC       := g++

# To create the targets
.PHONY: all	clean # We tell cmake not to create a file called all

all: $(EXE)

# For the executable we need the object which we get bu compliling them
$(EXE): $(OBJ) | $(BIN_DIR)
	$(CC) $(LDFLAGS) $^ $(LDLIBS) -o $@

# To create the object files from the c files
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp | $(OBJ_DIR)
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

# Create the relevant directories
$(BIN_DIR) $(OBJ_DIR):
	mkdir -p $@

# To clean everything
clean:
	@$(RM) -rv $(BIN_DIR) $(OBJ_DIR)

-include $(OBJ:.o=.d)