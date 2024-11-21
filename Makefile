.PHONY: all distclean

CFLAGS += -Wall -std=c2x -pedantic -Wno-newline-eof -I./include -march=native -Ofast
CFLAGS += -DVERBOSE

# Headers.
HEADERS = ./include/*.h

# Executables.
TESTS = $(subst src/,executables/,$(subst .c,.out,$(shell find src -name "Test_*.c")))

# Objects.
OBJECTS = $(subst src/,objects/,$(subst .c,.o,$(shell find src -name "NAQRA_*.c")))

# Directories.
DIRECTORIES = ./objects ./executables

# All.
all: $(DIRECTORIES) $(TESTS)
	@echo "Compiled everything!"

# Tests.
$(TESTS): executables/Test_%.out: objects/Test_%.o $(OBJECTS) 
	@echo "Linking to $@"
	@$(CC) $^ -o $@

# Objects.
$(OBJECTS): objects/%.o: src/%.c $(HEADERS)
	@echo "Compiling $<"
	@$(CC) $(CFLAGS) -c $< -o $@

$(subst src/,objects/,$(subst .c,.o,$(shell find src -name "Test_*.c"))): objects/%.o: src/%.c $(HEADERS)
	@echo "Compiling $<"
	@$(CC) $(CFLAGS) -c $< -o $@

# Directories.
$(DIRECTORIES):
	@mkdir -p $(DIRECTORIES)

# Clean.
distclean:
	@echo "Cleaning the repo."
	@$(RM) -r $(DIRECTORIES)