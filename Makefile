CXX := g++
CXXFLAGS := -std=c++17 -Wall -Wextra -pedantic -I.
BUILD_DIR := build

.PHONY: test clean

test: $(BUILD_DIR)/test_cfg_parser
	./$(BUILD_DIR)/test_cfg_parser

$(BUILD_DIR)/test_cfg_parser: tests/test_cfg_parser.cpp cfg/CFGParser.cpp cfg/CFGParser.hpp cfg/CFG.hpp | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) tests/test_cfg_parser.cpp cfg/CFGParser.cpp -o $@

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

clean:
	rm -rf $(BUILD_DIR)
