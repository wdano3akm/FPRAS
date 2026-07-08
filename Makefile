CXX := g++
CXXFLAGS := -std=c++17 -Wall -Wextra -pedantic -I.
BUILD_DIR := build

.PHONY: test clean

test: $(BUILD_DIR)/test_cfg_parser $(BUILD_DIR)/test_ptfromcfg $(BUILD_DIR)/test_algorithm
	./$(BUILD_DIR)/test_cfg_parser
	./$(BUILD_DIR)/test_ptfromcfg
	./$(BUILD_DIR)/test_algorithm

$(BUILD_DIR)/test_cfg_parser: tests/test_cfg_parser.cpp cfg/CFGParser.cpp cfg/CFGParser.hpp cfg/CFG.hpp | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) tests/test_cfg_parser.cpp cfg/CFGParser.cpp -o $@

$(BUILD_DIR)/test_ptfromcfg: tests/test_ptfromcfg.cpp tests/cyk.cpp tests/cyk.hpp plustimes/ptfromcfg.cpp plustimes/plustimes.hpp cfg/CFG.hpp | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) tests/test_ptfromcfg.cpp tests/cyk.cpp plustimes/ptfromcfg.cpp -o $@

$(BUILD_DIR)/test_algorithm: tests/test_algorithm.cpp algorithm/algorithm.cpp algorithm/algorithm.hpp plustimes/plustimes.hpp | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) tests/test_algorithm.cpp -o $@

$(BUILD_DIR)/main: main.cpp cfg/CFGParser.cpp cfg/CFGParser.hpp cfg/CFG.hpp plustimes/ptfromcfg.cpp plustimes/plustimes.hpp algorithm/algorithm.cpp algorithm/algorithm.hpp | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) main.cpp cfg/CFGParser.cpp plustimes/ptfromcfg.cpp algorithm/algorithm.cpp -o $@

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

clean:
	rm -rf $(BUILD_DIR)
