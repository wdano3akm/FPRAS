CXX := g++
CXXFLAGS := -std=c++17 -Wall -Wextra -pedantic -pthread -I.
OPT_LEVEL ?= 3
BUILD_DIR := build

.PHONY: test optimized clean

test: $(BUILD_DIR)/test_cfg_parser $(BUILD_DIR)/test_dnnf_parser $(BUILD_DIR)/test_dnnf_ops $(BUILD_DIR)/test_ptfromcfg $(BUILD_DIR)/test_algorithm
	./$(BUILD_DIR)/test_cfg_parser
	./$(BUILD_DIR)/test_dnnf_parser
	./$(BUILD_DIR)/test_dnnf_ops
	./$(BUILD_DIR)/test_ptfromcfg
	./$(BUILD_DIR)/test_algorithm

$(BUILD_DIR)/test_cfg_parser: tests/test_cfg_parser.cpp cfg/CFGParser.cpp cfg/CFGParser.hpp cfg/CFG.hpp | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) tests/test_cfg_parser.cpp cfg/CFGParser.cpp -o $@

$(BUILD_DIR)/test_dnnf_parser: tests/test_dnnf_parser.cpp dnnf/DNNFParser.cpp dnnf/DNNFParser.hpp dnnf/DNNF.hpp | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) tests/test_dnnf_parser.cpp dnnf/DNNFParser.cpp -o $@

$(BUILD_DIR)/test_dnnf_ops: tests/test_dnnf_ops.cpp utils/dnnf_ops.cpp utils/dnnf_ops.hpp dnnf/DNNFParser.cpp dnnf/DNNFParser.hpp dnnf/DNNF.hpp | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) tests/test_dnnf_ops.cpp utils/dnnf_ops.cpp dnnf/DNNFParser.cpp -o $@

$(BUILD_DIR)/test_ptfromcfg: tests/test_ptfromcfg.cpp tests/cyk.cpp tests/cyk.hpp plustimes/ptfromcfg.cpp plustimes/plustimes.cpp plustimes/plustimes.hpp cfg/CFG.hpp | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) tests/test_ptfromcfg.cpp tests/cyk.cpp plustimes/ptfromcfg.cpp plustimes/plustimes.cpp -o $@

$(BUILD_DIR)/test_algorithm: tests/test_algorithm.cpp algorithm/algorithm.cpp algorithm/algorithm.hpp cfg/CFGParser.cpp cfg/CFGParser.hpp cfg/CFG.hpp plustimes/ptfromcfg.cpp plustimes/plustimes.cpp plustimes/plustimes.hpp utils/alg_checks.cpp utils/alg_checks.hpp tests/fixtures/randomized_tiny.cfg tests/fixtures/randomized_theorem.cfg | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) tests/test_algorithm.cpp cfg/CFGParser.cpp plustimes/ptfromcfg.cpp plustimes/plustimes.cpp utils/alg_checks.cpp -o $@

$(BUILD_DIR)/main: main.cpp cfg/CFGParser.cpp cfg/CFGParser.hpp cfg/CFG.hpp dnnf/DNNFParser.cpp dnnf/DNNFParser.hpp dnnf/DNNF.hpp plustimes/ptfromcfg.cpp plustimes/ptfromdnnf.cpp plustimes/ptfromdnnf.hpp plustimes/plustimes.cpp plustimes/plustimes.hpp algorithm/algorithm.cpp algorithm/algorithm.hpp utils/alg_checks.cpp utils/alg_checks.hpp utils/dnnf_ops.cpp utils/dnnf_ops.hpp | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) main.cpp cfg/CFGParser.cpp dnnf/DNNFParser.cpp plustimes/ptfromcfg.cpp plustimes/ptfromdnnf.cpp plustimes/plustimes.cpp algorithm/algorithm.cpp utils/alg_checks.cpp utils/dnnf_ops.cpp -o $@

optimized: $(BUILD_DIR)/main-O$(OPT_LEVEL)

$(BUILD_DIR)/main-O%: main.cpp cfg/CFGParser.cpp cfg/CFGParser.hpp cfg/CFG.hpp dnnf/DNNFParser.cpp dnnf/DNNFParser.hpp dnnf/DNNF.hpp plustimes/ptfromcfg.cpp plustimes/ptfromdnnf.cpp plustimes/ptfromdnnf.hpp plustimes/plustimes.cpp plustimes/plustimes.hpp algorithm/algorithm.cpp algorithm/algorithm.hpp utils/alg_checks.cpp utils/alg_checks.hpp utils/dnnf_ops.cpp utils/dnnf_ops.hpp | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -O$* main.cpp cfg/CFGParser.cpp dnnf/DNNFParser.cpp plustimes/ptfromcfg.cpp plustimes/ptfromdnnf.cpp plustimes/plustimes.cpp algorithm/algorithm.cpp utils/alg_checks.cpp utils/dnnf_ops.cpp -o $@

$(BUILD_DIR)/maing: main.cpp cfg/CFGParser.cpp cfg/CFGParser.hpp cfg/CFG.hpp dnnf/DNNFParser.cpp dnnf/DNNFParser.hpp dnnf/DNNF.hpp plustimes/ptfromcfg.cpp plustimes/ptfromdnnf.cpp plustimes/ptfromdnnf.hpp plustimes/plustimes.cpp plustimes/plustimes.hpp algorithm/algorithm.cpp algorithm/algorithm.hpp utils/alg_checks.cpp utils/alg_checks.hpp utils/dnnf_ops.cpp utils/dnnf_ops.hpp | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) main.cpp cfg/CFGParser.cpp dnnf/DNNFParser.cpp plustimes/ptfromcfg.cpp plustimes/ptfromdnnf.cpp plustimes/plustimes.cpp algorithm/algorithm.cpp utils/alg_checks.cpp utils/dnnf_ops.cpp -g -o $@

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

clean:
	rm -rf $(BUILD_DIR)
