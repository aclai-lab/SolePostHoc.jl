CXX = g++
CXXFLAGS = -O2 -std=c++11 -static-libstdc++ -static-libgcc

SRCS = main.cpp BornAgainDecisionTree.cpp FSpace.cpp RandomForest.cpp
OBJS = $(SRCS:.cpp=.o)
TARGET = bornAgain

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS) $(TARGET)