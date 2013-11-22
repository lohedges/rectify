# Rectify makefile
#
# Author : Lester. O. Hedges
# Email  : lester.hedges+rectify@gmail.com
# Date   : January 21st 2013

# Load the common configuration file
include ../config.mk

all: rectify

rectify: lodepng.cpp rectify.cpp
	$(CXX) $(CXXFLAGS) lodepng.cpp rectify.cpp -o rectify

clean:
	rm -f rectify

.PHONY: all
