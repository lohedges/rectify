# Rectify makefile
#
# Author : Lester. O. Hedges
# Email  : lester.hedges+rectify@gmail.com
# Date   : January 21st 2013

# Common Makefile configurations

# C++ compiler
CXX=g++

# Flags for the C++ compiler
CXXFLAGS=-O3

# Installation directory
PREFIX=/usr/local

# Install command
INSTALL=install

# Flags for install command for executable
IFLAGS_EXEC=-m 0755

# Flags for install command for non-executable files
IFLAGS=-m 0644
