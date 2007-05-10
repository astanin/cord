#!/bin/sh
cat | grep -v ^# | perl -p -e 's/^[ \t]*\n//' 
