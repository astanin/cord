#!/bin/sh

cat log | grep xsize | awk '{print($2,$(NF-3),$(NF-2),$(NF-1),$(NF));}' | tr -d 'a-df-z=' | sed 's/ e/ /g' > size.dat

