#!/bin/bash


grep -E "(s__)|(^ID)" $1 | grep -v "t__" | sed 's/^.*s__//g' > $2
