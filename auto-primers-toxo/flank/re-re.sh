#!/bin/sh

./redump.pl slic/* | sort | uniq > re_avoid.txt;
echo "PacI" >> re_avoid.txt
