#!/bin/bash

awk '{print $1"\t"int(($2+$3)/2)"\t"int(($2+$3)/2)"\t"$4}' $1 > $1.centers
