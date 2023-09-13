#!/bin/bash
[[ $# -ne 1 ]] && echo Usage: $0 [TSV_FN] && exit -1

TSV_FN=$1

echo "<table>"
head -n 1 $TSV_FN | \
    sed -e 's/^/<tr><th>/' -e 's/[\t]/<\/th><th>/g' -e 's/$/<\/th><\/tr>/'
tail -n +2 $TSV_FN | \
    sed -e 's/^/<tr><td>/' -e 's/[\t]/<\/td><td>/g' -e 's/$/<\/td><\/tr>/'
echo "</table>"

