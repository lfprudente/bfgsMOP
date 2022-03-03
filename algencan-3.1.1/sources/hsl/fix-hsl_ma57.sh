#!/usr/bin/env bash

LOC="hsl_lma57d.f90"
TMP="hsl_lma57d.tmp"

cd $(dirname "$0")
cat "$1" > "$LOC"

LINE1=$(egrep -in "end\ +type\ +ma57_finfo" "$LOC" | cut -d: -f1)
LINE2=$(egrep -in "finfo%more\ +=\ +info\(2\)" "$LOC" | cut -d: -f1)

head -n $((LINE1-1)) "$LOC"                            > "$TMP"
echo '      real(wp) :: pivot'                        >> "$TMP"
tail -n +$((LINE1)) "$LOC" | head -n $((LINE2-LINE1)) >> "$TMP"
echo '      finfo%pivot = rinfo(20)'                  >> "$TMP"
tail -n +$((LINE2)) "$LOC"                            >> "$TMP"
cat "$TMP"                                             > "$LOC"

rm -f "$TMP"
exit 0
