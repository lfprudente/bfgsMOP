#!/usr/bin/env bash

LOC="lma57ad.f"
TMP="lma57ad.tmp"

cd $(dirname "$0")
cat "$1" > "$LOC"

BEGIN=$(egrep -in 'SUBROUTINE\ +MA57OD' "$LOC" | cut -d: -f1)

for LINE in $(egrep -in 'INFO\(1\)\ +=\ +-5' "$LOC" | cut -d: -f1); do
  if [ $LINE -gt $BEGIN ]; then
    head -n  $((LINE+1)) "$LOC"     > "$TMP"
    echo '      RINFO(20) = PIVOT' >> "$TMP"
    tail -n +$((LINE+2)) "$LOC"    >> "$TMP"
    cat "$TMP"                      > "$LOC"
    break
  fi
done

for LINE in $(egrep -in 'INFO\(1\)\ +=\ +-6' "$LOC" | cut -d: -f1); do
  if [ $LINE -gt $BEGIN ]; then
    head -n  $((LINE+1)) "$LOC"     > "$TMP"
    echo '      RINFO(20) = PIVOT' >> "$TMP"
    tail -n +$((LINE+2)) "$LOC"    >> "$TMP"
    cat "$TMP"                      > "$LOC"
    break
  fi
done

for LINE in $(egrep -in 'INFO\(1\)\ +=\ +4' "$LOC" | cut -d: -f1); do
  if [ $LINE -gt $BEGIN ]; then
    head -n $LINE "$LOC"              > "$TMP"
    echo '        RINFO(20) = PIVOT' >> "$TMP"
    tail -n +$((LINE+1)) "$LOC"      >> "$TMP"
    cat "$TMP"                        > "$LOC"
    break
  fi
done

rm -f "$TMP"
exit 0
