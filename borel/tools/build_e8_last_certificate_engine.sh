#!/bin/sh
# Reproducible verification dependency; does not install anything globally.
set -eu
tools=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
dest=${E8_F4_BUILD_DIR:-"$HOME/.cache/sage-borel-e8-f4"}
if [ ! -d "$dest/.git" ]; then
 git clone --depth 1 --branch v0.4.9 https://github.com/algebraic-solving/msolve.git "$dest"
fi
cd "$dest"
test "$(git rev-parse HEAD)" = c547cb83ac3e876a2d72a2cc8ee191b543287d85
python3 - <<'PY'
from pathlib import Path
p=Path('src/neogb/f4.c');s=p.read_text()
needle='''      rrt1 = realtime();
      if (st->info_level > 1) {
        printf("%13.2f sec\\n", rrt1-rrt0);
      }
'''
patch='''      if (getenv("E8_F4_BOUND") && st->current_deg >= atoi(getenv("E8_F4_BOUND"))) {
        printf("E8 partial stop at degree %d\\n",st->current_deg); fflush(stdout); break;
      }
'''
if 'E8_F4_BOUND' not in s:
 assert s.count(needle)==1
 p.write_text(s.replace(needle,needle+patch))
else:
 assert patch in s
PY
./autogen.sh
./configure
make -j4
gcc -O2 -I"$dest/src/neogb" "$tools/e8_last_f4_driver.c" -L"$dest/src/neogb/.libs" -Wl,-rpath,"$dest/src/neogb/.libs" -lneogb -lgmp -o "$dest/e8-last-f4"
printf '\nCertificate engine: %s/e8-last-f4\n' "$dest"
