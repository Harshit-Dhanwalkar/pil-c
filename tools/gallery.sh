#!/bin/sh
# Render every example into site/: an mp4 + gif + still per figure and
# an index.html gallery
# Needs ffmpeg and a built tree (make)
# Run: make gallery

set -e

OUT=${1:-site}
BIN=${BIN:-build}
REPO=${REPO:-https://github.com/Harshit-Dhanwalkar/pil-c}
mkdir -p "$OUT/assets"

# name | seconds | title
LIST='
halfcar|6|Half-car on a wavy road
incline|6|Cylinder and block on a wedge
atwood|6|Atwood machine
blocks|6|Two blocks with pendulums
bowl|6|Circle rolling inside a circle
crank|6|Disc driving a crank and coupled sliders
frame|6|Sliding frame with sprung blocks and a two-link arm
ring|6|Ring of masses
pulleys|6|Pulley system
swing|8|Swinging Atwood machine with energy plot
impact|8|Impact oscillator with a pendulum
loop|5|Loop-the-loop
rlc|5|Series RLC driven by an AC source
rectifier|5|Half-wave rectifier
motor|5|DC motor driving a flywheel
'

CARDS=
for row in $(echo "$LIST" | tr ' ' '_'); do
  [ -n "$row" ] || continue
  name=${row%%|*}; rest=${row#*|}; secs=${rest%%|*}; title=$(echo "${rest#*|}" | tr '_' ' ')
  [ -x "$BIN/$name" ] || { echo "skip $name (not built)"; continue; }
  echo "== $name"
  "$BIN/$name" "$OUT/assets/$name.png" >/dev/null
  "$BIN/$name" --record "$OUT/assets/$name.mp4" "$secs" >/dev/null
  # gif: 20 fps, 560 px wide, one palette per figure keeps it small
  ffmpeg -loglevel error -y -i "$OUT/assets/$name.mp4" \
    -vf "fps=20,scale=560:-1:flags=lanczos,split[a][b];[a]palettegen=max_colors=64[p];[b][p]paletteuse=dither=bayer" \
    "$OUT/assets/$name.gif"
  CARDS="$CARDS
  <article>
    <video autoplay loop muted playsinline poster=\"assets/$name.png\"><source src=\"assets/$name.mp4\" type=\"video/mp4\"></video>
    <h2>$title</h2>
    <p><a href=\"$REPO/blob/main/examples/$name.c\">$name.c</a> &middot; <a href=\"assets/$name.gif\">gif</a> &middot; <a href=\"assets/$name.mp4\">mp4</a></p>
  </article>"
done

cat > "$OUT/index.html" <<HTML
<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>pil-c: Physics Illustration Library</title>
<style>
  :root { --bg1:#0f1b2d; --bg2:#1d3557; --card:#fff; --ink:#1a1a1a; --mute:#a9b8d0; }
  * { box-sizing: border-box; }
  body { margin:0; min-height:100vh; color:#fff; font-family: "Latin Modern Roman", Georgia, serif;
         background: radial-gradient(1200px 600px at 20% -10%, var(--bg2), var(--bg1)) fixed; }
  header { max-width: 980px; margin: 0 auto; padding: 48px 20px 8px; }
  h1 { margin:0 0 8px; font-size: 2.2rem; font-weight: 600; }
  header p { margin:0; color: var(--mute); max-width: 62ch; line-height: 1.5; }
  header a { color:#9cc3ff; }
  main { max-width: 1100px; margin: 0 auto; padding: 24px 20px 64px;
         display:grid; grid-template-columns: repeat(auto-fill, minmax(320px, 1fr)); gap: 22px; }
  article { background: rgba(255,255,255,.06); border:1px solid rgba(255,255,255,.12);
            border-radius: 14px; padding: 12px; backdrop-filter: blur(2px); }
  video { width:100%; display:block; border-radius: 8px; background: var(--card); }
  h2 { font-size: 1.05rem; margin: 12px 4px 2px; font-weight: 600; }
  article p { margin: 0 4px 4px; font-size:.9rem; color: var(--mute); }
  article a { color:#9cc3ff; text-decoration:none; }
  article a:hover { text-decoration: underline; }
  footer { text-align:center; color: var(--mute); padding: 0 0 36px; font-size:.9rem; }
</style>
</head>
<body>
<header>
  <h1>pil-c</h1>
  <p>A small C library that simulates mechanical systems and circuits and draws them the way a
  textbook does. Every figure below is an example in the
  <a href="$REPO">repository</a>, rendered by CI on each push.</p>
</header>
<main>$CARDS
</main>
<footer>built $(date -u +%Y-%m-%d) from <a href="$REPO" style="color:#9cc3ff">Harshit-Dhanwalkar/pil-c</a></footer>
</body>
</html>
HTML
echo "wrote $OUT/index.html"

