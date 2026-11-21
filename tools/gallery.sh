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
atwood|6|Atwood machine
blocks|6|Two blocks with pendulums
bowl|6|Circle rolling inside a circle
collision|10|Two-disc collision in lab and CM frames, momentum arrows
crank|6|Disc driving a crank and coupled sliders
cycloid|5|Rolling wheel: cycloid, curtate and prolate traces
doublependulum|20|Double pendulum from its Lagrangian, with a diverging twin
frame|6|Sliding frame with sprung blocks and a two-link arm
freebody|13|Torque-free asymmetric body: tumbling about the intermediate axis
halfcar|6|Half-car on a wavy road
hoop|6|Bead on a rotating hoop: energy changes, Jacobi integral does not
impact|8|Impact oscillator with a pendulum
incline|6|Cylinder and block on a wedge
kepler|7|Kepler orbit with equal-area wedges
modes|11|Coupled oscillators: normal modes and energy beats
motor|5|DC motor driving a flywheel
linkage|6|Parallelogram linkage with a rolling disc, slider and double pendulum
loop|5|Loop-the-loop
phase|14|Damped pendulum with a phase-portrait plot
projectile|5|Projectile with and without drag: traces and vectors
pulleys|6|Pulley system
race|16|Rolling down an incline: hoop, disc, sphere and a frictionless slider
resonance|20|Swept-frequency drive: resonance curve of a damped oscillator
ring|6|Ring of masses
spherical|14|Spherical pendulum: conserved L_z, shadow rosette on the floor
shaker|6|Shaker with a driven block, rolling disc and pendulum
swing|8|Swinging Atwood machine with energy plot
tower|8|Rotating tower with a sliding block and a pendulum (Lagrangian)
top|4|Heavy symmetric top: precession and nutation
waves|15|Standing waves on a bead chain: modes 1, 2 and 3
bjt|5|Common-emitter transistor amplifier
bridge|10|Full-wave bridge rectifier with smoothing capacitor
cyclotron|8|Charges in a magnetic field: circles and equal periods
exb|9|E x B drift of a charge from rest
faraday|12|Loop pulled through a magnetic field: flux, emf and Lenz braking
lc|5|LC tank: energy exchange between C and L
rc|5|RC charging and discharging, time constant
opamp|5|Inverting op-amp amplifier that saturates
rectifier|5|Half-wave rectifier
rlc|5|Series RLC driven by an AC source
transformer|5|Step-up transformer with dot convention
zener|5|Zener shunt regulator against supply ripple
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

