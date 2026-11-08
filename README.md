# Physics Illustration library in C (pil-c)

A small C99 library that simulates 2D mechanical systems and electric circuits and
draws them the way a textbook does: grey bodies, coil springs, dashpots, hatched
ground, hinges, dimension lines, circuit symbols, italic serif labels with
`x_1`-style subscripts and Greek letters. You describe a problem by dropping
elements into a world; the library solves the dynamics and renders it to PNG /
PDF / SVG, to a GIF/MP4, or to an interactive window.

See all rendered animations in [gallery](https://harshit-dhanwalkar.github.io/pil-c/)

<table>
  <tr>
    <td align="center"><img src="https://harshit-dhanwalkar.github.io/pil-c/assets/halfcar.gif" width="250"><br><b>Half-car on a wavy road</b></td>
    <td align="center"><img src="https://harshit-dhanwalkar.github.io/pil-c/assets/loop.gif" width="250"><br><b>Loop-the-loop</b></td>
  </tr>
  <tr>
    <td align="center"><img src="https://harshit-dhanwalkar.github.io/pil-c/assets/impact.gif" width="250"><br><b>Impact oscillator with a pendulum</b></td>
    <td align="center"><img src="https://harshit-dhanwalkar.github.io/pil-c/assets/swing.gif" width="250"><br><b>Swinging Atwood machine</b></td>
  </tr>
  <tr>
    <td align="center"><img src="https://harshit-dhanwalkar.github.io/pil-c/assets/rlc.gif" width="250"><br><b>Series RLC, current and node voltages</b></td>
    <td align="center"><img src="https://harshit-dhanwalkar.github.io/pil-c/assets/motor.gif" width="250"><br><b>DC motor driving a flywheel</b></td>
  </tr>
</table>

## Features

- **Rigid body dynamics:** full 2D translation and rotation with impulse-based constraints.
- **Textbook renderer:** Cairo-powered backend that draws hatching, zig-zag springs, dashpots and $\LaTeX$-style labels.
- **Drop-in elements:** springs, dampers, pulleys, ropes, slots, gears/belts (`pil_add_link`) and contacts: floor lines, terrain, circular tracks (loops, bowls), impact stops with restitution.
- **Circuits:** resistors, capacitors, inductors, AC/DC sources, diodes, switches, motors and voice coils, solved by modified nodal analysis and coupled to the mechanics. Currents are animated as arrows followed by dots and node voltages appear at the joints while it runs.
- **Generalised coordinates:** give the kinetic and potential energy of a 3D problem and the library forms and integrates the Euler-Lagrange equations (`pil_lagrange.h`), with a small 3D projection for drawing.
- **Plots:** time-series panels next to the figure with a cursor that follows the animation and an energy bar.
- **Interactive:** built-in SDL2 window for live simulation (grab bodies with your mouse, pause, step and change gravity).
- **Export:** output directly to `.png`, `.svg`, or record `.gif`/`.mp4` animations via `ffmpeg`.

## Build

```sh
sudo apt install build-essential pkg-config libcairo2-dev   # required
sudo apt install libsdl2-dev ffmpeg                         # optional: window, gif/mp4
make            # library + examples + tests   (SDL2 is detected automatically)
make test       # numerical regression tests against analytic solutions
sudo make install        # headers, libpil.a, pil.pc
```

Use it from your own program with `cc prog.c $(pkg-config --cflags --libs pil)`.
Text uses _Latin Modern_ (Roman + Math) if installed, as it is with $\TeX$ Live.

## Try it

```sh
./build/halfcar   half.png            # render the first frame
./build/halfcar   --record half.gif 5 # simulate 5 s into a gif (mp4/webm also work)
./build/halfcar   --live              # interactive window

# other examples
./build/atwood      ./build/blocks     ./build/bowl      ./build/crank
./build/frame       ./build/incline    ./build/ring      ./build/pulleys
./build/swing       ./build/impact     ./build/loop
./build/rlc         ./build/rectifier  ./build/motor     # circuits

make gallery
```

Window keys: `space` pause · `r` reset · `→` step · `+`/`-` speed · `p`/`f`/`g` save
PNG/PDF/SVG · `h` HUD · `esc`. Mouse: left-drag pulls a body, wheel zooms,
right-drag pans.

## A complete program

```c
#include <pil/pil_core.h>
#include <pil/pil_draw.h>

int main(void)
{
    pil_world *w = pil_world_new();                       // gravity defaults to (0,-9.81)

    pil_add_ground(w, (pil_vec2){-1, 0}, (pil_vec2){1, 0});
    pil_add_ground(w, (pil_vec2){1, 1.8}, (pil_vec2){-1, 1.8});    // ceiling: hatching on the other side
    pil_body *m = pil_add_body(w, &(pil_body_desc){
        .shape = pil_rect(0.4, 0.25), .mass = 2.0, .pos = {0, 0.9}, .label = "m" });
    pil_add_prismatic(w, pil_on(m, 0, 0), pil_at(0, 0), (pil_vec2){0, 1});   // vertical guide

    pil_element *k = pil_add_spring(w, pil_at(0, 1.8), pil_on(m, 0, 0.125), 150, 0.6);
    pil_element_set_label(k, "k");
    pil_add_damper(w, pil_at(0.1, 1.8), pil_on(m, 0.1, 0.125), 3.0);

    pil_world_settle(w, 5.0);                             // optional: start at equilibrium
    pil_view v;
    pil_view_fit(&v, -1, -0.1, 1, 2, 600, 600, 20);
    pil_style st = pil_style_textbook();
    pil_record(w, &v, &st, "out.gif", 4.0, 30.0, 5e-4);   // or pil_export_png/pdf/svg, or pil_app_run
    pil_world_free(w);
}
```

## Concepts

- **Bodies** have a shape (`pil_rect`, `pil_disc`, `pil_point_mass`, `pil_poly`,
  `pil_arc_band`), mass, inertia (computed from the shape if you give 0), a state,
  and flags (`PIL_STATIC`, `PIL_SPOKE`, `PIL_HIDDEN`, ...). `pos` is the centre of mass.
  For polygons and arc bands the shape is re-centred on its centroid automatically.
- **Anchors** say _where_ an element attaches: `pil_on(body, lx, ly)` is a point fixed
  in a body, `pil_at(x, y)` is a point fixed in the world. Every element takes anchors.
- **Elements** are what you drop in. Force elements add forces; constraints remove
  degrees of freedom.

| you want                             | call                                                        |
| ------------------------------------ | ----------------------------------------------------------- |
| spring / damper                      | `pil_add_spring`, `pil_add_damper`                          |
| torsional spring (+damper)           | `pil_add_torsion_spring`                                    |
| constant force / "pull"              | `pil_add_force`                                             |
| any force law you like               | `pil_add_custom_force(fn)` - call `pil_apply_force` inside  |
| pin / hinge                          | `pil_add_revolute`                                          |
| rigid rod / rope                     | `pil_add_distance(..., 0)` / `(..., PIL_ROPE)`              |
| rope over a pulley                   | `pil_add_pulley_rope` (+ a static disc body for the pulley) |
| sliding guide, rail                  | `pil_add_slot`, `pil_add_prismatic` (slot + fixed angle)    |
| lock an angle / motor                | `pil_add_fixed_angle`, `pil_add_rotation_drive`             |
| roll without slipping on a line      | `pil_add_rolling` (discs, cylinders, arc bands)             |
| block/wedge/floor contact + friction | `pil_add_floor_contact` (line may belong to a moving body)  |
| uneven moving road                   | `pil_add_terrain` + `pil_add_terrain_contact`               |
| belt, gears, rack: two coordinates tied by a ratio | `pil_add_link` |
| disc on a circular track, leaves it when N = 0 | `pil_add_circle_contact` (+ `_set_arc`) |
| stop with a coefficient of restitution | `pil_add_stop` (lost energy in `pil_world_dissipated`) |
| ramp or ledge that ends | `pil_floor_set_range` |


- **Annotations** (visual only): `pil_add_dimension`, `_label`, `_arrow`, `_line`,
  `_angle_marker`, `_ground`, `_gravity_indicator`. For anything else register an overlay
  with `pil_world_add_overlay` and draw with Cairo plus the public symbol primitives
  (`pil_draw_spring`, `pil_draw_dashpot`, `pil_draw_hinge`, `pil_draw_arrow`, `pil_draw_text`, ...).
- **Text** understands `x_1`, `x_{12}`, `x^2` and UTF‑8 (`"\xcf\x86"` $\rightarrow$ φ).
- **Plots and energy:** `pil_add_plot` + `pil_plot_add_signal` put time-series panels
- beside the figure (curves are precomputed, a cursor follows the animation);
- `pil_add_energy_bar` draws the E0 / -ΔE / T / V column used for impacts.
- **Text** understands `x_1`, `x_{12}`, `x^2` and UTF‑8 (`"\xcf\x86"` $\rightarrow$ φ).

## Circuits

```c
#include <pil/pil_circuit.h>

pil_circuit *c = pil_add_circuit(w);
pil_node a = pil_circuit_node(c, 0, 1), b = pil_circuit_node(c, 2, 1), g = pil_circuit_node(c, 0, 0);
pil_circuit_ground(c, g);
pil_vsource_dc(c, a, g, 5.0, "V");
pil_resistor(c, a, b, 1e3, "R");
pil_capacitor(c, b, g, 1e-6, 0.0, "C");
```

Nodes sit at points of the figure and `pil_wire` joins them. The circuit is stepped with the
world (also when there are no bodies), so a motor or voice coil (`pil_motor`, `pil_voice_coil`)
couples current to torque or force and back-EMF to speed. Read results with
`pil_celement_voltage/current/power`, `pil_circuit_node_voltage`, and check Kirchhoff's current law
with `pil_circuit_node_current_sum`.

While a circuit runs, every branch shows its current as an arrow followed by dots that stream with
the current (faster for the strongest branch, reversing with the sign), and each joint is labelled with
its voltage. Style knobs: `show_current`, `current_dot_speed`, `current_dot_r`, `current_color`
(light grey by default), `show_voltages`. In the examples `PIL_EXAMPLE_TIME_SCALE` slows the video
down (circuits live on millisecond scales) and `PIL_EXAMPLE_IMAGE_TIME` advances the still image so
it shows currents.

## 3D and generalised coordinates

Some figures (a rotating tower with a sliding block and a pendulum, a spinning body) are easier to
state as energies than as constraints:

```c
static double T(const double *q, const double *qd, void *u) { /* kinetic energy */ }
static double V(const double *q, void *u)                   { /* potential energy */ }

pil_lsys *s = pil_add_lagrangian(w, 2, T, V, &params);
pil_lsys_set_state(s, (double[]){0.9, 0.7}, (double[]){0.3, -0.5});
```

The Euler-Lagrange equations are formed numerically and integrated along with the world, so plots,
recording and reset work as usual. `T` may contain terms linear in the velocities (a rotating frame),
which gives the Coriolis and gyroscopic terms; `pil_lsys_jacobi` returns the conserved Jacobi
integral. Draw with an overlay and `pil_proj3_point`, an orthographic z-up projection.



## How it works

Each body carries its own coordinates $(x, y, \theta)$ and every joint is an equation
$\phi(q, t) = 0$. With $J = \partial \phi/\partial q$, the motion is $M a + J^{T} \lambda = F$, $J a = \gamma$.
Eliminating $a$ leaves a small dense system for the multipliers $\lambda$. $\gamma$ includes the
velocity-dependent terms and Baumgarte stabilisation, and the state is advanced
with RK4. A new constraint therefore only has to supply $\phi$, $J$ and $\gamma$ - no equations of motion are derived per problem. Ropes add a one-line active set (tension $\lambda \ge 0$, otherwise the rope is slack).

## Tests

`make test` checks the engine against closed-form results: spring period, physical
pendulum period, rolling at `g sin α / 1.5`, Atwood acceleration and tension (to
machine precision), Coulomb deceleration, frictionless slide, double-pendulum energy
drift ($< 5\cdot 10^{-4}$ of mg.L over 10 s), rope slack/taut, and mouse picking.
The circuit solver is checked against RC/RL/RLC closed forms, energy in LC, AC gain,
a diode rectifier, a switch, a motor with flywheel, and Kirchhoff's laws on a resistor
network; impact stops against `e·v` and the lost energy; the circular track against the
analytic detachment angle ($\cos \phi = −1/3$) and rolling speed; the Lagrangian solver against
a pendulum, the centrifugal equilibrium of the rotating arm and the Jacobi integral of
the gyroscopic tower. `make check` runs every example for a few seconds and fails on
constraint drift or NaN.

## Adding your own element

For a one-off, `pil_add_custom_force` needs no changes to the library. For a reusable
element: add a kind to `pil_el_kind` (`src/pil_internal.h`), a creator in `elements.c`,
a case in `pil_element_apply` (force) or `pil_constraint_eval` (fill `phi`, `vel`, `ja`,
`jb`, and `gamma` from the derivation `Φ̈ = J a − γ₀`), a case in `draw_element`, and a
test in `tests/test_core.c`. `EL_PULLEY` is a compact worked example of each step.

# TODO

- [x] Electrodynamics and circuits
- [x] Impacts, circular tracks, generalised coordinates
- [ ] Make controller `w/a/s/d` to control actuator
- [ ] Spinning 3D bodies: the Lagrangian solver covers problems stated as energies; a rigid-body 3D renderer (quaternions) is still missing
- [ ] Rolling is implemented on straight lines and circular tracks. Two curved bodies rolling on each other are not yet available.
- [ ] No body–body collision. Contact exists between a body and a line (`pil_add_floor_contact`) or a circular track (`pil_add_circle_contact`).
- [ ] Contacts are penalty springs: use `max_dt ≲ 5·10⁻⁴` with stiff contacts (default in the examples).
- [ ] The circuit solver assumes the wire graph is a forest when it assigns currents to wires; wire loops carry zero.
- [ ] Mouse interaction in the window is covered by a unit test of the pick/grab physics and a
      headless run of the window; to be exercised with a real pointer

## License

This project is licensed under the [Apache License 2.0](http://www.apache.org/licenses/LICENSE-2.0).

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
