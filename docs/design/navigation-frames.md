# Navigation reference frames

Boss's design for nested navigation frames, recorded 2026-09-30 (TODO
items 33 and 34). Implemented in `stellarObjects/navigation.py`:
`course_between` is Boss's `compute_course` below, `format_course` writes
"000 mark 000", and `warp_speed_c`/`fold_speed_c` are the travel speeds.

Decisions Boss approved on 2026-09-30:

- Mark is the elevation mod 360: 000-090 is up, 270-359 is down (270 is
  straight down), and nothing between 091 and 269 appears.
- 0 mark 0 points toward the current frame's center (the galactic core
  between sectors, the sector's center inside one).
- The quasar sits at the galactic center, so it can't set the zero
  meridian; +X stays the galaxy's existing axis (ring slot 0). It is only
  used as North's fallback when the ship sits on the frame's up axis.

How NAV picks the frame (`queryDb.nav_between`): two systems in the same
sector use the Sector Local Frame, centered on the sector's center (the
origin of sector-local positions); anything else uses the Galactic
Standard Frame. East is North x Up, as below, so with North toward the
core and Up galactic north, East runs clockwise seen from above. NAV
endpoints are whole systems and phenomena, so a course always leaves the
heliopause and the System Local Frame is never chosen yet; `course_between`
takes a center and up vector for when in-system navigation exists.

The design below is Boss's text as given.

Course projections are written as `0-359 mark 0-359`, with 0 mark 0
pointing toward the galactic core.

---

To revamp the coordinate hierarchy where "North" dynamically points
toward the local dominant center of mass (the central star in a system,
the sector core in a sector, and the galactic nucleus between sectors),
the navigation architecture requires a **nested hierarchical reference
frame**.

Instead of rotating absolute space, every local frame is defined as a
rigid transformation (translation + rotation) derived from an absolute
Cartesian standard.

## Frame hierarchy and reference planes

1. **Galactic Standard Frame (GSF):**
   - **Origin ($O_{\text{gal}}$):** Galactic Core (Sagittarius A*), $(0, 0, 0)$.
   - **Plane:** Invariable plane of the galactic disk.
   - **Orientation:** $+Z_{\text{gal}}$ points toward Galactic North
     (perpendicular to disk), $+X_{\text{gal}}$ points along an arbitrary
     zero-meridian (e.g., standard baseline toward a reference quasar),
     and $+Y_{\text{gal}} = Z_{\text{gal}} \times X_{\text{gal}}$.
   - **Inter-Sector Course "North":** Bearing $000$ points directly at
     $O_{\text{gal}}$.
2. **Sector Local Frame (SLF):**
   - **Origin ($O_{\text{sec}}$):** Barycenter of a cubic/spherical
     sector cell ($\sim 4\text{ pc}$ on a side or arc length).
   - **Orientation ("Sector North"):** $+X_{\text{sec}}$ points directly
     from the ship toward $O_{\text{sec}}$.
   - **Plane:** Aligned parallel to the galactic plane
     ($+Z_{\text{sec}} = +Z_{\text{gal}}$).
3. **System Local Frame (SysLF):**
   - **Origin ($O_{\text{sys}}$):** Central Star / Stellar Barycenter.
   - **Orientation ("System North"):** At any point $\vec{P}$, the local
     horizontal "North" vector points directly from $\vec{P}$ inward
     toward $O_{\text{sys}}$.
   - **Plane:** System invariable (ecliptic) plane, with
     $+Z_{\text{sys}}$ aligned with the star's net angular momentum
     vector.

## Mathematical model

Every navigation bearing consists of a vector $\vec{D}$ pointing toward
the destination:

$$\vec{D} = \vec{P}_{\text{target}} - \vec{P}_{\text{ship}}$$

To compute `Bearing (Yaw)` and `Mark (Pitch)` in any frame:

1. Define the local orthonormal basis $\{\hat{N}, \hat{E}, \hat{U}\}$:
   - **Up ($\hat{U}$):** Unit normal to the operational reference plane
     (e.g., stellar spin axis $\hat{\omega}_{\text{star}}$ or Galactic
     $+Z$).
   - **Inward Radial ($\hat{R}_{\text{in}}$):** Directed toward the
     center $O_{\text{center}}$:

     $$\vec{r} = O_{\text{center}} - \vec{P}_{\text{ship}}, \quad \vec{r}_{xy} = \vec{r} - (\vec{r} \cdot \hat{U})\hat{U}, \quad \hat{N} = \frac{\vec{r}_{xy}}{\Vert\vec{r}_{xy}\Vert}$$

     $\hat{N}$ is the instantaneous local horizontal North.
   - **East ($\hat{E}$):** Orthogonal horizontal vector matching the
     right-hand rule: $\hat{E} = \hat{N} \times \hat{U}$
2. Project target displacement $\vec{D}$ into this local basis:

   $$D_N = \vec{D} \cdot \hat{N}, \quad D_E = \vec{D} \cdot \hat{E}, \quad D_U = \vec{D} \cdot \hat{U}$$
3. Calculate spherical angles:
   - **Horizontal Bearing ($\theta$):** $\theta = \text{atan2}(D_E, D_N)$,
     normalized to $[0^\circ, 360^\circ)$ where $000^\circ$ is North
     toward the center and $090^\circ$ is East.
   - **Vertical Angle / Mark ($\phi$):**
     $\phi = \text{atan2}\left(D_U, \sqrt{D_N^2 + D_E^2}\right)$, ranging
     from $-90^\circ$ straight down to $+90^\circ$ straight up, or
     expressed as $[0^\circ, 360^\circ)$ where $270^\circ = -90^\circ$.

## Implementation pseudocode

```python
import math

class Vector3:
    def __init__(self, x, y, z):
        self.x, self.y, self.z = float(x), float(y), float(z)

    def __sub__(self, other):
        return Vector3(self.x - other.x, self.y - other.y, self.z - other.z)

    def dot(self, other):
        return self.x * other.x + self.y * other.y + self.z * other.z

    def cross(self, other):
        return Vector3(
            self.y * other.z - self.z * other.y,
            self.z * other.x - self.x * other.z,
            self.x * other.y - self.y * other.x
        )

    def magnitude(self):
        return math.sqrt(self.dot(self))

    def normalize(self):
        mag = self.magnitude()
        if mag == 0:
            return Vector3(0, 0, 0)
        return Vector3(self.x / mag, self.y / mag, self.z / mag)


class CelestialFrame:
    SYSTEM = "SYSTEM"
    SECTOR = "SECTOR"
    GALACTIC = "GALACTIC"


def compute_course(ship_pos_gal, target_pos_gal, frame_type, frame_center_gal, plane_up_gal=Vector3(0, 0, 1)):
    """
    Computes (bearing, mark, distance) where:
    - Bearing 000 is horizontal 'North' (toward frame_center_gal)
    - Bearing 090 is 'East'
    - Mark is vertical elevation (+90 = straight up, -90 = straight down)
    All input coordinates are absolute Galactic Standard Frame (GSF) vectors.
    """
    # 1. Target vector relative to ship
    displacement = target_pos_gal - ship_pos_gal
    distance = displacement.magnitude()
    if distance == 0:
        return {"bearing": 0.0, "mark": 0.0, "distance": 0.0}

    # 2. Establish Up Vector (normal to reference plane)
    u_hat = plane_up_gal.normalize()

    # 3. Vector from ship toward the frame center (Star, Sector Core, or Galactic Core)
    to_center = frame_center_gal - ship_pos_gal

    # Flatten 'to_center' against the reference plane to establish Horizontal North
    radial_proj = to_center - Vector3(u_hat.x * to_center.dot(u_hat),
                                     u_hat.y * to_center.dot(u_hat),
                                     u_hat.z * to_center.dot(u_hat))

    # Singularity handling: directly over center pole
    if radial_proj.magnitude() < 1e-9:
        # Fallback arbitrary reference vector if ship is directly on the z-axis
        fallback = Vector3(1, 0, 0)
        radial_proj = fallback - Vector3(u_hat.x * fallback.dot(u_hat),
                                         u_hat.y * fallback.dot(u_hat),
                                         u_hat.z * fallback.dot(u_hat))

    n_hat = radial_proj.normalize()      # North (Inward toward center)
    e_hat = n_hat.cross(u_hat).normalize() # East (Right-handed horizontal)

    # 4. Project displacement vector into (North, East, Up) local basis
    d_north = displacement.dot(n_hat)
    d_east  = displacement.dot(e_hat)
    d_up    = displacement.dot(u_hat)

    # 5. Calculate Bearing (Yaw) and Mark (Pitch)
    horizontal_dist = math.sqrt(d_north**2 + d_east**2)

    # Heading around plane: 0 deg = North, 90 deg = East
    bearing_rad = math.atan2(d_east, d_north)
    bearing_deg = math.degrees(bearing_rad) % 360.0

    # Mark: Elevation above/below plane (-90 to +90)
    mark_rad = math.atan2(d_up, horizontal_dist)
    mark_deg = math.degrees(mark_rad)

    return {
        "frame": frame_type,
        "bearing": round(bearing_deg, 6),
        "mark": round(mark_deg, 6),
        "distance": distance
    }
```

## Frame transition rules

- **System Departure to Sector Space:** When a vessel exits the
  heliopause ($\approx 120\text{ AU}$), the navigational computer
  automatically hands off the reference center from the primary star to
  the sector barycenter ($O_{\text{sec}}$).
- **Cross-Sector Transit:** When navigating across sector boundaries
  ($>4\text{ pc}$), the navigation console toggles to Galactic Standard
  Frame, setting North toward Sagittarius A* ($0, 0, 0$).

## Travel speeds

- Warp factor $w$: speed in c
  $= w^{10/3} + \frac{1}{1 + e^{-9.3575(w - 9.5)}} \left( \frac{198.9}{(10 - w)^{0.75}} + 1721.7 - w^{10/3} \right)$
- Dimensional fold factor $F$: speed in c $= \frac{6F^4}{10 - F}$

Every coefficient is a named constant in `program_constants` (`WARP_*`
and `FOLD_*`), and `src/tests/test_navigation.py` pins these values
(1 ly per 365.25 days at 1c; 1 kpc = 3,261.56 ly):

| Warp | Speed (c) | ly/day | Days per ly | Days per kpc |
|---:|---:|---:|---:|---:|
| 1 | 1.0 | 0.003 | 365.25 | 1,191,286 |
| 2 | 10.1 | 0.028 | 36.24 | 118,191 |
| 4 | 101.6 | 0.278 | 3.60 | 11,726 |
| 8 | 1,024.0 | 2.804 | 0.357 | 1,163 |
| 9 | 1,520.1 | 4.162 | 0.240 | 784 |
| 9.5 | 1,936.0 | 5.301 | 0.189 | 615 |
| 9.9 | 2,822.7 | 7.728 | 0.129 | 422 |
| 9.995 | 12,201.9 | 33.41 | 0.030 | 98 |

| Fold | Speed (c) | ly/day | Days per ly | Days per kpc |
|---:|---:|---:|---:|---:|
| 4 | 256.0 | 0.701 | 1.43 | 4,653 |
| 5 | 750.0 | 2.053 | 0.487 | 1,588 |
| 6 | 1,944.0 | 5.322 | 0.188 | 613 |
| 6.5 | 3,060.1 | 8.378 | 0.119 | 389 |
| 7 | 4,802.0 | 13.15 | 0.076 | 248 |
| 7.5 | 7,593.8 | 20.79 | 0.048 | 157 |
| 8 | 12,288.0 | 33.64 | 0.030 | 97 |
| 8.5 | 20,880.2 | 57.17 | 0.017 | 57 |
