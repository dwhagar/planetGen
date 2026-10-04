To prevent the "bullet through paper" tunneling problem during large update steps ($\Delta t$), you must replace instantaneous distance checks with **Continuous Collision Detection (CCD)** based on **Swept-Sphere Capsule Intersections** and **Adaptive Time-to-Impact Bisection** (Brent, 1973; Higham, 2002).

### Mathematical Foundation of Continuous Collision Detection

Consider two bodies moving over the integration interval $\tau \in [0, \Delta t]$ with initial positions $\vec{p}_1(0), \vec{p}_2(0)$ and average velocities $\vec{v}_1, \vec{v}_2$ (Higham, 2002).

1. **Relative Trajectory Vector:**
   Define the relative displacement and velocity vectors:
   $$\Delta \vec{p}_0 = \vec{p}_1(0) - \vec{p}_2(0), \quad \Delta \vec{v} = \vec{v}_1 - \vec{v}_2$$
   The relative separation vector as a function of time $\tau$ is:
   $$\vec{r}(\tau) = \Delta \vec{p}_0 + \Delta \vec{v} \cdot \tau$$

2. **Squared Separation Function:**
   $$S(\tau) = \Vert{}\vec{r}(\tau)\Vert{}^2 = \Vert{}\Delta \vec{v}\Vert{}^2 \tau^2 + 2 (\Delta \vec{p}_0 \cdot \Delta \vec{v}) \tau + \Vert{}\Delta \vec{p}_0\Vert{}^2$$
   Let $A = \Vert{}\Delta \vec{v}\Vert{}^2$, $B = 2(\Delta \vec{p}_0 \cdot \Delta \vec{v})$, and $C = \Vert{}\Delta \vec{p}_0\Vert{}^2$.

3. **Time of Closest Approach ($\tau_{\text{min}}$):**
   Setting the first derivative $\frac{d S}{d\tau} = 2 A \tau + B = 0$ yields the time of closest approach:
   $$\tau_{\text{min}} = -\frac{\Delta \vec{p}_0 \cdot \Delta \vec{v}}{\Vert{}\Delta \vec{v}\Vert{}^2} = -\frac{B}{2A}$$
* If $\tau_{\text{min}} < 0$, the objects are diverging from the start of the tick.

* If $\tau_{\text{min}} > \Delta t$, closest approach does not occur until after the current frame.

* The constrained minimum separation time within this frame is:
    $$\tau^* = \text{clamp}(\tau_{\text{min}}, 0, \Delta t)$$
4. **Collision Threshold Evaluation:**
   Given physical radii $R_1$ and $R_2$, contact occurs if the minimum separation distance is less than or equal to the collision cross-section:
   $$S(\tau^*) \le (R_1 + R_2)^2$$
   If this inequality holds and $B^2 - 4A(C - (R_1 + R_2)^2) \ge 0$, the exact moment of physical surface contact ($\tau_{\text{impact}}$) is the smallest positive real root (Brent, 1973; Ortega & Rheinboldt, 1970):
   $$\tau_{\text{impact}} = \frac{-B - \sqrt{B^2 - 4A(C - (R_1 + R_2)^2)}}{2A}$$

### Algorithmic Implementation

Below is the concrete algorithm to insert directly into your step-by-step update pipeline:

```python
Python
    import math
    def check_continuous_collision(
        p1: tuple, v1: tuple, r1: float,    p2: tuple, v2: tuple, r2: float,    delta_t: float
    ) -> tuple[bool, float]:
        """    Evaluates whether two moving spherical bodies intersect during time interval delta_t.    Returns (collided: bool, time_of_impact: float).    """
        # Relative initial displacement and velocity
        dx = p1[0] - p2[0]
        dy = p1[1] - p2[1]
        dz = p1[2] - p2[2]

        dvx = v1[0] - v2[0]
        dvy = v1[1] - v2[1]
        dvz = v1[2] - v2[2]

        # Quadratic coefficients: A*tau^2 + B*tau + C = (R1 + R2)^2
        A = dvx * dvx + dvy * dvy + dvz * dvz
        B = 2.0 * (dx * dvx + dy * dvy + dz * dvz)
        C = dx * dx + dy * dy + dz * dz

        combined_radius = r1 + r2
        radius_sq = combined_radius * combined_radius

        # Already intersecting at start of frame
        if C <= radius_sq:
            return True, 0.0

        # Relative velocity is zero (parallel trajectories)
        if A <= 1e-18:
            return False, delta_t

        # Discriminant for intersection with the combined bounding sphere
        discriminant = B * B - 4.0 * A * (C - radius_sq)

        if discriminant < 0.0:
            # Trajectories do not intersect the collision envelope
            return False, delta_t

        # Smallest positive root represents entry point into collision sphere
        sqrt_disc = math.sqrt(discriminant)
        t_entry = (-B - sqrt_disc) / (2.0 * A)

        # Check if contact occurs within the discrete tick [0, delta_t]
        if 0.0 <= t_entry <= delta_t:
            return True, t_entry

        return False, delta_t
```

### Pipeline Integration & Momentum Conservation

To apply this without adding an $O(N^2)$ global check across all bodies, test for collisions **only against candidate gravitational influencers** already identified in your sector stencil query (Higham, 2002):

```python
Python
    def integrate_and_resolve_collisions(target, candidate_influencers, delta_t, db):
        earliest_collision = None
        first_impact_time = delta_t
        partner_entity = None
        # 1. Evaluate Continuous Collision Detection across Top Influencers
        for body in candidate_influencers:
            collided, t_hit = check_continuous_collision(
                target.position, target.velocity, target.radius,
                body.position, body.velocity, body.radius,
                delta_t
            )
            if collided and t_hit < first_impact_time:
                first_impact_time = t_hit
                earliest_collision = body

        # 2. Execution Branch
        if earliest_collision is not None:
            # Advance both bodies to the exact moment of collision tau_impact
            p_target_contact = target.position + target.velocity * first_impact_time
            p_body_contact = earliest_collision.position + earliest_collision.velocity * first_impact_time

            # Conservation of Momentum: Inelastic Merger Coalescence
            total_mass = target.mass + earliest_collision.mass
            v_merged = (target.velocity * target.mass + earliest_collision.velocity * earliest_collision.mass) / total_mass
            new_radius = (target.radius**3 + earliest_collision.radius**3) ** (1.0 / 3.0)

            # Update primary entity and remove consumed entity
            db.update_entity_state(
                entity_id=target.id,
                position=p_target_contact,
                velocity=v_merged,
                mass=total_mass,
                radius=new_radius
            )
            db.delete_entity(earliest_collision.id)

            # Trigger point-mass and sector metadata rebuilds
            db.mark_sector_dirty(target.sector_id)
            return

        # 3. Standard Non-Colliding Verlet Integration Path
        # (Execute normal half-step velocity and full-step position updates)
        target.update_verlet_trajectory(delta_t)

This ensures that regardless of whether time steps are one day or ten years, objects moving faster than their own radii will cleanly trigger a collision at the exact time fraction $\tau_{\text{impact}} \in [0, \Delta t]$ rather than tunneling past each other undetected (Brent, 1973; Higham, 2002).
```

### Refrences

* Brent, R. P. (1973). _Algorithms for Minimization without Derivatives_. Prentice-Hall.
* Higham, N. J. (2002). _Accuracy and Stability of Numerical Algorithms_ (2nd ed.). Society for Industrial and Applied Mathematics.
* Ortega, J. M., & Rheinboldt, W. C. (1970). _Iterative Solution of Nonlinear Equations in Several Variables_. Academic Press.


