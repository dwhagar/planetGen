### Simulation Architecture & Data Structures

```
CONSTANT G = 6.67430e-11                  # Gravitational constant (m^3 kg^-1 s^-2)[cite: 1]
CONSTANT LIGHT_YEAR_METERS = 9.461e15      # 1 Light Year in meters[cite: 1]
CONSTANT DELTA_R_LY = 11.5                # Sector radial thickness[cite: 1]
CONSTANT DELTA_Z_LY = 11.5                # Sector vertical thickness[cite: 1]
CONSTANT TARGET_ARC_LY = 11.5             # Sector arc length along meridian[cite: 1]
CONSTANT MAX_MACRO_INFLUENCERS = 10        # Top gravitational sources summed[cite: 1]

ENUM MassTier:
    TIER_0_SUPERMASSIVE                    # Galactic anchors / SMBH (> 10,000 M_sun)[cite: 1]
    TIER_1_MAJOR                           # Giant stars (10 - 10,000 M_sun)[cite: 1]
    TIER_2_STANDARD                        # Main sequence stars, brown dwarfs (0.001 - 10 M_sun)[cite: 1]
    TIER_3_MINOR                           # Rogue planets, bound planets, asteroids (< 0.001 M_sun)[cite: 1]

STRUCT Entity:
    id: UUID
    sector_id: String                      # Cylindrical sector key (ring, azimuth, col)[cite: 1]
    parent_system_id: UUID                 # NULL if free/rogue; Star ID if bound[cite: 1]
    mass: float                            # kg[cite: 1]
    mu: float                              # G * mass[cite: 1]
    radius: float                          # Collision radius in meters[cite: 1, 3]
    r_hill_max: float                      # Conservative upper bounding ceiling[cite: 1]
    r_soi_max: float                       # Conservative upper bounding ceiling[cite: 1]
    position: Vector3                      # Global Cartesian meters (x, y, z)[cite: 1, 3]
    velocity: Vector3                      # Global Cartesian meters/second (vx, vy, vz)[cite: 1]
    last_indexed_position: Vector3         # Position at last occupant synchronization[cite: 1]
    hill_occupants: List[UUID]             # IDs within r_hill_max[cite: 1]
    is_deleted: Boolean                    # Flag for collision/merger removal[cite: 1]
    is_dirty: Boolean                      # Flag for Pass 2 occupancy refresh[cite: 1]

```

---

### Master Simulation Workflow (Nightly Batch Execution)

```python
FUNCTION RunNightlyGalacticUpdate(delta_t):
    """
    Executes the full nightly batch update across the entire simulated galaxy.
    Hierarchical 3-Pass Galactic Engine coupled with a 2-Phase Star System Engine.
    """
    # ==========================================================================
    # STAGE 1: GALACTIC PASS 1 - MACRO INTEGRATION & TRAJECTORY ADVANCEMENT
    # ==========================================================================
    staged_poses = {}             # Map<EntityID, {position, velocity, sector_id}>
    close_encounter_ledger = []   # List of pairs requiring micro-pass resolution
    dirty_entities = Set()        # Entities requiring Pass 2 metadata/occupant sync
    affected_sectors = Set()      # Sectors requiring mass/density recomputation

    # Query all unparented root entities (Stars, Remnants, and Rogue Bodies)
    root_entities = DB.GetEntitiesWhere("parent_system_id IS NULL AND is_deleted = FALSE")
    tier_0_registry = DB.GetGlobalTier0Supermassives()

    FOR EACH entity IN root_entities:
        # 1. Fetch multi-sector stencil point-masses
        active_sector_keys = ResolveNeighborSectors(entity.position, entity.a_min)
        point_masses = DB.FetchPointMassesForSectors(active_sector_keys)

        # 2. Evaluate net acceleration a(t)
        a_current = ComputeNetAcceleration(entity.position, entity.id, tier_0_registry, point_masses, entity.a_min)

        # 3. Half-step velocity & full-step position update (Velocity Verlet)
        v_half = entity.velocity + (a_current * (0.5 * delta_t))
        p_candidate = entity.position + (v_half * delta_t)

        # 4. Continuous swept-sphere collision & broadphase encounter sweep
        earliest_hit = None
        earliest_hit_time = delta_t
        encounter_partner = None

        FOR EACH other IN point_masses:
            IF other.id == entity.id OR other.is_deleted:
                CONTINUE

            # Gravitational focusing collision cross-section
            v_rel = (entity.velocity - other.velocity).Length()
            v_esc = Sqrt((2.0 * (entity.mu + other.mu)) / (entity.radius + other.radius))
            r_effective = (entity.radius + other.radius) * Sqrt(1.0 + (v_esc * v_esc) / Max(v_rel * v_rel, 1.0))

            hit, t_hit = CheckContinuousSphereOverlap(
                entity.position, v_half,
                other.position, other.velocity,
                r_effective, delta_t
            )

            IF hit AND t_hit < earliest_hit_time:
                earliest_hit = other
                earliest_hit_time = t_hit

            # Check if entering true mutual Hill Sphere
            dist_min = ComputeMinimumSeparation(entity.position, v_half, other.position, other.velocity, delta_t)
            r_encounter = entity.r_hill_max + other.r_hill_max

            IF dist_min <= r_encounter:
                close_encounter_ledger.Append({"entity_a": entity, "entity_b": other})

        # 5. Handle direct physical impacts (Coalescence)
        IF earliest_hit IS NOT NULL:
            merged_record = ExecuteInelasticMerger(entity, earliest_hit, earliest_hit_time)
            staged_poses[entity.id] = merged_record
            dirty_entities.Add(entity.id)
            affected_sectors.Add(entity.sector_id)
            CONTINUE

        # 6. Evaluate a(t + dt) and complete full-step velocity update
        a_new = ComputeNetAcceleration(p_candidate, entity.id, tier_0_registry, point_masses, entity.a_min)
        v_new = v_half + (a_new * (0.5 * delta_t))

        # Check for sector migration
        new_sector_id = CalculateCylindricalSectorKey(p_candidate)
        IF new_sector_id != entity.sector_id:
            affected_sectors.Add(entity.sector_id)
            affected_sectors.Add(new_sector_id)

        # Stage updated state
        staged_poses[entity.id] = {
            "position": p_candidate,
            "velocity": v_new,
            "sector_id": new_sector_id,
            "star_delta_pos": p_candidate - entity.position  # Track displacement vector
        }

    # ==========================================================================
    # STAGE 2: GALACTIC PASS 1.5 - LOCALIZED MICRO-PASS (ENCOUNTERS)
    # ==========================================================================
    # Resolve pairs entering mutual spheres of influence (Bypasses macro step)
    FOR EACH encounter IN close_encounter_ledger:
        a = encounter["entity_a"]
        b = encounter["entity_b"]
        IF a.is_deleted OR b.is_deleted:
            CONTINUE

        ResolveCloseEncounterMicroPass(a, b, delta_t, staged_poses, dirty_entities, affected_sectors)

    # ==========================================================================
    # STAGE 3: STAR SYSTEM ENGINE (PHASES 1 & 2)
    # ==========================================================================
    # Execute planetary and moon system updates coupled to parent stars
    ExecuteStarSystemUpdates(staged_poses, dirty_entities, delta_t)

    # Atomic commit of all staging positions to DB (Everything updates simultaneously)
    DB.BatchApplyStagedPoses(staged_poses)

    # ==========================================================================
    # STAGE 4: GALACTIC PASS 2 - TOPOLOGY & OCCUPANT SYNCHRONIZATION
    # ==========================================================================
    # A. Recompute metadata for affected sectors
    FOR EACH sector_id IN affected_sectors:
        RecomputeSectorPointMassTable(sector_id, trigger_event="SYSTEM_UPDATE")

    # B. Synchronize Max Hill Sphere Occupants for Dirty/Moved Entities
    FOR EACH entity_id IN staged_poses.Keys():
        entity = DB.GetEntity(entity_id)
        IF entity.is_deleted:
            DB.DeleteEntity(entity_id)
            CONTINUE

        disp = (entity.position - entity.last_indexed_position).Length()
        IF disp >= entity.delta_x_min OR entity_id IN dirty_entities:
            SynchronizeHillOccupants(entity)

```

---

### Star System Engine (System-Level 3 Passes)

```python
FUNCTION ExecuteStarSystemUpdates(staged_poses, dirty_entities, delta_t):
    """
    Updates bound planets, moons, and system objects across two phases:
    Phase A: Rigid spatial translation matching star displacement.
    Phase B: Keplerian / perturbed orbit evaluation and Hill boundary crossings.
    """
    all_stars = DB.GetEntitiesWhere("parent_system_id IS NULL AND tier <= MassTier.TIER_2_STANDARD")

    FOR EACH star IN all_stars:
        IF star.id NOT IN staged_poses:
            CONTINUE

        star_update = staged_poses[star.id]
        star_delta_pos = star_update["star_delta_pos"]

        # ----------------------------------------------------------------------
        # PHASE A: STAR DISPLACEMENT PROPAGATION (Pass 1 - Parent Translation)
        # ----------------------------------------------------------------------
        # Fetch all child entities (planets, asteroids, moons) bound to this star
        system_children = DB.GetEntitiesWhere(Format("parent_system_id = '{0}'", star.id))

        FOR EACH child IN system_children:
            # Shift the child's position by the identical displacement of its parent star
            child.position = child.position + star_delta_pos
            child.last_indexed_position = child.last_indexed_position + star_delta_pos

        # ----------------------------------------------------------------------
        # PHASE B: IN-SYSTEM ORBITAL UPDATES & INTERACTION CHECKS (Passes 2 & 3)
        # ----------------------------------------------------------------------
        # Pass 2: Evaluate local two-body/perturbed orbits around parent star
        FOR EACH child IN system_children:
            r_rel = child.position - star_update["position"]
            dist_to_star = r_rel.Length()

            # Dynamic local Hill Radius relative to the host star
            r_hill_local = dist_to_star * ((child.mass / (3.0 * star.mass)) ** (1.0 / 3.0))

            # Advance orbit relative to star using semi-analytical Kepler propagator
            p_child_new, v_child_new = PropagateKeplerianOrbit(
                primary_pos=star_update["position"],
                primary_mu=star.mu,
                child_pos=child.position,
                child_vel=child.velocity,
                dt=delta_t
            )

            # Check if internal child displacement exceeds threshold
            internal_disp = (p_child_new - child.position).Length()
            IF internal_disp >= child.delta_x_min:
                dirty_entities.Add(child.id)

            staged_poses[child.id] = {
                "position": p_child_new,
                "velocity": v_child_new,
                "sector_id": star_update["sector_id"]
            }

            # Pass 3: Planetary Hill Sphere Crossing & System Foreign Interloper Check
            # Check if rogue bodies or other planets enter/exit child's Hill Sphere
            CheckPlanetaryHillCrossings(child, r_hill_local, system_children, staged_poses, dirty_entities)

```

---

### Encounters & Planetary Hill Crossings (Pass 1.5 & Phase B Detail)

```python
FUNCTION CheckPlanetaryHillCrossings(child, r_hill_local, system_children, staged_poses, dirty_entities):
    """
    Detects if an object enters or exits a planet/moon's Hill sphere.
    Handles captures, ejections, Roche disruptions, and surface impacts.
    """
    p_child = staged_poses[child.id]["position"]

    FOR EACH other IN system_children:
        IF other.id == child.id:
            CONTINUE

        p_other = staged_poses[other.id]["position"]
        v_other = staged_poses[other.id]["velocity"]
        dist = (p_child - p_other).Length()

        was_in_hill = other.id IN child.hill_occupants
        is_in_hill = dist <= r_hill_local

        # 1. NEW OBJECT ENTERS PLANETARY HILL SPHERE
        IF is_in_hill AND NOT was_in_hill:
            dirty_entities.Add(child.id)
            dirty_entities.Add(other.id)

            # Compute relative energy to determine capture vs hyperbolic pass
            v_rel = (v_other - staged_poses[child.id]["velocity"]).Length()
            spec_energy = 0.5 * (v_rel ** 2) - (child.mu / dist)

            IF spec_energy < 0.0:
                # Bound capture: Object becomes a captured moon
                DB.UpdateEntityField(other.id, "parent_system_id", child.id)
                child.hill_occupants.Append(other.id)
            ELSE:
                # Hyperbolic slingshot: Trigger localized encounter deflection
                v_defl_child, v_defl_other = ComputeTwoBodyHyperbolicDeflection(child, other)
                staged_poses[child.id]["velocity"] = v_defl_child
                staged_poses[other.id]["velocity"] = v_defl_other

            # Check Fluid Roche Limit Disruption
            d_roche = 2.44 * child.radius * ((child.density / other.density) ** (1.0 / 3.0))
            IF dist <= d_roche AND child.mass > other.mass:
                # Tidal fragmentation: Destroy secondary entity and stage ring formation
                DB.UpdateEntityField(other.id, "is_deleted", TRUE)
                child.mass += other.mass
                child.mu = G * child.mass

        # 2. OBJECT EXITS PLANETARY HILL SPHERE (Ejection)
        ELSE IF NOT is_in_hill AND was_in_hill:
            dirty_entities.Add(child.id)
            dirty_entities.Add(other.id)

            # Sever moon hierarchy, promote object back to star-level child
            DB.UpdateEntityField(other.id, "parent_system_id", child.parent_system_id)
            child.hill_occupants.Remove(other.id)


FUNCTION SynchronizeHillOccupants(entity):
    """
    Galactic Pass 2: Refreshes the entity's hill_occupants list in SQL
    and updates reciprocal entries in neighboring bodies.
    """
    active_sector_keys = ResolveNeighborSectors(entity.position, entity.a_min)
    candidate_bodies = DB.FetchPointMassesForSectors(active_sector_keys)

    new_occupants = []

    FOR EACH other IN candidate_bodies:
        IF other.id == entity.id:
            CONTINUE

        dist = (entity.position - other.position).Length()

        # Update entity's occupant list
        IF dist <= entity.r_hill_max:
            new_occupants.Append(other.id)

        # Reciprocal asymmetric update: check if entity is inside other's max Hill sphere
        IF dist <= other.r_hill_max:
            IF entity.id NOT IN other.hill_occupants:
                DB.AppendHillOccupant(primary_id=other.id, occupant_id=entity.id)
        ELSE:
            IF entity.id IN other.hill_occupants:
                DB.RemoveHillOccupant(primary_id=other.id, occupant_id=entity.id)

    # Persist updated occupants and reset indexing baseline
    DB.UpdateEntity(entity.id, {
        "hill_occupants": new_occupants,
        "last_indexed_position": entity.position,
        "is_dirty": FALSE
    })

```

---

### Helper Routines: Topology, Acceleration, & Collisions

```python
FUNCTION ResolveNeighborSectors(eval_pos, a_min):
    """
    Cylindrical arc sector stencil mapping using multiples of 6.
    Resolves 25-33 intersecting sectors across non-aligned radial shells.
    """
    r = Sqrt(eval_pos.x ** 2 + eval_pos.y ** 2) / LIGHT_YEAR_METERS
    theta = Modulo(atan2(eval_pos.y, eval_pos.x), 2.0 * PI)
    z = eval_pos.z / LIGHT_YEAR_METERS

    ring_idx = Floor(r / DELTA_R_LY)
    col_idx = Floor(z / DELTA_Z_LY)

    # Radial shell sector quantization (Multiples of 6)
    r_mid = (ring_idx + 0.5) * DELTA_R_LY
    n_theta = Max(6, 6 * Round((2.0 * PI * r_mid) / (6.0 * TARGET_ARC_LY)))
    d_theta = (2.0 * PI) / n_theta

    # Search window buffer
    sector_keys = []
    FOR r_i IN [Max(0, ring_idx - 1), ring_idx, ring_idx + 1]:
        n_divs = Max(6, 6 * Round((2.0 * PI * (r_i + 0.5) * DELTA_R_LY) / (6.0 * TARGET_ARC_LY)))
        d_th = (2.0 * PI) / n_divs

        j_center = Floor(theta / d_th)
        FOR j_offset IN [-1, 0, 1]:
            j_wrapped = Modulo(j_center + j_offset, n_divs)
            FOR c_k IN [col_idx - 1, col_idx, col_idx + 1]:
                sector_keys.Append(Format("{0}:{1}:{2}", r_i, j_wrapped, c_k))

    RETURN Unique(sector_keys)


FUNCTION CheckContinuousSphereOverlap(p1, v1, p2, v2, r_combined, dt):
    """
    Continuous Collision Detection (CCD): Ray vs Capsule intersection over [0, dt].
    Prevents tunneling during large integration steps.
    """
    dp = p1 - p2
    dv = v1 - v2

    A = dv.Dot(dv)
    B = 2.0 * dp.Dot(dv)
    C = dp.Dot(dp) - (r_combined * r_combined)

    IF C <= 0.0:
        RETURN TRUE, 0.0  # Already overlapping

    IF A <= 1e-15:
        RETURN FALSE, dt  # Trajectories parallel

    discriminant = (B * B) - (4.0 * A * C)
    IF discriminant < 0.0:
        RETURN FALSE, dt  # No intersection

    t_entry = (-B - Sqrt(discriminant)) / (2.0 * A)
    IF 0.0 <= t_entry <= dt:
        RETURN TRUE, t_entry

    RETURN FALSE, dt


FUNCTION ExecuteInelasticMerger(primary, secondary, t_hit):
    """
    Conserves linear momentum during an inelastic physical impact.
    Consumes secondary body into primary body.
    """
    p_impact_prim = primary.position + (primary.velocity * t_hit)
    p_impact_sec = secondary.position + (secondary.velocity * t_hit)

    total_mass = primary.mass + secondary.mass
    v_merged = (primary.velocity * primary.mass + secondary.velocity * secondary.mass) / total_mass
    r_merged = ((primary.radius ** 3) + (secondary.radius ** 3)) ** (1.0 / 3.0)

    secondary.is_deleted = TRUE
    DB.UpdateEntityField(secondary.id, "is_deleted", TRUE)

    RETURN {
        "position": p_impact_prim,
        "velocity": v_merged,
        "mass": total_mass,
        "mu": G * total_mass,
        "radius": r_merged,
        "sector_id": primary.sector_id,
        "star_delta_pos": p_impact_prim - primary.position
    }

```