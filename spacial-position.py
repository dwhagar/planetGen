import math
from typing import Dict, Optional, Tuple, Any

class SpatialPosition3D:
    """
    3D hierarchical position and kinematic tracker maintaining synchronized
    Galactic, Sector, and Star/System frames in Cartesian, Cylindrical,
    and Spherical coordinates.
    """

    # Minimum displacement thresholds required to register observable movement
    MIN_OBSERVABLE_GALACTIC: float = 1e-4  # e.g., in parsecs or light-years
    MIN_OBSERVABLE_SECTOR: float = 1e-5
    MIN_OBSERVABLE_SYSTEM: float = 1e-6    # e.g., in AU or km

    def __init__(
        self,
        galactic_cartesian: Tuple[float, float, float],
        sector_center_galactic: Tuple[float, float, float],
        velocity_vector_cartesian: Tuple[float, float, float] = (0.0, 0.0, 0.0),
        star_center_galactic: Optional[Tuple[float, float, float]] = None,
        is_star: bool = False,
        point_mass_data: Optional[Dict[str, Any]] = None,
    ):
        # Frame anchor origins in galactic Cartesian coordinates
        self._sector_center_gal = tuple(sector_center_galactic)
        self._star_center_gal = tuple(star_center_galactic) if star_center_galactic else None
        self._is_star = is_star

        # Internal state stores
        self._coords: Dict[str, Dict[str, Tuple[float, float, float]]] = {
            "galactic": {},
            "sector": {},
            "system": {},
        }
        self._velocity: Dict[str, Any] = {}
        self._time_to_observable_movement: Dict[str, float] = {}
        self._point_mass_table: Dict[str, Any] = point_mass_data if point_mass_data is not None else {}

        # Initial full sync
        self._set_state(
            galactic_cartesian=galactic_cartesian,
            velocity_vector=velocity_vector_cartesian,
        )

    # -------------------------------------------------------------------------
    # Coordinate Conversion Helpers (Private)
    # -------------------------------------------------------------------------

    @staticmethod
    def _cartesian_to_cylindrical(x: float, y: float, z: float) -> Tuple[float, float, float]:
        r = math.sqrt(x**2 + y**2)
        theta = math.atan2(y, x)
        return (r, theta, z)

    @staticmethod
    def _cartesian_to_spherical(x: float, y: float, z: float) -> Tuple[float, float, float]:
        r = math.sqrt(x**2 + y**2 + z**2)
        if r == 0.0:
            return (0.0, 0.0, 0.0)
        theta = math.atan2(y, x)
        phi = math.acos(max(min(z / r, 1.0), -1.0))
        return (r, theta, phi)

    @staticmethod
    def _cylindrical_to_cartesian(r: float, theta: float, z: float) -> Tuple[float, float, float]:
        x = r * math.cos(theta)
        y = r * math.sin(theta)
        return (x, y, z)

    @staticmethod
    def _spherical_to_cartesian(r: float, theta: float, phi: float) -> Tuple[float, float, float]:
        x = r * math.sin(phi) * math.cos(theta)
        y = r * math.sin(phi) * math.sin(theta)
        z = r * math.cos(phi)
        return (x, y, z)

    # -------------------------------------------------------------------------
    # Comprehensive Recalculation Pipeline
    # -------------------------------------------------------------------------

    def _sync_frame(self, frame_name: str, x: float, y: float, z: float) -> None:
        self._coords[frame_name] = {
            "cartesian": (x, y, z),
            "cylindrical": self._cartesian_to_cylindrical(x, y, z),
            "spherical": self._cartesian_to_spherical(x, y, z),
        }

    def _update_kinematics(self, vx: float, vy: float, vz: float) -> None:
        speed = math.sqrt(vx**2 + vy**2 + vz**2)
        if speed > 0.0:
            unit_dir = (vx / speed, vy / speed, vz / speed)
        else:
            unit_dir = (0.0, 0.0, 0.0)

        self._velocity = {
            "cartesian": (vx, vy, vz),
            "speed": speed,
            "direction_unit": unit_dir,
        }

        # Observable movement recalculation thresholds
        if speed > 0.0:
            self._time_to_observable_movement = {
                "galactic": self.MIN_OBSERVABLE_GALACTIC / speed,
                "sector": self.MIN_OBSERVABLE_SECTOR / speed,
                "system": self.MIN_OBSERVABLE_SYSTEM / speed,
                "minimum_next_recalc": min(
                    self.MIN_OBSERVABLE_GALACTIC,
                    self.MIN_OBSERVABLE_SECTOR,
                    self.MIN_OBSERVABLE_SYSTEM,
                ) / speed,
            }
        else:
            inf = float("inf")
            self._time_to_observable_movement = {
                "galactic": inf,
                "sector": inf,
                "system": inf,
                "minimum_next_recalc": inf,
            }

    def _set_state(
        self,
        galactic_cartesian: Tuple[float, float, float],
        velocity_vector: Optional[Tuple[float, float, float]] = None,
    ) -> None:
        gx, gy, gz = galactic_cartesian

        # 1. Galactic Frame
        self._sync_frame("galactic", gx, gy, gz)

        # 2. Sector Frame (Relative to Sector Origin)
        sx = gx - self._sector_center_gal[0]
        sy = gy - self._sector_center_gal[1]
        sz = gz - self._sector_center_gal[2]
        self._sync_frame("sector", sx, sy, sz)

        # 3. System / Nearest Star Frame
        if self._is_star or self._star_center_gal is None:
            # Stars have no system coordinates; objects with no star reference omit it
            self._coords["system"] = {}
        else:
            tx = gx - self._star_center_gal[0]
            ty = gy - self._star_center_gal[1]
            tz = gz - self._star_center_gal[2]
            self._sync_frame("system", tx, ty, tz)

        # 4. Kinematics & Recalculation Thresholds
        if velocity_vector is not None:
            self._update_kinematics(*velocity_vector)
        else:
            curr_v = self._velocity.get("cartesian", (0.0, 0.0, 0.0))
            self._update_kinematics(*curr_v)

    # -------------------------------------------------------------------------
    # Frame Mutators (Triggers Full Synchronization)
    # -------------------------------------------------------------------------

    def set_galactic_cartesian(self, x: float, y: float, z: float) -> None:
        self._set_state(galactic_cartesian=(x, y, z))

    def set_galactic_cylindrical(self, r: float, theta: float, z: float) -> None:
        cartesian = self._cylindrical_to_cartesian(r, theta, z)
        self._set_state(galactic_cartesian=cartesian)

    def set_galactic_spherical(self, r: float, theta: float, phi: float) -> None:
        cartesian = self._spherical_to_cartesian(r, theta, phi)
        self._set_state(galactic_cartesian=cartesian)

    def set_sector_cartesian(self, x: float, y: float, z: float) -> None:
        gx = x + self._sector_center_gal[0]
        gy = y + self._sector_center_gal[1]
        gz = z + self._sector_center_gal[2]
        self._set_state(galactic_cartesian=(gx, gy, gz))

    def set_sector_cylindrical(self, r: float, theta: float, z: float) -> None:
        cx, cy, cz = self._cylindrical_to_cartesian(r, theta, z)
        self.set_sector_cartesian(cx, cy, cz)

    def set_sector_spherical(self, r: float, theta: float, phi: float) -> None:
        cx, cy, cz = self._spherical_to_cartesian(r, theta, phi)
        self.set_sector_cartesian(cx, cy, cz)

    def set_system_cartesian(self, x: float, y: float, z: float) -> None:
        if self._is_star or self._star_center_gal is None:
            raise ValueError("Cannot set system coordinates for a star or unbound object without a star reference.")
        gx = x + self._star_center_gal[0]
        gy = y + self._star_center_gal[1]
        gz = z + self._star_center_gal[2]
        self._set_state(galactic_cartesian=(gx, gy, gz))

    def set_system_cylindrical(self, r: float, theta: float, z: float) -> None:
        cx, cy, cz = self._cylindrical_to_cartesian(r, theta, z)
        self.set_system_cartesian(cx, cy, cz)

    def set_system_spherical(self, r: float, theta: float, phi: float) -> None:
        cx, cy, cz = self._spherical_to_cartesian(r, theta, phi)
        self.set_system_cartesian(cx, cy, cz)

    def set_velocity_cartesian(self, vx: float, vy: float, vz: float) -> None:
        curr_pos = self._coords["galactic"]["cartesian"]
        self._set_state(galactic_cartesian=curr_pos, velocity_vector=(vx, vy, vz))

    def set_nearest_star_center(self, star_galactic_cartesian: Optional[Tuple[float, float, float]]) -> None:
        self._star_center_gal = tuple(star_galactic_cartesian) if star_galactic_cartesian else None
        curr_pos = self._coords["galactic"]["cartesian"]
        self._set_state(galactic_cartesian=curr_pos)

    def set_sector_center(self, sector_galactic_cartesian: Tuple[float, float, float]) -> None:
        self._sector_center_gal = tuple(sector_galactic_cartesian)
        curr_pos = self._coords["galactic"]["cartesian"]
        self._set_state(galactic_cartesian=curr_pos)

    # -------------------------------------------------------------------------
    # Getters (Read-Only Access via Functions)
    # -------------------------------------------------------------------------

    def get_coordinates(self, frame: str, coord_type: str) -> Optional[Tuple[float, float, float]]:
        frame_key = frame.lower()
        type_key = coord_type.lower()
        if frame_key not in self._coords:
            raise KeyError(f"Invalid frame '{frame}'. Valid frames: 'galactic', 'sector', 'system'.")
        if type_key not in ("cartesian", "cylindrical", "spherical"):
            raise KeyError(f"Invalid coordinate type '{coord_type}'. Valid types: 'cartesian', 'cylindrical', 'spherical'.")
        
        return self._coords[frame_key].get(type_key)

    def get_velocity_vector(self) -> Tuple[float, float, float]:
        return self._velocity["cartesian"]

    def get_velocity_direction(self) -> Tuple[float, float, float]:
        return self._velocity["direction_unit"]

    def get_speed(self) -> float:
        return self._velocity["speed"]

    def get_time_to_observable_movement(self, frame: Optional[str] = None) -> float:
        if frame is None:
            return self._time_to_observable_movement["minimum_next_recalc"]
        frame_key = frame.lower()
        if frame_key not in self._time_to_observable_movement:
            raise KeyError(f"Invalid frame '{frame}'. Valid frames: 'galactic', 'sector', 'system'.")
        return self._time_to_observable_movement[frame_key]

    def get_point_mass_data(self) -> Dict[str, Any]:
        return dict(self._point_mass_table)

    def update_point_mass_data(self, key: str, value: Any) -> None:
        self._point_mass_table[key] = value

    def is_star(self) -> bool:
        return self._is_star