"""
config.py
=========
Metaclass-driven configuration system for inelastic-rebound-md simulations.

Design highlights:
- ``FieldDescriptor``: a data descriptor that owns type coercion and validation.
- ``ConfigMeta``: metaclass that collects all FieldDescriptors declared on a
  class body, builds ``__slots__``, and wires validation into ``__init__``.
- ``SimulationConfig`` / ``EquilibrationConfig`` / ``ImpactConfig``:
  concrete configuration dataclasses built with the metaclass.
- ``ConfigRegistry``: maps string keys to config classes (metaprogramming
  application of ``__init_subclass__``).
- JSON/TOML round-trip serialisation built in.

Usage
-----
>>> cfg = ImpactConfig(v0=2.5, angle=30.0)
>>> cfg.v0
2.5
>>> cfg.to_dict()
{'v0': 2.5, 'angle': 30.0, ...}
>>> cfg2 = ImpactConfig.from_dict({'v0': 3.0})
"""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any, Callable, ClassVar, Optional, Type, TypeVar

from .exceptions import MissingParameterError, ValidationError

T = TypeVar("T")


# ==============================================================================
# Field descriptor
# ==============================================================================

class FieldDescriptor:
    """Data descriptor for a single configuration parameter.

    Owns type coercion, bounds checking, and default values.
    Instances are created at class-definition time (not instance time).

    Parameters
    ----------
    dtype : type
        Python type to coerce the value to.
    default : Any
        Sentinel ``REQUIRED`` if the field has no default, else the default value.
    validator : Callable[[Any], None] | None
        Extra validation callable.  Receives the coerced value; should raise
        ``ValidationError`` on failure.
    doc : str
        Human-readable description of the field.
    unit : str
        Physical unit string, for documentation and serialisation.
    """

    REQUIRED: ClassVar[object] = object()  # sentinel

    def __init__(
        self,
        dtype: type,
        default: Any = REQUIRED,
        *,
        validator: Optional[Callable[[Any], None]] = None,
        doc: str = "",
        unit: str = "",
    ) -> None:
        self.dtype = dtype
        self.default = default
        self.validator = validator
        self.doc = doc
        self.unit = unit
        self.name: str = ""          # set by ConfigMeta.__new__
        self.private: str = ""       # mangled storage name

    def __set_name__(self, owner: type, name: str) -> None:
        self.name = name
        self.private = f"_field_{name}"

    def __get__(self, obj: Any, objtype: type = None) -> Any:
        if obj is None:
            return self
        try:
            return getattr(obj, self.private)
        except AttributeError:
            if self.default is FieldDescriptor.REQUIRED:
                raise MissingParameterError(self.name)
            return self.default

    def __set__(self, obj: Any, value: Any) -> None:
        if value is None and self.default is not FieldDescriptor.REQUIRED:
            setattr(obj, self.private, self.default)
            return
        try:
            coerced = self.dtype(value)
        except (TypeError, ValueError) as exc:
            raise ValidationError(
                self.name, value, f"convertible to {self.dtype.__name__}"
            ) from exc
        if self.validator is not None:
            self.validator(coerced)
        setattr(obj, self.private, coerced)

    def __repr__(self) -> str:
        return (
            f"FieldDescriptor({self.name!r}, dtype={self.dtype.__name__}, "
            f"default={self.default!r}, unit={self.unit!r})"
        )


# ==============================================================================
# Convenience validator factories
# ==============================================================================

def positive(name: str = "") -> Callable[[float], None]:
    """Return a validator that requires value > 0."""
    def _validate(v: float) -> None:
        if v <= 0:
            raise ValidationError(name or "value", v, "positive (> 0)")
    return _validate


def in_range(lo: float, hi: float, name: str = "") -> Callable[[float], None]:
    """Return a validator that requires lo <= value <= hi."""
    def _validate(v: float) -> None:
        if not (lo <= v <= hi):
            raise ValidationError(name or "value", v, f"in [{lo}, {hi}]")
    return _validate


def non_negative(name: str = "") -> Callable[[float], None]:
    """Return a validator that requires value >= 0."""
    def _validate(v: float) -> None:
        if v < 0:
            raise ValidationError(name or "value", v, "non-negative (>= 0)")
    return _validate


# ==============================================================================
# Metaclass
# ==============================================================================

class ConfigMeta(type):
    """Metaclass for simulation configuration classes.

    Responsibilities:
    1. Collects all ``FieldDescriptor`` instances declared in the class body.
    2. Stores them in ``cls._fields`` (ordered dict: name -> descriptor).
    3. Generates ``__init__`` that accepts keyword arguments for all fields.
    4. Generates ``__repr__``, ``__eq__``, ``to_dict``, ``from_dict``.
    5. Registers concrete subclasses in ``ConfigRegistry``.
    """

    def __new__(
        mcs,
        name: str,
        bases: tuple[type, ...],
        namespace: dict[str, Any],
        **kwargs: Any,
    ) -> "ConfigMeta":
        fields: dict[str, FieldDescriptor] = {}

        # Inherit fields from base classes
        for base in reversed(bases):
            if hasattr(base, "_fields"):
                fields.update(base._fields)

        # Collect new fields declared in this class body
        for attr, val in namespace.items():
            if isinstance(val, FieldDescriptor):
                val.__set_name__(None, attr)   # ensure name is set
                fields[attr] = val

        namespace["_fields"] = fields
        cls = super().__new__(mcs, name, bases, namespace)

        # Auto-generate __init__
        def __init__(self: Any, **kwargs: Any) -> None:
            for field_name, descriptor in self._fields.items():
                if field_name in kwargs:
                    setattr(self, field_name, kwargs.pop(field_name))
                elif descriptor.default is FieldDescriptor.REQUIRED:
                    raise MissingParameterError(field_name)
                # Else: descriptor.__get__ will return default lazily
            if kwargs:
                unknown = ", ".join(f"'{k}'" for k in kwargs)
                raise TypeError(f"{name}.__init__() got unexpected keyword argument(s): {unknown}")

        cls.__init__ = __init__

        # Auto-generate __repr__
        def __repr__(self: Any) -> str:
            parts = [
                f"{k}={getattr(self, k)!r}"
                for k in self._fields
            ]
            return f"{type(self).__name__}({', '.join(parts)})"

        cls.__repr__ = __repr__

        # Auto-generate __eq__
        def __eq__(self: Any, other: Any) -> bool:
            if type(self) is not type(other):
                return NotImplemented
            return all(
                getattr(self, k) == getattr(other, k)
                for k in self._fields
            )

        cls.__eq__ = __eq__

        # to_dict
        def to_dict(self: Any) -> dict[str, Any]:
            return {k: getattr(self, k) for k in self._fields}

        cls.to_dict = to_dict

        # from_dict (classmethod)
        @classmethod
        def from_dict(cls_inner: type, data: dict[str, Any]) -> Any:
            return cls_inner(**data)

        cls.from_dict = from_dict

        # to_json / from_json
        def to_json(self: Any, path: Optional[Path | str] = None) -> str:
            payload = {"_type": type(self).__name__, **self.to_dict()}
            text = json.dumps(payload, indent=2)
            if path is not None:
                Path(path).write_text(text)
            return text

        cls.to_json = to_json

        @classmethod
        def from_json(cls_inner: type, path: Path | str) -> Any:
            data = json.loads(Path(path).read_text())
            data.pop("_type", None)
            return cls_inner.from_dict(data)

        cls.from_json = from_json

        return cls


# ==============================================================================
# Registry via __init_subclass__
# ==============================================================================

class ConfigRegistry:
    """Global registry mapping config class names to their types.

    Populated automatically via ``BaseConfig.__init_subclass__``.
    """
    _registry: ClassVar[dict[str, Type["BaseConfig"]]] = {}

    @classmethod
    def register(cls, name: str, config_class: type) -> None:
        cls._registry[name] = config_class

    @classmethod
    def get(cls, name: str) -> Type["BaseConfig"]:
        if name not in cls._registry:
            raise KeyError(f"No config class registered under '{name}'. "
                           f"Available: {list(cls._registry)}")
        return cls._registry[name]

    @classmethod
    def from_json_auto(cls, path: Path | str) -> "BaseConfig":
        """Load a config JSON and auto-detect its type from the '_type' key."""
        data = json.loads(Path(path).read_text())
        type_name = data.pop("_type", None)
        if type_name is None:
            raise ValueError(f"Config file '{path}' missing '_type' key.")
        config_cls = cls.get(type_name)
        return config_cls.from_dict(data)


class BaseConfig(metaclass=ConfigMeta):
    """Abstract base for all configuration classes.

    All concrete subclasses are auto-registered in ``ConfigRegistry``.
    """

    def __init_subclass__(cls, **kwargs: Any) -> None:
        super().__init_subclass__(**kwargs)
        ConfigRegistry.register(cls.__name__, cls)


# ==============================================================================
# Concrete configuration classes
# ==============================================================================

class EquilibrationConfig(BaseConfig):
    """Configuration for the surface equilibration phase (in.equilibrate).

    All values in LJ reduced units unless noted.
    """

    T_target   = FieldDescriptor(float, 0.1,    validator=positive("T_target"),
                                  doc="Target temperature", unit="eps/kB")
    T_damp     = FieldDescriptor(float, 1.0,    validator=positive("T_damp"),
                                  doc="Nose-Hoover damping time", unit="tau")
    a0         = FieldDescriptor(float, 1.5496, validator=positive("a0"),
                                  doc="FCC lattice constant", unit="sigma")
    Nx         = FieldDescriptor(int,   8,      validator=positive("Nx"),
                                  doc="Unit cells in x")
    Ny         = FieldDescriptor(int,   8,      validator=positive("Ny"),
                                  doc="Unit cells in y")
    Nz         = FieldDescriptor(int,   6,      validator=positive("Nz"),
                                  doc="Unit cells in z (slab thickness)")
    steps_equil = FieldDescriptor(int, 50000,   validator=positive("steps_equil"),
                                   doc="NVT equilibration steps")
    steps_prod  = FieldDescriptor(int, 30000,   validator=positive("steps_prod"),
                                   doc="Production run steps")
    dt         = FieldDescriptor(float, 0.005,  validator=positive("dt"),
                                  doc="Timestep", unit="tau")
    seed       = FieldDescriptor(int,   12345,
                                  doc="Random seed for velocity initialisation")

    @property
    def total_atoms(self) -> int:
        """Number of atoms in the FCC slab (4 atoms/unit cell)."""
        return 4 * self.Nx * self.Ny * self.Nz

    @property
    def slab_height(self) -> float:
        """Approximate height of the slab in sigma."""
        return self.Nz * self.a0


class ImpactConfig(BaseConfig):
    """Configuration for the projectile impact phase (in.impact).

    All values in LJ reduced units unless noted.
    """

    v0       = FieldDescriptor(float, 2.0,  validator=positive("v0"),
                                doc="Incident speed", unit="sigma/tau")
    angle    = FieldDescriptor(float, 0.0,  validator=in_range(0.0, 89.9, "angle"),
                                doc="Impact angle from surface normal", unit="degrees")
    T_target = FieldDescriptor(float, 0.1,  validator=positive("T_target"),
                                doc="Surface temperature", unit="eps/kB")
    T_damp   = FieldDescriptor(float, 1.0,  validator=positive("T_damp"),
                                doc="Nose-Hoover damping time", unit="tau")
    dt       = FieldDescriptor(float, 0.002, validator=positive("dt"),
                                doc="Timestep (smaller than equilibration)", unit="tau")
    z_gap    = FieldDescriptor(float, 8.0,  validator=positive("z_gap"),
                                doc="Vacuum gap above surface", unit="sigma")
    steps_approach = FieldDescriptor(int, 10000, validator=positive("steps_approach"),
                                      doc="Steps in approach phase")
    steps_impact   = FieldDescriptor(int, 80000, validator=positive("steps_impact"),
                                      doc="Steps in impact + rebound phase")
    thermo_freq    = FieldDescriptor(int, 50,    validator=positive("thermo_freq"),
                                      doc="Thermo output frequency", unit="steps")
    dump_freq      = FieldDescriptor(int, 100,   validator=positive("dump_freq"),
                                      doc="Trajectory dump frequency", unit="steps")

    @property
    def vz0(self) -> float:
        """Normal (z) component of incident velocity (negative = downward)."""
        return -self.v0 * math.cos(math.radians(self.angle))

    @property
    def vx0(self) -> float:
        """Tangential (x) component of incident velocity."""
        return self.v0 * math.sin(math.radians(self.angle))

    @property
    def total_time(self) -> float:
        """Total simulated time [tau]."""
        return (self.steps_approach + self.steps_impact) * self.dt


class ParametricSweepConfig(BaseConfig):
    """Configuration for a parametric study over v0 and/or angle.

    Used by the runner to launch multiple ImpactConfig instances.
    """

    v0_min   = FieldDescriptor(float, 0.5,  validator=positive("v0_min"),  unit="sigma/tau")
    v0_max   = FieldDescriptor(float, 5.0,  validator=positive("v0_max"),  unit="sigma/tau")
    v0_steps = FieldDescriptor(int,   10,   validator=positive("v0_steps"))
    angle_list = FieldDescriptor(str, "0.0",
                                  doc="Comma-separated list of impact angles [degrees]")
    T_surface  = FieldDescriptor(float, 0.1, validator=positive("T_surface"), unit="eps/kB")
    T_surface_list = FieldDescriptor(str, "",
                                      doc="Comma-separated list of surface temperatures [eps/kB]. "
                                          "If non-empty, overrides T_surface; runs one equilibration "
                                          "per temperature.")
    n_replicas = FieldDescriptor(int,   3,   validator=positive("n_replicas"),
                                  doc="Number of independent replicas per condition")
    parallel_jobs = FieldDescriptor(int, 4,  validator=positive("parallel_jobs"),
                                     doc="Max concurrent LAMMPS processes")

    @property
    def v0_values(self) -> list[float]:
        """Linearly spaced v0 values for the sweep."""
        if self.v0_steps == 1:
            return [self.v0_min]
        step = (self.v0_max - self.v0_min) / (self.v0_steps - 1)
        return [self.v0_min + i * step for i in range(self.v0_steps)]

    @property
    def angles(self) -> list[float]:
        """Parsed list of impact angles."""
        return [float(a.strip()) for a in self.angle_list.split(",")]

    @property
    def temperatures(self) -> list[float]:
        """Parsed list of surface temperatures.

        If ``T_surface_list`` is non-empty, returns that list.
        Otherwise returns ``[T_surface]`` (single temperature, backward-compatible).
        """
        s = self.T_surface_list.strip()
        if s:
            return [float(t.strip()) for t in s.split(",")]
        return [self.T_surface]
