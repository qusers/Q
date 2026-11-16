"""Find, layer and resolve QligFEP configuration files.

Layers, lowest to highest priority:

1. the packaged ``defaults.yaml`` (historical cluster profiles and the default protocol);
2. the user file, ``~/.config/qligfep/config.yaml`` or ``$QLIGFEP_CONFIG`` (profiles only);
3. the project file, ``./qligfep.yaml`` or an explicit path (profiles and protocols).

Command-line flags are applied on top by the CLIs. A name defined in several layers is
merged key by key, so a higher layer only states what it changes.
"""

from __future__ import annotations

import difflib
import os
from pathlib import Path

import yaml

from .schema import (
    PROFILE_BASE,
    PROFILE_KEYS,
    RESOURCE_KEYS,
    SLURM_KEYS,
    SYSTEMS,
    TOP_KEYS,
    ConfigError,
    Profile,
    Protocol,
    build_profile,
    build_protocol,
    check_keys,
    deep_merge,
    validate_protocol_definition,
)

CONFIG_VERSION = 1
PACKAGED_DEFAULTS = Path(__file__).with_name("defaults.yaml")
PROJECT_FILE_NAMES = ("qligfep.yaml", "qligfep.yml")
PACKAGED_LABEL = "packaged defaults"


def user_config_path() -> Path:
    """Location of the per-user configuration file (it may not exist)."""
    explicit = os.environ.get("QLIGFEP_CONFIG")
    if explicit:
        return Path(explicit).expanduser()
    base = os.environ.get("XDG_CONFIG_HOME") or "~/.config"
    return Path(base).expanduser() / "qligfep" / "config.yaml"


def find_project_config(cwd: str | os.PathLike | None = None) -> Path | None:
    """Return the project file in ``cwd`` if there is one."""
    directory = Path(cwd) if cwd else Path.cwd()
    for name in PROJECT_FILE_NAMES:
        candidate = directory / name
        if candidate.is_file():
            return candidate
    return None


def _read_yaml(path: Path) -> dict:
    try:
        data = yaml.safe_load(path.read_text())
    except yaml.YAMLError as error:
        raise ConfigError(f"{path}: not valid YAML.\n{error}") from error
    except OSError as error:
        raise ConfigError(f"{path}: cannot be read ({error.strerror}).") from error
    if data is None:
        return {}
    if not isinstance(data, dict):
        raise ConfigError(f"{path}: the top level must be a mapping.")
    return data


def _validate_layer(data: dict, label: str, allow_protocols: bool) -> None:
    check_keys(data, TOP_KEYS, label)
    version = data.get("version", CONFIG_VERSION)
    if version != CONFIG_VERSION:
        raise ConfigError(
            f"{label}: config version {version!r} is not supported by this QligFEP "
            f"(expects {CONFIG_VERSION})."
        )
    profiles = data.get("profiles") or {}
    if not isinstance(profiles, dict):
        raise ConfigError(f"{label}: 'profiles' must be a mapping of name -> profile.")
    for name, raw in profiles.items():
        _validate_profile_definition(raw or {}, f"{label}: profiles.{name}")
    protocols = data.get("protocols") or {}
    if not isinstance(protocols, dict):
        raise ConfigError(f"{label}: 'protocols' must be a mapping of name -> protocol.")
    if protocols and not allow_protocols:
        raise ConfigError(
            f"{label}: protocols ({', '.join(protocols)}) are not read from the user file. "
            "They define the science of a project, so keep them in the project's qligfep.yaml "
            "where they travel with the results."
        )
    for name, raw in protocols.items():
        validate_protocol_definition(name, raw or {}, label)


def _validate_profile_definition(raw: dict, where: str, nested: bool = False) -> None:
    allowed = tuple(k for k in PROFILE_KEYS if k not in ("extends", "systems")) if nested else PROFILE_KEYS
    check_keys(raw, allowed, where)
    if "resources" in raw:
        check_keys(raw["resources"], RESOURCE_KEYS, f"{where}.resources")
    if "slurm" in raw:
        check_keys(raw["slurm"], SLURM_KEYS, f"{where}.slurm")
    for system, override in (raw.get("systems") or {}).items():
        if system not in SYSTEMS:
            raise ConfigError(f"{where}.systems: '{system}' is not one of {', '.join(SYSTEMS)}.")
        _validate_profile_definition(override or {}, f"{where}.systems.{system}", nested=True)


class Config:
    """The layered configuration, resolved on demand into profiles and protocols."""

    def __init__(self, layers: list[tuple[str, dict]]):
        self.layers = layers

    # -- lookup helpers -------------------------------------------------------------

    @property
    def sources(self) -> list[str]:
        return [label for label, _ in self.layers]

    def _names(self, kind: str) -> list[str]:
        names: dict[str, None] = {}
        for _, data in self.layers:
            names.update(dict.fromkeys(data.get(kind) or {}))
        return list(names)

    def profile_names(self) -> list[str]:
        return self._names("profiles")

    def protocol_names(self) -> list[str]:
        return self._names("protocols")

    def _setting(self, key: str):
        value = None
        for _, data in self.layers:
            if data.get(key) is not None:
                value = data[key]
        return value

    @property
    def default_profile(self) -> str | None:
        return self._setting("default_profile")

    @property
    def default_protocol(self) -> str:
        return self._setting("default_protocol") or "default"

    def _definition(self, kind: str, name: str, skip_packaged: bool = False) -> tuple[dict, tuple[str, ...]]:
        """Merge every layer's entry for ``name`` (without following ``extends``)."""
        merged: dict = {}
        sources = []
        for label, data in self.layers:
            entry = (data.get(kind) or {}).get(name)
            if entry is None and name not in (data.get(kind) or {}):
                continue
            sources.append(label)
            if skip_packaged and label == PACKAGED_LABEL:
                continue
            merged = deep_merge(merged, entry or {})
        return merged, tuple(sources)

    def _unknown(self, kind: str, name: str) -> ConfigError:
        names = self._names(kind)
        close = difflib.get_close_matches(str(name), names, n=1)
        hint = f" Did you mean '{close[0]}'?" if close else ""
        singular = kind[:-1]
        return ConfigError(f"Unknown {singular} '{name}'.{hint} Available {kind}: {', '.join(names) or 'none'}.")

    def _inherited(self, kind: str, name: str, chain: tuple[str, ...] = ()) -> tuple[dict, tuple[str, ...]]:
        """Resolve ``name`` through its ``extends`` chain into one mapping."""
        if name in chain:
            raise ConfigError(f"{kind}: circular 'extends': {' -> '.join(chain + (name,))}.")
        if name not in self._names(kind):
            if chain:
                raise ConfigError(
                    f"{kind}.{chain[-1]}.extends: {self._unknown(kind, name)}"
                )
            raise self._unknown(kind, name)
        own, sources = self._definition(kind, name)
        own = dict(own)
        parent = own.pop("extends", None)
        if parent is None and kind == "protocols" and name != "default":
            parent = "default"
        if parent is None:
            base, base_sources = (PROFILE_BASE if kind == "profiles" else {}), ()
        else:
            base, base_sources = self._inherited(kind, str(parent), chain + (name,))
        return deep_merge(base, own), tuple(dict.fromkeys(base_sources + sources))

    # -- public resolution ----------------------------------------------------------

    def select_profile_name(self, name: str | None) -> str:
        """Return the profile to use, falling back to ``default_profile``."""
        chosen = name or self.default_profile
        if not chosen:
            raise ConfigError(
                "No execution profile selected. Pass -c/--profile NAME, or set `default_profile` in "
                f"{user_config_path()} (create it with `qligfep_config init`). "
                f"Available profiles: {', '.join(self.profile_names())}."
            )
        return chosen

    def get_profile(
        self,
        name: str | None = None,
        system: str | None = None,
        q_bin: str | os.PathLike | None = None,
    ) -> Profile:
        """Resolve a profile, applying its per-system overrides when ``system`` is given."""
        chosen = self.select_profile_name(name)
        raw, sources = self._inherited("profiles", chosen)
        raw = dict(raw)
        systems = raw.pop("systems", None) or {}
        if system and systems.get(system):
            raw = deep_merge(raw, systems[system])
        raw["systems"] = {}
        return build_profile(chosen, raw, sources, q_bin_override=q_bin)

    def get_protocol(self, name: str | None = None) -> Protocol:
        chosen = name or self.default_protocol
        raw, sources = self._inherited("protocols", chosen)
        own, _ = self._definition("protocols", chosen, skip_packaged=True)
        return build_protocol(chosen, raw, own, sources)


def load_config(
    project_file: str | os.PathLike | None = None,
    cwd: str | os.PathLike | None = None,
    include_user: bool = True,
) -> Config:
    """Load the packaged defaults, the user file and the project file.

    Args:
        project_file: Explicit project file. Defaults to ``qligfep.yaml`` in ``cwd``.
        cwd: Directory searched for the project file. Defaults to the working directory.
        include_user: Read the per-user file. Disable for hermetic tests.
    """
    layers = [(PACKAGED_LABEL, _read_yaml(PACKAGED_DEFAULTS))]
    _validate_layer(layers[0][1], PACKAGED_LABEL, allow_protocols=True)

    user_file = user_config_path()
    if include_user and user_file.is_file():
        data = _read_yaml(user_file)
        _validate_layer(data, str(user_file), allow_protocols=False)
        layers.append((str(user_file), data))
    elif include_user and os.environ.get("QLIGFEP_CONFIG"):
        raise ConfigError(
            f"$QLIGFEP_CONFIG points to {user_file}, which does not exist. "
            "Create it with `qligfep_config init`, or unset the variable."
        )

    if project_file is not None:
        project = Path(project_file).expanduser()
        if not project.is_file():
            raise ConfigError(f"Config file {project} does not exist.")
    else:
        project = find_project_config(cwd)
    if project is not None:
        data = _read_yaml(project)
        _validate_layer(data, str(project), allow_protocols=True)
        layers.append((str(project), data))
    return Config(layers)
