"""Settings shared by the command line and the Python API."""

import os
from dataclasses import dataclass, field
from pathlib import Path

DEFAULT_ICSD = Path(os.environ["ICSD_PATH"]).expanduser() if os.environ.get("ICSD_PATH") else None
# The search keeps dara's defaults: looser values let a main phase absorb minor ones. See `auto-dara params`.
DEFAULT_SEARCH_PHASE_PARAMS = {
    "lattice_range": 0.01,
    "k1": "0_0^0.01",
    "k2": "fixed",
    "b1": "0_0^0.005",
    "gewicht": "0_0",
    "rp": 4,
}
DEFAULT_SEARCH_REFINEMENT_PARAMS = {"eps1": 0, "eps2": "0_-0.05^0.05"}
DEFAULT_PHASE_PARAMS = {
    "lattice_range": 0.03,
    "k1": "0_0^0.005",
    "k2": "0_0^0.0001",
    "b1": "0_0^0.005",
    "gewicht": "SPHAR0",
    "rp": 4,
}
DEFAULT_REFINEMENT_PARAMS = {"eps1": 0, "eps2": "0_-0.1^0.1"}

COD_CACHE = Path(os.environ.get("AUTO_DARA_COD_CACHE") or Path.home() / ".cache/auto-dara/cod")


@dataclass
class Settings:
    precursors: str = ""
    system: str = ""
    extra_elements: str = ""
    target: str = ""
    temperature: str = ""
    atmosphere: str = "air"
    method: str = "solid-state synthesis in a box furnace"
    characterization: str = "ex-situ powder XRD at room temperature, after cooling in air"

    databases: tuple = ("icsd", "cod")
    icsd_path: Path | None = DEFAULT_ICSD
    custom_cifs: Path | None = None

    max_e_hull: float | None = 0.1
    skip_filter: bool = False
    use_uniqueness: bool = True
    use_llm: bool = True
    llm_model: str = "gemini-3.1-pro-preview"
    llm_keep_threshold: float = 0.3

    wavelength: str = "Cu"
    instrument: str = "Aeris-fds-Pixcel1d-Medipix3"
    max_phases: int = 5
    n_cases: int = 3
    rpb_threshold: float | None = None
    search_phase_params: dict = field(default_factory=lambda: dict(DEFAULT_SEARCH_PHASE_PARAMS))
    search_refinement_params: dict = field(default_factory=lambda: dict(DEFAULT_SEARCH_REFINEMENT_PARAMS))
    phase_params: dict = field(default_factory=lambda: dict(DEFAULT_PHASE_PARAMS))
    refinement_params: dict = field(default_factory=lambda: dict(DEFAULT_REFINEMENT_PARAMS))
    phase_overrides: dict = field(default_factory=dict)
    cpus: int | None = None

    fom_gof_cutoff: float = 12.0
    fom_weight_threshold: float = 0.01
    fom_lattice_catastrophe: float = 0.5
    check_elements: bool = True
    optional_elements: str = ""

    def _precursor_elements(self) -> set:
        from pymatgen.core import Composition

        return {e.symbol for f in self.precursors.split(",") if f.strip() for e in Composition(f.strip()).elements}

    def _system_elements(self, text) -> set:
        return {t.strip().capitalize() for t in text.replace("-", ",").split(",") if t.strip()}

    def elements(self) -> list:
        return sorted(self._precursor_elements() | self._system_elements(self.system + "," + self.extra_elements))

    def chemsys(self) -> str:
        return "-".join(self.elements())

    def required_elements(self) -> set:
        if not self.check_elements:
            return set()
        req = self._precursor_elements() if self.precursors.strip() else self._system_elements(self.system)
        return req - {"C", "H", "N", "O"} - self._system_elements(self.optional_elements)
