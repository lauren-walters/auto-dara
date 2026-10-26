"""Local copy of the COD structures dara can select, so searches never download from COD.

dara's COD index lists about 62,000 entries (about 1.5 GB); those are the only COD structures a search can
use. This copies exactly those with rsync, which is how COD asks bulk users to fetch files.
"""

import shutil
import subprocess
from pathlib import Path

RSYNC_URL = "rsync://www.crystallography.net/cif/"


def wanted_ids() -> list:
    from dara.structure_db import CODDatabase

    ids = {entry[1] for entries in CODDatabase().preparsed_info.values() for entry in entries}
    return sorted(i for i in ids if i.isdigit() and len(i) == 7)


def relative_path(cod_id: str) -> str:
    """COD's own layout, the one dara reads: 1/00/00/1000000.cif."""
    return f"{cod_id[0]}/{cod_id[1:3]}/{cod_id[3:5]}/{cod_id}.cif"


def mirror(dest: Path, bwlimit: int = 4000, ids=None) -> int:
    """Copies the COD files dara can use into `dest`; returns the number still missing afterwards."""
    if shutil.which("rsync") is None:
        raise RuntimeError("rsync is not installed (on Windows, run this inside WSL)")
    ids = wanted_ids() if ids is None else ids
    dest = Path(dest).expanduser()
    dest.mkdir(parents=True, exist_ok=True)
    listing = dest / ".cod_files.txt"
    listing.write_text("".join(f"{relative_path(i)}\n" for i in ids))
    print(f"Copying {len(ids)} COD structures into {dest} (about 1.5 GB for the full set)")
    # rsync exits non-zero when some listed files are not on COD; the count below is what matters
    rsync = subprocess.run(
        ["rsync", "-a", f"--bwlimit={bwlimit}", f"--files-from={listing}", RSYNC_URL, f"{dest}/"],
        check=False,
    )
    missing = [i for i in ids if not (dest / relative_path(i)).is_file()]
    if len(missing) == len(ids):
        raise RuntimeError(
            f"rsync copied nothing (exit {rsync.returncode}); check that this machine can reach "
            "www.crystallography.net on port 873, which some networks block"
        )
    print(f"{len(ids) - len(missing)}/{len(ids)} present in {dest}")
    if missing:
        print(f"{len(missing)} not copied (removed from COD, or the transfer stopped early; rerun to retry)")
    return len(missing)
