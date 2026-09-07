"""Download only the pinned reproduction case; never write to Snellius."""
import json
import subprocess

from reproduce import ROOT, verify_download


def fetch():
    selection = json.loads((ROOT / "selection.json").read_text())
    raw = ROOT / "raw"
    if raw.exists():
        raise FileExistsError("raw/ already exists; verify/preserve it rather than overwrite it")
    raw.mkdir()
    host = selection["host"]
    transport = "ssh -T -o BatchMode=yes -o ConnectTimeout=10 -o RemoteCommand=none"

    def copy(remote, destination, patterns):
        destination.mkdir(parents=True)
        command = ["rsync", "-a", "-e", transport]
        command += [f"--include={pattern}" for pattern in patterns]
        command += ["--exclude=*", f"{host}:{remote}/", str(destination)+"/"]
        subprocess.run(command, check=True)

    for item in selection["selection"]:
        copy(item["remote_input"], raw / item["label"] / "inputfiles", ["*.inp", "*.fep", "*.json"])
        copy(item["remote_energy"], raw / item["label"] / "energy",
             ["*.en", "qfep.out", "md_0000_1000.log", "eq5.re", "dualtop.top"])
    copy(selection["campaign"]+"/"+selection["reference_engine_subdir"], raw / "engine",
         ["*.f90", "makefile", "qfep", "qfep.md5"])
    (raw / "metadata").mkdir()
    paths = ["campaign_manifest.json", "qfep_safebar_validation.json", "scripts/analyze_cmet_results.py"]
    subprocess.run(["rsync", "-a", "-e", transport,
                    *(f"{host}:{selection['campaign']}/{path}" for path in paths),
                    str(raw / "metadata")+"/"], check=True)
    print(json.dumps(verify_download(selection), indent=2))


if __name__ == "__main__":
    fetch()
