"""
Preuzimanje podataka sa Zenodo repozitorijuma.

Dataset: https://zenodo.org/records/15866724
Preuzima sve .h5ad, .csv i .rds fajlove u folder data/raw/.
"""

import json
from pathlib import Path

import requests
import yaml


def load_config(path: str):
    """
    Učitava config.yaml da dobije Zenodo record ID i putanju za raw podatke.
    """
    with open(path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def download_file(url: str, output_path: Path):
    """
    Preuzima jedan fajl sa URL-a i čuva ga lokalno.

    Ako fajl već postoji, preskače download (korisno za ponovno pokretanje).
    Stream download (chunk po chunk) štedi memoriju za velike .h5ad fajlove (~200 MB).
    """
    output_path.parent.mkdir(parents=True, exist_ok=True)
    if output_path.exists():
        print(f"[SKIP] {output_path.name} already exists")
        return

    print(f"[DOWNLOADING] {output_path.name}")
    with requests.get(url, stream=True, timeout=120) as r:
        r.raise_for_status()
        with open(output_path, "wb") as f:
            for chunk in r.iter_content(chunk_size=8192):
                if chunk:
                    f.write(chunk)
    print(f"[DONE] {output_path}")


def main():
    """
    Glavna funkcija za download.

    Koraci:
      1. Poziva Zenodo API i dobija listu fajlova
      2. Čuva metadata u zenodo_record.json (za dokumentaciju)
      3. Preuzima samo relevantne ekstenzije: .h5ad, .csv, .rds
    """
    cfg = load_config("config/config.yaml")
    record_id = cfg["data"]["zenodo_record_id"]
    raw_dir = Path(cfg["data"]["raw_dir"])

    api_url = f"https://zenodo.org/api/records/{record_id}"
    response = requests.get(api_url, timeout=60)
    response.raise_for_status()
    payload = response.json()

    with open(raw_dir / "zenodo_record.json", "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2)

    allowed_ext = (".h5ad", ".csv", ".rds", ".annoy")
    for file_obj in payload.get("files", []):
        key = file_obj["key"]
        if key.lower().endswith(allowed_ext):
            url = file_obj["links"]["self"]
            download_file(url, raw_dir / key)

    azimuth_cfg = cfg.get("azimuth", {})
    azimuth_record = azimuth_cfg.get("record_id")
    if azimuth_record:
        print(f"Downloading Azimuth reference record {azimuth_record}...")
        api_url = f"https://zenodo.org/api/records/{azimuth_record}"
        response = requests.get(api_url, timeout=60)
        response.raise_for_status()
        payload = response.json()
        with open(raw_dir / f"zenodo_record_{azimuth_record}.json", "w", encoding="utf-8") as f:
            json.dump(payload, f, indent=2)
        for file_obj in payload.get("files", []):
            key = file_obj["key"]
            if key.lower().endswith(allowed_ext):
                url = file_obj["links"]["self"]
                download_file(url, raw_dir / key)

    print("Download complete.")


if __name__ == "__main__":
    main()
