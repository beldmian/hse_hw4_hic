from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
from pathlib import Path
import shutil
import subprocess

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "data"
ACCESSIONS = ("4DNFIJ2JKO7D", "4DNFI8I9LN74")
CHUNK = 8 * 1024 * 1024


def digest(path):
    with path.open("rb") as f:
        return hashlib.file_digest(f, "md5").hexdigest()


def download(accession):
    metadata = json.loads((DATA / "sources.json").read_text())[accession]
    size, md5 = metadata["file_size"], metadata["md5sum"]
    dest = DATA / f"{accession}.mcool"
    if dest.exists() and dest.stat().st_size == size and digest(dest) == md5:
        print(f"{accession}: verified cached file", flush=True)
        return
    chunks = DATA / ".downloads" / accession
    chunks.mkdir(parents=True, exist_ok=True)

    def fetch(start):
        end = min(start + CHUNK, size) - 1
        part = chunks / str(start)
        if part.exists() and part.stat().st_size == end - start + 1:
            return part
        headers = part.with_suffix(".headers")
        subprocess.run(["curl", "--silent", "--show-error", "--location", "--fail",
            "--retry", "4", "--retry-all-errors", "--connect-timeout", "30", "--max-time", "1200",
            "--range", f"{start}-{end}", "--dump-header", str(headers),
            "--output", str(part.with_suffix(".tmp")), metadata["open_data_url"].replace(".s3.amazonaws.com", ".s3.us-east-1.amazonaws.com") + f"?chunk={start}"], check=True)
        if f"content-range: bytes {start}-{end}/{size}" not in headers.read_text().lower():
            raise ValueError("Server did not honor the requested byte range")
        if part.with_suffix(".tmp").stat().st_size != end - start + 1:
            raise ValueError("Truncated download")
        part.with_suffix(".tmp").replace(part)
        print(f"{accession}: chunk {start // CHUNK + 1}/{(size + CHUNK - 1) // CHUNK}", flush=True)
        return part

    with ThreadPoolExecutor(max_workers=8) as pool:
        parts = list(pool.map(fetch, range(0, size, CHUNK)))
    temporary = dest.with_suffix(".mcool.part")
    with temporary.open("wb") as out:
        for part in parts:
            with part.open("rb") as inp:
                shutil.copyfileobj(inp, out)
    if temporary.stat().st_size != size or digest(temporary) != md5:
        raise ValueError(f"{accession}: size/MD5 mismatch; do not use this file")
    temporary.replace(dest)
    shutil.rmtree(chunks)
    print(f"{accession}: downloaded and MD5 verified", flush=True)


if __name__ == "__main__":
    for accession in ACCESSIONS:
        download(accession)
