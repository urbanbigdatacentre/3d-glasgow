#!/usr/bin/env python3
"""
fetch_dtm.py: download the UBDC Glasgow 0.5 m DTM tiles (Zenodo record 13273124, DTM_5x5km.zip, 5.35 GB)
without fetching the whole archive at once: each zip member is read with one HTTP range request and
inflated straight to <out>/<tile>.tif. Small tiles first, several members in parallel, resumable
(existing complete files are skipped).

usage: python3 fetch_dtm.py --out ../DTM_5x5km [--workers 3] [--only NS56SE,NS56NE]
"""
import argparse
import io
import os
import struct
import sys
import time
import urllib.request
import zipfile
import zlib
from concurrent.futures import ThreadPoolExecutor, as_completed

URL = "https://zenodo.org/api/records/13273124/files/DTM_5x5km.zip/content"


class RemoteFile(io.RawIOBase):
    """Seekable read-only view of an HTTP resource (used only to read the zip directory)."""
    def __init__(self, url):
        r = urllib.request.urlopen(urllib.request.Request(url, method="HEAD"), timeout=60)
        self.size = int(r.headers["Content-Length"])
        self.url = r.geturl()
        self.pos = 0
    def seekable(self): return True
    def readable(self): return True
    def tell(self): return self.pos
    def seek(self, off, whence=0):
        self.pos = {0: off, 1: self.pos + off, 2: self.size + off}[whence]
        return self.pos
    def read(self, n=-1):
        if n < 0:
            n = self.size - self.pos
        if n == 0:
            return b""
        req = urllib.request.Request(self.url, headers={"Range": "bytes=%d-%d" % (self.pos, self.pos + n - 1)})
        data = urllib.request.urlopen(req, timeout=120).read()
        self.pos += len(data)
        return data


def fetch_member(url, info, dest):
    """Stream one deflated member to dest (inflating on the fly), verifying the CRC."""
    # local file header: 30 bytes + name + extra (lengths come from the local header itself)
    req = urllib.request.Request(url, headers={"Range": "bytes=%d-%d" % (info.header_offset, info.header_offset + 29)})
    hdr = urllib.request.urlopen(req, timeout=120).read()
    sig, _, _, method, _, _, _, csize, usize, nlen, xlen = struct.unpack("<IHHHHHIIIHH", hdr)
    assert sig == 0x04034b50, "bad local header for %s" % info.filename
    start = info.header_offset + 30 + nlen + xlen
    end = start + info.compress_size - 1
    tmp = dest + ".part"
    for attempt in range(5):
        try:
            req = urllib.request.Request(url, headers={"Range": "bytes=%d-%d" % (start, end)})
            resp = urllib.request.urlopen(req, timeout=300)
            d = zlib.decompressobj(-15) if info.compress_type == zipfile.ZIP_DEFLATED else None
            crc = 0
            t0 = time.time()
            got = 0
            with open(tmp, "wb") as f:
                while True:
                    chunk = resp.read(1 << 20)
                    if not chunk:
                        break
                    got += len(chunk)
                    out = d.decompress(chunk) if d else chunk
                    crc = zlib.crc32(out, crc)
                    f.write(out)
                if d:
                    out = d.flush()
                    crc = zlib.crc32(out, crc)
                    f.write(out)
            if got != info.compress_size or (crc & 0xFFFFFFFF) != info.CRC:
                raise IOError("size/crc mismatch (%d of %d bytes)" % (got, info.compress_size))
            os.replace(tmp, dest)
            return "%s %.0f MB in %.0fs" % (os.path.basename(dest), info.file_size / 1e6, time.time() - t0)
        except Exception as e:
            if attempt == 4:
                raise
            time.sleep(5 * (attempt + 1))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--workers", type=int, default=3)
    ap.add_argument("--only", default=None, help="comma-separated tile names")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)
    rf = RemoteFile(URL)
    z = zipfile.ZipFile(rf)
    members = [i for i in z.infolist() if i.filename.lower().endswith(".tif")]
    if args.only:
        keep = set(args.only.split(","))
        members = [m for m in members if os.path.basename(m.filename).split("_")[0] in keep]
    members.sort(key=lambda m: m.compress_size)          # small tiles first
    todo = []
    for m in members:
        dest = os.path.join(args.out, os.path.basename(m.filename))
        if os.path.exists(dest) and os.path.getsize(dest) == m.file_size:
            print("skip (done):", os.path.basename(dest), flush=True)
        else:
            todo.append((m, dest))
    print("downloading %d tiles, %.2f GB" % (len(todo), sum(m.compress_size for m, _ in todo) / 1e9), flush=True)
    with ThreadPoolExecutor(max_workers=args.workers) as ex:
        futs = {ex.submit(fetch_member, rf.url, m, dest): m for m, dest in todo}
        for fut in as_completed(futs):
            try:
                print("done:", fut.result(), flush=True)
            except Exception as e:
                print("FAILED:", futs[fut].filename, e, flush=True)
                sys.exit(1)
    print("all done", flush=True)


if __name__ == "__main__":
    main()
