import json
import struct
import sqlite3
import subprocess

batch_size = 10000

sbs = subprocess.Popen(
    [
        "sbs-observation",
        "add",
        "--batch-size",
        str(batch_size),
        "-D",
        "ztf.db",
        "-i",
        "-",
    ],
    stdin=subprocess.PIPE,
)

db = sqlite3.connect("/home/msk/data/ztf/zdata/zchecker.db")
db.row_factory = sqlite3.Row

results = db.execute(
    "SELECT pid,jd_start - 2400000.5 AS mjd_start,jd_stop - 2400000.5 AS mjd_stop, fov "
    "FROM ztf INNER JOIN obs USING (obsid);"
)

observations = []
while True:
    rows = results.fetchmany(3000)
    if len(rows) == 0:
        break

    for row in rows:
        corners = [
            a * 57.29577951308232087685 for a in struct.unpack("10d", row["fov"])
        ]  # degrees
        fov = []
        for i in range(4):
            fov.append(f"{corners[2 * i]}:{corners[2 * i + 1]}")
        fov = ",".join(fov)

        observations.append(
            {
                "source": "ztf",
                "observatory": "I41",
                "product_id": str(row["pid"]),
                "mjd_start": row["mjd_start"],
                "mjd_stop": row["mjd_stop"],
                "fov": fov,
            }
        )

    if len(observations) >= batch_size:
        sbs.stdin.write(json.dumps(observations).encode())
        observations = []

# left overs
if len(observations) > 0:
    sbs.stdin.write(json.dumps(observations).encode())

sbs.stdin.close()
sbs.wait()
