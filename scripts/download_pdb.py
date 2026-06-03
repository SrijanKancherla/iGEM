# downloads HBsAg structure (9UBQ)

import requests
import os

os.makedirs("data/raw", exist_ok=True)

url = "https://files.rcsb.org/download/9UBQ.pdb"
out = "data/raw/9UBQ.pdb"

r = requests.get(url)

with open(out, "w") as f:
    f.write(r.text)

print("Downloaded 9UBQ.pdb")